library(simcausal)
library(tibble)
library(dplyr)
library(purrr)
library(mgcv)
library(readr)
library(tidyr)
library(ggplot2)
simLIV <- function(n,D){
  mydat <- tibble(sim(D,n=n,wide = F)) %>%
    select(-contains("OPA"), 
           -contains("deltaA"), 
           -contains("pA0"),
           -contains("U"))
  return(mydat)  
}

# PIV-style nuisance formulas. smoothness controls k = smoothness^length(C),
# the basis dimension; mgcv still estimates the smoothing penalties itself.
auto_smooth_term <- function(vars, data, discrete_threshold = 5L,
                             smoothness = 3L, bs_type = "tp") {
  for (name in c("smoothness", "discrete_threshold")) {
    value <- get(name)
    lower <- if (name == "smoothness") 3L else 1L
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value != floor(value) || value < lower ||
        value > .Machine$integer.max) {
      stop(sprintf("`%s` must be an integer >= %d.", name, lower), call. = FALSE)
    }
  }
  if (length(bs_type) != 1L || !is.character(bs_type) || is.na(bs_type) ||
      !nzchar(bs_type)) stop("`bs_type` must be a non-empty string.", call. = FALSE)
  if (!is.data.frame(data) || !is.character(vars) || anyNA(vars) ||
      any(!vars %in% names(data))) {
    stop("`vars` must name columns in `data`.", call. = FALSE)
  }
  vars <- unique(vars)
  observed <- function(x) x[!is.na(x)]
  # Constant columns cannot define contrasts or smooths; the intercept absorbs them.
  vars <- vars[vapply(vars, function(v) length(unique(observed(data[[v]]))) > 1L,
                      logical(1L))]
  is_discrete <- function(x) {
    is.factor(x) || is.character(x) ||
      length(unique(observed(x))) <= discrete_threshold
  }
  D <- vars[vapply(vars, function(v) is_discrete(data[[v]]), logical(1L))]
  C <- setdiff(vars, D)
  if (!length(D) && !length(C)) return("1")
  
  quote_vars <- function(x) vapply(x, function(v) {
    paste(deparse(as.name(v), backtick = TRUE), collapse = "")
  }, character(1L))
  discrete_term <- if (length(D)) {
    sprintf("interaction(%s, drop = TRUE)", paste(quote_vars(D), collapse = ", "))
  } else NULL
  if (!length(C)) return(discrete_term)
  
  k <- smoothness ^ length(C)
  if (!is.finite(k) || k > .Machine$integer.max) {
    stop("Requested smooth basis is too large; reduce `smoothness`.", call. = FALSE)
  }
  continuous_vars <- paste(quote_vars(C), collapse = ", ")
  basis <- encodeString(bs_type, quote = '"')
  if (!length(D)) {
    return(sprintf("s(%s, k = %d, bs = %s)", continuous_vars, as.integer(k), basis))
  }
  c(sprintf("s(%s, by = %s, k = %d, bs = %s)",
            continuous_vars, discrete_term, as.integer(k), basis), discrete_term)
}

smooth_all <- function(response, data, smoothness = 3L,
                       discrete_threshold = 5L, bs_type = "tp") {
  if (length(response) != 1L || !is.character(response) || is.na(response) ||
      !nzchar(response)) stop("`response` must be a non-empty string.", call. = FALSE)
  terms <- auto_smooth_term(setdiff(names(data), response), data,
                            discrete_threshold = discrete_threshold,
                            smoothness = smoothness, bs_type = bs_type)
  lhs <- paste(deparse(as.name(response), backtick = TRUE), collapse = "")
  stats::as.formula(paste(lhs, "~", paste(terms, collapse = " + ")))
}

# K=2; H=1; mydat=data_LIV; action_list=c(0,0); time=TM;k=2
# Debiased machine learning with longituidinal IV (prespecify weighting function)
# K subject folds share one split across all time points; K=1 uses all subjects.
# workers=K runs folds in parallel; workers=1 runs the same folds serially.
# Call set.seed() before fitting to reproduce the split and fold-specific RNGs.
# smoothness sets k=smoothness^number_of_continuous_variables; discrete_threshold
# classifies numeric predictors by their distinct training-fold values.
DML_LIV <- function(mydat,action_list,K=2, H=2,print_state=F, workers=K,
                    smoothness = 3L, discrete_threshold = 5L){
  # Validate smoothing controls before creating folds or starting workers.
  auto_smooth_term(character(), data.frame(), discrete_threshold, smoothness)
  required_columns <- c("ID", "t", "L", "Z", "A", "Y")
  if (!is.data.frame(mydat) || !all(required_columns %in% names(mydat))) {
    stop("`mydat` must contain ID, t, L, Z, A, and Y columns.", call. = FALSE)
  }
  mydat <- tibble::as_tibble(mydat)
  if (is.factor(mydat$ID)) mydat$ID <- as.character(mydat$ID)
  if ((!is.numeric(mydat$ID) && !is.character(mydat$ID)) ||
      anyNA(mydat$ID) || !nrow(mydat) ||
      (is.numeric(mydat$ID) && any(!is.finite(mydat$ID)))) {
    stop("`ID` must contain non-missing numeric or character identifiers.", call. = FALSE)
  }
  if (!is.numeric(mydat$t) || anyNA(mydat$t) ||
      any(!is.finite(mydat$t)) || any(mydat$t < 0 | mydat$t != floor(mydat$t))) {
    stop("`t` must contain nonnegative integer time points.", call. = FALSE)
  }
  time_points <- sort(unique(mydat$t))
  if (length(time_points) < 2L || time_points[1L] != 0 ||
      any(diff(time_points) != 1)) {
    stop("Time points must run consecutively from 0 through the outcome time T+1.",
         call. = FALSE)
  }
  if (anyDuplicated(mydat[c("ID", "t")])) {
    stop("Each (ID, t) pair must occur exactly once.", call. = FALSE)
  }
  id_order <- sort(unique(mydat$ID))
  n <- length(id_order)
  if (length(K) != 1L || !is.numeric(K) || is.na(K) || !is.finite(K) ||
      K != floor(K) || K < 1L || K > n) {
    stop("`K` must be an integer between 1 and the number of individuals.",
         call. = FALSE)
  }
  if (length(workers) != 1L || !is.numeric(workers) || is.na(workers) ||
      !is.finite(workers) || workers < 1L || workers != floor(workers)) {
    stop("`workers` must be a positive integer.", call. = FALSE)
  }
  workers <- min(workers, K)
  
  # Align every time slice by ID before combining histories or initializing ETA.
  # ETA row i and every history row i always refer to id_order[i].
  aligned_data <- lapply(time_points, function(time_point) {
    current <- mydat[mydat$t == time_point, , drop = FALSE]
    index <- match(id_order, current$ID)
    if (nrow(current) != n || anyNA(index)) {
      stop(sprintf("Every individual must have one record at time %d.", time_point),
           call. = FALSE)
    }
    current[index, , drop = FALSE]
  })
  mydat <- dplyr::bind_rows(aligned_data)
  TM <- max(time_points) - 1L
  
  make_id_folds <- function(ids, K) {
    # K=1 fits and evaluates all individuals without cross-fitting.
    if (K == 1L) {
      index <- seq_along(ids)
      return(list(list(train_ids = ids, test_ids = ids,
                       train_idx = index, test_idx = index)))
    }
    assignment <- sample(rep(seq_len(K), length.out = length(ids)))
    lapply(seq_len(K), function(k) {
      train_idx <- which(assignment != k)
      test_idx <- which(assignment == k)
      list(train_ids = ids[train_idx], test_ids = ids[test_idx],
           train_idx = train_idx, test_idx = test_idx)
    })
  }
  
  build_fY_history <- function(mydat, time, H){
    fY <- mydat %>% filter(t == time) %>% select(L, Z, A)
    
    h<-1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L, Z, A) %>%
          rename_with(~ paste0(., "_lag", h))
        fY <- bind_cols(fY, temp)
      } else {
        break
      }
      h <- h+1
    }  
    return(fY)
  }
  
  build_fA_history <- function(mydat, time, H) {
    fA <- mydat[which(mydat$t==time),] %>% select(L,Z)
    
    h <- 1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L, Z) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    h <- 1
    while(h <=H) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(A) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    return(fA)
  }
  
  build_fZ_history <- function(mydat, time, H) {
    fA <- mydat[which(mydat$t==time),] %>% select(L)
    
    h <- 1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    h <- 1
    while(h <=H) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(Z, A) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    return(fA)
  }
  
  
  
  # Split once by individual; these same indices are used at every time point.
  folds <- make_id_folds(id_order, K = K)
  fold_seeds <- sample.int(.Machine$integer.max, K)
  master_rng <- get(".Random.seed", envir = globalenv())
  # Serial and parallel execution must leave the caller at the same RNG state.
  on.exit(assign(".Random.seed", master_rng, envir = globalenv()), add = TRUE)
  state <- list(
    n = n, TM = TM, id_order = id_order, folds = folds, fold_seeds = fold_seeds,
    rng_kind = RNGkind(),
    histories = lapply(seq.int(0L, TM), function(time) build_fZ_history(mydat, time, H)),
    instruments = lapply(aligned_data[seq_len(TM + 1L)], function(x) x$Z),
    treatments = lapply(aligned_data[seq_len(TM + 1L)], function(x) x$A),
    outcome = aligned_data[[TM + 2L]]$Y,
    action_list = action_list, print_state = print_state,
    smoothness = smoothness, discrete_threshold = discrete_threshold
  )
  
  # One job owns a complete backward recursion for a single, fixed fold.
  fit_fold <- function(k, state) {
    do.call(RNGkind, as.list(state$rng_kind))
    set.seed(state$fold_seeds[k])
    train <- state$folds[[k]]$train_idx
    ETA_fold <- matrix(0, nrow = state$n, ncol = state$TM + 2L,
                       dimnames = list(as.character(state$id_order), NULL))
    ETA_fold[, state$TM + 2L] <- state$outcome
    kappas <- if (state$print_state) vector("list", state$TM + 1L) else NULL
    time <- state$TM
    while (time >= 0L) {
      oY <- ETA_fold[, time + 2L]
      oZ <- state$instruments[[time + 1L]]
      his <- state$histories[[time + 1L]]
      trt <- tryCatch(
        {state$action_list[[time + 1]](his)},
        error = function(e){state$action_list[time + 1]}
      )
      
      oA <- state$treatments[[time + 1L]]
      oA_trt <- (oA==trt)
      
      his_A <- his %>% mutate(oA=oA)
      
      
      
      famil <- if (length(unique(oZ[train])) <= 2) {
        binomial()
      } else {
        gaussian()
      }
      
      fom <- smooth_all("o", his[train,], smoothness = state$smoothness,
                        discrete_threshold = state$discrete_threshold)
      fom_A <- smooth_all("o", his_A[train,], smoothness = state$smoothness,
                          discrete_threshold = state$discrete_threshold)
      
      fit_rho   <- gam(fom, data = his[train,] %>% mutate(o = oZ[train]),family = famil) # Z~L
      rho <- as.vector(predict(fit_rho, newdata = his,type='response'))
      
      fit_delta <- gam(fom, data = his[train,] %>% mutate(o = oA[train]),family = binomial()) # A~L
      delta1 <- as.vector(predict(fit_delta, newdata = his,type='response'))
      delta <- trt * delta1 + (1-trt) * (1-delta1)
      
      fit_eta   <- gam(fom_A, data = his_A[train,] %>% mutate(o = oY[train]),family = gaussian()) # Y~L,A
      eta1 <- as.vector(predict(fit_eta, newdata = his %>% mutate(oA=1),type='response')) * delta1
      eta0 <- as.vector(predict(fit_eta, newdata = his %>% mutate(oA=0),type='response')) * (1 - delta1)
      eta <- trt * eta1 + (1-trt) * eta0
      
      pseudo_kappa <- (oA_trt - delta) * (oZ - rho)
      fit_kappa <- gam(fom, data = his[train,] %>% mutate(o = pseudo_kappa[train]),family = gaussian())
      kappa <- as.vector(predict(fit_kappa, newdata = his,type='response'))
      
      
      pseudo_zeta <- (oY * oA_trt - eta) * (oZ - rho)
      fit_zeta  <- gam(fom, data = his[train,] %>% mutate(o = pseudo_zeta[train]),family = gaussian())
      zeta <- as.vector(predict(fit_zeta, newdata = his,type='response'))
      gamma <- zeta/kappa
      
      
      if (state$print_state) kappas[[time + 1L]] <- kappa
      est_aug <- pseudo_zeta/kappa + (1 - pseudo_kappa/kappa) * gamma
      ETA_fold[, time + 1L] <- est_aug
      time <- time - 1L
    }
    list(ETA = ETA_fold, kappas = kappas)
  }
  # Send only the worker's dependencies, not the caller's execution environment.
  worker_smooth_all <- smooth_all
  worker_auto_smooth_term <- auto_smooth_term
  environment(worker_auto_smooth_term) <- baseenv()
  environment(worker_smooth_all) <- list2env(
    list(auto_smooth_term = worker_auto_smooth_term), parent = globalenv()
  )
  environment(fit_fold) <- list2env(
    list(smooth_all = worker_smooth_all, gam = gam, predict = predict),
    parent = globalenv()
  )
  run_fold <- function(k, fit, state) {
    warnings <- list()
    result <- tryCatch(withCallingHandlers(fit(k, state), warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }), error = function(e) {
      stop(sprintf("Fold %d failed: %s", k, conditionMessage(e)), call. = FALSE)
    })
    result$warnings <- warnings
    result
  }
  environment(run_fold) <- baseenv()
  
  if (workers == 1L || K == 1L) {
    fold_results <- lapply(seq_len(K), run_fold, fit = fit_fold, state = state)
  } else {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    initialize_worker <- function(lib_paths) {
      .libPaths(lib_paths)
      suppressPackageStartupMessages({library(dplyr); library(mgcv)})
      NULL
    }
    environment(initialize_worker) <- baseenv()
    parallel::clusterCall(cluster, initialize_worker, .libPaths())
    
    # Retain global constants/helpers referenced by user-supplied treatment rules.
    globals <- list()
    seen <- list()
    collect_globals <- function(fun) {
      env <- environment(fun)
      if (is.null(env) || isNamespace(env) || identical(env, baseenv()) ||
          startsWith(environmentName(env), "package:")) return(invisible(NULL))
      if (any(vapply(seen, function(x) identical(x, fun), logical(1)))) return(invisible(NULL))
      seen[[length(seen) + 1L]] <<- fun
      for (name in codetools::findGlobals(fun, merge = TRUE)) {
        if (exists(name, envir = globalenv(), inherits = FALSE)) {
          globals[name] <<- list(get(name, envir = globalenv(), inherits = FALSE))
        }
        if (exists(name, envir = env, inherits = TRUE)) {
          value <- get(name, envir = env, inherits = TRUE)
          if (is.function(value)) collect_globals(value)
        }
      }
      invisible(NULL)
    }
    for (rule in as.list(action_list)) if (is.function(rule)) collect_globals(rule)
    if (length(globals)) {
      parallel::clusterExport(cluster, names(globals),
                              envir = list2env(globals, parent = emptyenv()))
    }
    fold_results <- parallel::parLapply(cluster, seq_len(K), run_fold,
                                        fit = fit_fold, state = state)
  }
  ETA <- array(0, dim = c(n, K, TM + 2L),
               dimnames = list(as.character(id_order), NULL, NULL))
  for (k in seq_len(K)) {
    ETA[, k, ] <- fold_results[[k]]$ETA
    for (w in fold_results[[k]]$warnings) warning(w)
    if (print_state) {
      for (time in seq.int(TM, 0L)) hist(fold_results[[k]]$kappas[[time + 1L]], breaks = 30)
    }
  }
  return(list(
    ETA = ETA,
    folds = folds,
    id_order = id_order,
    fold_seeds = fold_seeds,
    workers = workers,
    input_data = mydat,
    Action = action_list, 
    H = H,
    smoothness = smoothness, discrete_threshold = discrete_threshold
  ))
}













# K=2; H=1; mydat=data_LIV; action_list=c(0,0); time=TM;k=2
# Debiased machine learning with longituidinal IV (adaptive weighting function)
# K subject folds share one split across all time points; K=1 uses all subjects.
# workers=K runs folds in parallel; workers=1 runs the same folds serially.
# Call set.seed() before fitting to reproduce the split and fold-specific RNGs.
# smoothness sets k=smoothness^number_of_continuous_variables; discrete_threshold
# classifies numeric predictors by their distinct training-fold values.
DML_LIV_adap <- function(mydat,action_list,K=2, H=2,print_state=F, workers=K,
                         smoothness = 3L, discrete_threshold=5L){
  # Validate smoothing controls before creating folds or starting workers.
  auto_smooth_term(character(), data.frame(), discrete_threshold, smoothness)
  required_columns <- c("ID", "t", "L", "Z", "A", "Y")
  if (!is.data.frame(mydat) || !all(required_columns %in% names(mydat))) {
    stop("`mydat` must contain ID, t, L, Z, A, and Y columns.", call. = FALSE)
  }
  mydat <- tibble::as_tibble(mydat)
  if (is.factor(mydat$ID)) mydat$ID <- as.character(mydat$ID)
  if ((!is.numeric(mydat$ID) && !is.character(mydat$ID)) ||
      anyNA(mydat$ID) || !nrow(mydat) ||
      (is.numeric(mydat$ID) && any(!is.finite(mydat$ID)))) {
    stop("`ID` must contain non-missing numeric or character identifiers.", call. = FALSE)
  }
  if (!is.numeric(mydat$t) || anyNA(mydat$t) ||
      any(!is.finite(mydat$t)) || any(mydat$t < 0 | mydat$t != floor(mydat$t))) {
    stop("`t` must contain nonnegative integer time points.", call. = FALSE)
  }
  time_points <- sort(unique(mydat$t))
  if (length(time_points) < 2L || time_points[1L] != 0 ||
      any(diff(time_points) != 1)) {
    stop("Time points must run consecutively from 0 through the outcome time T+1.",
         call. = FALSE)
  }
  if (anyDuplicated(mydat[c("ID", "t")])) {
    stop("Each (ID, t) pair must occur exactly once.", call. = FALSE)
  }
  id_order <- sort(unique(mydat$ID))
  n <- length(id_order)
  if (length(K) != 1L || !is.numeric(K) || is.na(K) || !is.finite(K) ||
      K != floor(K) || K < 1L || K > n) {
    stop("`K` must be an integer between 1 and the number of individuals.",
         call. = FALSE)
  }
  
  if (length(workers) != 1L || !is.numeric(workers) || is.na(workers) ||
      !is.finite(workers) || workers < 1L || workers != floor(workers)) {
    stop("`workers` must be a positive integer.", call. = FALSE)
  }
  workers <- min(workers, K)
  
  # Align every time slice by ID before combining histories or initializing ETA.
  # ETA row i and every history row i always refer to id_order[i].
  aligned_data <- lapply(time_points, function(time_point) {
    current <- mydat[mydat$t == time_point, , drop = FALSE]
    index <- match(id_order, current$ID)
    if (nrow(current) != n || anyNA(index)) {
      stop(sprintf("Every individual must have one record at time %d.", time_point),
           call. = FALSE)
    }
    current[index, , drop = FALSE]
  })
  mydat <- dplyr::bind_rows(aligned_data)
  TM <- max(time_points) - 1L
  
  make_id_folds <- function(ids, K) {
    # K=1 fits and evaluates all individuals without cross-fitting.
    if (K == 1L) {
      index <- seq_along(ids)
      return(list(list(train_ids = ids, test_ids = ids,
                       train_idx = index, test_idx = index)))
    }
    assignment <- sample(rep(seq_len(K), length.out = length(ids)))
    lapply(seq_len(K), function(k) {
      train_idx <- which(assignment != k)
      test_idx <- which(assignment == k)
      list(train_ids = ids[train_idx], test_ids = ids[test_idx],
           train_idx = train_idx, test_idx = test_idx)
    })
  }
  
  build_fY_history <- function(mydat, time, H){
    fY <- mydat %>% filter(t == time) %>% select(L, Z, A)
    h<-1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L, Z, A) %>%
          rename_with(~ paste0(., "_lag", h))
        fY <- bind_cols(fY, temp)
      } else {
        break
      }
      h <- h+1
    }  
    return(fY)
  }
  
  build_fA_history <- function(mydat, time, H) {
    fA <- mydat[which(mydat$t==time),] %>% select(L,Z)
    h <- 1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L, Z) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    h <- 1
    while(h <=H) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(A) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    return(fA)
  }
  
  build_fZ_history <- function(mydat, time, H) {
    fA <- mydat[which(mydat$t==time),] %>% select(L)
    
    h <- 1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    h <- 1
    while(h <=H) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(Z, A) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    return(fA)
  }
  
  
  
  # Split once by individual; these same indices are used at every time point.
  folds <- make_id_folds(id_order, K = K)
  fold_seeds <- sample.int(.Machine$integer.max, K)
  master_rng <- get(".Random.seed", envir = globalenv())
  # Serial and parallel execution must leave the caller at the same RNG state.
  on.exit(assign(".Random.seed", master_rng, envir = globalenv()), add = TRUE)
  state <- list(
    n = n, TM = TM, id_order = id_order, folds = folds, fold_seeds = fold_seeds,
    rng_kind = RNGkind(),
    histories = lapply(seq.int(0L, TM), function(time) build_fZ_history(mydat, time, H)),
    instruments = lapply(aligned_data[seq_len(TM + 1L)], function(x) x$Z),
    treatments = lapply(aligned_data[seq_len(TM + 1L)], function(x) x$A),
    outcome = aligned_data[[TM + 2L]]$Y,
    action_list = action_list, print_state = print_state,
    smoothness = smoothness, discrete_threshold = discrete_threshold
  )
  
  # One job owns a complete backward recursion for a single, fixed fold.
  fit_fold <- function(k, state) {
    do.call(RNGkind, as.list(state$rng_kind))
    set.seed(state$fold_seeds[k])
    train <- state$folds[[k]]$train_idx
    ETA_fold <- matrix(0, nrow = state$n, ncol = state$TM + 2L,
                       dimnames = list(as.character(state$id_order), NULL))
    ETA_fold[, state$TM + 2L] <- state$outcome
    kappas <- if (state$print_state) vector("list", state$TM + 1L) else NULL
    time <- state$TM
    while (time >= 0L) {
      oY <- ETA_fold[, time + 2L]
      oZ <- state$instruments[[time + 1L]]
      his <- state$histories[[time + 1L]]
      trt <- tryCatch({
        state$action_list[[time + 1]](his)
      }, error = function(e) {
        state$action_list[time + 1]
      })
      oA <- state$treatments[[time + 1L]]
      oA_trt <- (oA==trt)
      
      his_A <- his %>% mutate(oA=oA)
      his_Z <- his %>% mutate(oZ=oZ)
      his_ZA <- his_A %>% mutate(oZ=oZ)
      
      fom <- smooth_all("o", his[train,], smoothness = state$smoothness,
                        discrete_threshold = state$discrete_threshold)
      fom_A <- smooth_all("o", his_A[train,], smoothness = state$smoothness,
                          discrete_threshold = state$discrete_threshold)
      fom_Z <- smooth_all("o", his_Z[train,], smoothness = state$smoothness,
                          discrete_threshold = state$discrete_threshold)
      fom_ZA <- smooth_all("o", his_ZA[train,], smoothness = state$smoothness,
                           discrete_threshold = state$discrete_threshold)
      
      fit_delta <- gam(fom, data = his[train,] %>% mutate(o = oA[train]),family = binomial()) # A ~ L
      delta1 <- as.vector(predict(fit_delta, newdata = his,type='response'))
      delta <- trt * delta1 + (1-trt) * (1-delta1)
      
      fit_pi <- gam(fom_Z, data = his_Z[train,] %>% mutate(o = oA[train]),family = binomial()) # A ~ L,Z
      pi1 <- as.vector(predict(fit_pi, newdata = his_Z,type='response'))
      pi <- trt * pi1 + (1-trt) * (1-pi1)
      
      fit_eta <- gam(fom_A, data = his_A[train,] %>% mutate(o = oY[train]),family = gaussian()) # Y ~ L,A
      eta1  <- as.vector(predict(fit_eta,  newdata = his  %>% mutate(oA = 1), type='response')) * delta1
      eta0  <- as.vector(predict(fit_eta,  newdata = his  %>% mutate(oA = 0), type='response')) * (1- delta1)
      eta   <- trt * eta1 + (1-trt) * eta0
      
      fit_xi  <- gam(fom_ZA, data = his_ZA[train,] %>% mutate(o = oY[train]),family = gaussian()) # Y~L,A,Z
      xi1   <- as.vector(predict(fit_xi, newdata = his_Z %>% mutate(oA = 1), type='response')) * pi1
      xi0   <- as.vector(predict(fit_xi, newdata = his_Z %>% mutate(oA = 0), type='response')) * (1-pi1)
      xi    <- trt * xi1 + (1-trt) * xi0 
      
      pseudo_kappa <- (oA - delta1)**2 - (oA - pi1)**2 
      fit_kappa <- gam(fom, data = his[train,] %>% mutate(o = pseudo_kappa[train]),family = gaussian(link = 'log')) 
      kappa <- as.vector(predict(fit_kappa, newdata = his,type='response'))
      
      pseudo_zeta <- (pi - delta) * oA_trt * oY + (oA_trt - pi) * xi - (oA_trt - delta) * eta
      fit_zeta <- gam(fom, data = his[train,] %>% mutate(o = pseudo_zeta[train]),family = gaussian())  
      zeta <- as.vector(predict(fit_zeta, newdata = his, type='response'))
      gamma <- zeta/kappa
      
      if (state$print_state) kappas[[time + 1L]] <- kappa
      est_aug <- (pseudo_zeta + gamma * (kappa + (oA_trt-pi)^2 - (oA_trt-delta)^2)) / kappa
      ETA_fold[, time + 1L] <- est_aug
      time <- time - 1L
    }
    list(ETA = ETA_fold, kappas = kappas)
  }
  # Send only the worker's dependencies, not the caller's execution environment.
  worker_smooth_all <- smooth_all
  worker_auto_smooth_term <- auto_smooth_term
  environment(worker_auto_smooth_term) <- baseenv()
  environment(worker_smooth_all) <- list2env(
    list(auto_smooth_term = worker_auto_smooth_term), parent = globalenv()
  )
  environment(fit_fold) <- list2env(
    list(smooth_all = worker_smooth_all, gam = gam, predict = predict),
    parent = globalenv()
  )
  run_fold <- function(k, fit, state) {
    warnings <- list()
    result <- tryCatch(withCallingHandlers(fit(k, state), warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }), error = function(e) {
      stop(sprintf("Fold %d failed: %s", k, conditionMessage(e)), call. = FALSE)
    })
    result$warnings <- warnings
    result
  }
  environment(run_fold) <- baseenv()
  
  if (workers == 1L || K == 1L) {
    fold_results <- lapply(seq_len(K), run_fold, fit = fit_fold, state = state)
  } else {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    initialize_worker <- function(lib_paths) {
      .libPaths(lib_paths)
      suppressPackageStartupMessages({library(dplyr); library(mgcv)})
      NULL
    }
    environment(initialize_worker) <- baseenv()
    parallel::clusterCall(cluster, initialize_worker, .libPaths())
    
    # Retain global constants/helpers referenced by user-supplied treatment rules.
    globals <- list()
    seen <- list()
    collect_globals <- function(fun) {
      env <- environment(fun)
      if (is.null(env) || isNamespace(env) || identical(env, baseenv()) ||
          startsWith(environmentName(env), "package:")) return(invisible(NULL))
      if (any(vapply(seen, function(x) identical(x, fun), logical(1)))) return(invisible(NULL))
      seen[[length(seen) + 1L]] <<- fun
      for (name in codetools::findGlobals(fun, merge = TRUE)) {
        if (exists(name, envir = globalenv(), inherits = FALSE)) {
          globals[name] <<- list(get(name, envir = globalenv(), inherits = FALSE))
        }
        if (exists(name, envir = env, inherits = TRUE)) {
          value <- get(name, envir = env, inherits = TRUE)
          if (is.function(value)) collect_globals(value)
        }
      }
      invisible(NULL)
    }
    for (rule in as.list(action_list)) if (is.function(rule)) collect_globals(rule)
    if (length(globals)) {
      parallel::clusterExport(cluster, names(globals),
                              envir = list2env(globals, parent = emptyenv()))
    }
    fold_results <- parallel::parLapply(cluster, seq_len(K), run_fold,
                                        fit = fit_fold, state = state)
  }
  ETA <- array(0, dim = c(n, K, TM + 2L),
               dimnames = list(as.character(id_order), NULL, NULL))
  for (k in seq_len(K)) {
    ETA[, k, ] <- fold_results[[k]]$ETA
    for (w in fold_results[[k]]$warnings) warning(w)
    if (print_state) {
      for (time in seq.int(TM, 0L)) hist(fold_results[[k]]$kappas[[time + 1L]], breaks = 30)
    }
  }
  return(list(
    ETA = ETA,
    folds = folds,
    id_order = id_order,
    fold_seeds = fold_seeds,
    workers = workers,
    input_data = mydat,
    Action = action_list, 
    H = H,
    smoothness = smoothness, discrete_threshold = discrete_threshold
  ))
}












predict_LIV <- function(fit_LIV){
  
  predict_LIV_sub <- function(fit_LIV, time){
    n <- dim(fit_LIV$ETA)[1L]
    id_order <- fit_LIV$id_order
    if (is.null(id_order)) id_order <- seq_len(n)
    if (length(id_order) != n || anyNA(id_order) || anyDuplicated(id_order)) {
      stop("The fit must identify each ETA row with a unique ID.", call. = FALSE)
    }
    K <- length(fit_LIV$folds)
    ETA_merge <- rep(0,n)
    coverage <- integer(n)
    for (k in 1:K) {
      index <- match(fit_LIV$folds[[k]]$test_ids, id_order)
      if (!length(index) || anyNA(index) || anyDuplicated(index)) {
        stop("Test-fold IDs do not match the ETA row identifiers.", call. = FALSE)
      }
      ETA_merge[index] <- fit_LIV$ETA[index,k,time+1]
      coverage[index] <- coverage[index] + 1L
    }
    if (any(coverage != 1L)) {
      stop("Every individual must belong to exactly one test fold.", call. = FALSE)
    }
    ATE <- mean(ETA_merge)
    SD <- sd(ETA_merge)/sqrt(n)
    return(list(time=time,
                est=ATE, 
                std=SD
    ))
  }
  
  res <- unlist(predict_LIV_sub(fit_LIV,time=0))
  for (time in 1:(dim(fit_LIV$ETA)[3]-1)) {
    res <- rbind(res,unlist(predict_LIV_sub(fit_LIV,time)))  
  }
  res <- as.data.frame(res)
  rownames(res) <- NULL
  res <- res %>%
    mutate(
      lower = est - 1.96 * std,
      upper = est + 1.96 * std
    )
  
  return(res)
}





ATE <- function(D,TM,PROB,accuracy=1e5){
  generate_action_from_vector <- function(treatment_vec, varname = "A") {
    nodes=list()
    for (i in 1: length(treatment_vec)) {
      temp=treatment_vec[i]
      if(temp==1){
        nodes = c(nodes, node(varname, t = i-1, distr = "rconst", const = 1))  
      }else{
        nodes = c(nodes, node(varname, t = i-1, distr = "rconst", const = 0))  
      }
    }
    
    return(nodes)
  }
  
  PROB<-generate_action_from_vector(PROB)
  ate <- rep(0,TM+2) 
  for(s in 0:(TM)){
    prob = PROB[(s+1):length(PROB)]
    
    Dact <- D + action(name = "do", nodes = prob)
    dat_do <- sim(Dact, n = accuracy, actions = "do",wide = T)$do
    ate[s+1] <- mean(dat_do[[paste('Y',TM+1,sep = '_')]])
    # ate[s+1] <- eval.target(DAG = set.targetE(Dact, outcome=paste("Y", TM+1, sep = "_"), param="do"), 
    #                         n=accuracy,actions = "do")$res
  }
  ate[TM+2] <- mean(simLIV(accuracy,D)$Y,na.rm = T)
  return(ate)
}


plot_LIV <- function(res,D,PROB){
  TM<-nrow(res)-2
  res$ate=ATE(D,TM,PROB)
  fig <- ggplot(res, aes(x = time)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CI"), alpha = 0.2) +
    geom_line(aes(y = est, color = "Estimate"), linewidth = 1) +
    geom_point(aes(y = est, color = "Estimate")) +
    geom_line(aes(y = ate, color = "True ATE"), linewidth = 1) +
    geom_point(aes(y = ate, color = "True ATE")) +
    scale_color_manual(
      name = "Lines",
      values = c("Estimate" = "blue", "True ATE" = "red")
    ) +
    scale_fill_manual(
      name = "Confidence Band",
      values = c("95% CI" = "blue")
    ) +
    labs(
      x = "Time", y = "Estimate", 
      title = "Estimate with 95% Confidence Band"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  return(fig)
}
