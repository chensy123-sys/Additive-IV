rm(list = ls())
library(dplyr)
library(tibble)
library(splines)
library(mgcv)



sim_AIV <- function(n, index, print_state = TRUE) {
  if (length(index) != 1L || !is.numeric(index) || is.na(index) ||
      !is.finite(index) || index != as.integer(index) || !index %in% 1:6) {
    stop("`index` must be one of 1, 2, 3, 4, 5, or 6.", call. = FALSE)
  }
  index <- as.integer(index)
  
  # index = 1: (Y1, A1); index = 2: (Y2, A1)
  # index = 3: (Y1, A2); index = 4: (Y2, A2)
  # index = 5: signed treatment with Y = 2U - 2L + 4AL
  # index = 6: signed treatment with Y2
  treatment_model <- c(1L, 1L, 2L, 2L, 3L, 3L)[index]
  outcome_model <- c(1L, 2L, 1L, 2L, 1L, 2L)[index]
  
  L <- runif(n, -1, 1)
  U <- runif(n, -1, 1)
  Z <- (L + sin(3 * L)) + 4*runif(n,min = -1,max = 1)
  
  if (treatment_model == 1L) {
    prob_A <- 0.7 * pnorm(- ((Z - 2)^2 - 6) + 2 * L) + 0.3 * pnorm(3 * U - L)
    A <- rbinom(n, 1, prob_A)
  } else if (treatment_model == 2L) {
    prob_A <- plogis(0.5 * Z * (2 + U) - L + 2 * U)
    A <- rbinom(n, 1, prob_A)
  } else {
    e <- runif(n)
    A <- (as.integer(3*Z + L >= 4) - as.integer(3*Z + L <= -4)) *
      as.integer(e <= 0.8) +
      (as.integer(U >= 1 / 3) - as.integer(U <= -1 / 3)) *
      as.integer(e >= 0.8)
  }
  
  if (outcome_model == 1L) {
    Y0 <- 2 * U - 2 * L
    Y1 <- 2 * U + 2 * L
  } else {
    Y0 <- 3 * cos(2 * U) - 3 * cos(2 * L)
    Y1 <- 3 * sin(2 * U) + 2 * L
  }
  Y <- Y0 + A * (Y1 - Y0)
  
  mydat <- tibble(L = L, Z = Z, A = A, Y = Y)
  
  if (print_state) {
    outcome_fit <- lm(Y ~ A * ns(L, df = 10), data = mydat)
    adjusted_ate <- mean(predict(outcome_fit,
                                 newdata = data.frame(A = 1, L = L))) -
      mean(predict(outcome_fit,
                   newdata = data.frame(A = 0, L = L)))
    
    cat("DGP index:", index,
        "\nCor(U, A):", cor(U, A),
        "\nCor(Z, A):", cor(Z, A),
        "\nTrue ATE:", mean(Y1 - Y0),
        "\nMean Y(1):", mean(Y1),
        "\nMean Y(0):", mean(Y0),
        "\nTreat as RCT:", mean(Y[A == 1]) - mean(Y[A == 0]),
        "\nAdjust for L:", adjusted_ate, "\n")
  }
  
  return(mydat)
}











.aiv_parallel_folds <- function(fit_fold, K) {
  worker_fun <- fit_fold
  environment(worker_fun) <- list2env(list(
    mydat = get("mydat", envir = environment(fit_fold)),
    split_indices = get("split_indices", envir = environment(fit_fold))
  ), parent = globalenv())
  
  cluster <- parallel::makePSOCKcluster(K)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(lib_paths) {
    .libPaths(lib_paths)
    suppressPackageStartupMessages({
      library(dplyr)
      library(mgcv)
    })
    NULL
  }, .libPaths())
  
  results <- parallel::parLapply(cluster, seq_len(K), function(k, fold_fun) {
    tryCatch(
      list(ok = TRUE, result = fold_fun(k)),
      error = function(e) list(ok = FALSE, error = conditionMessage(e))
    )
  }, fold_fun = worker_fun)
  
  failed <- which(!vapply(results, function(x) isTRUE(x$ok), logical(1)))
  if (length(failed)) {
    messages <- vapply(failed, function(k) {
      paste0("Fold ", k, ": ", results[[k]]$error)
    }, character(1))
    stop(paste(messages, collapse = "\n"), call. = FALSE)
  }
  lapply(results, function(x) x$result)
}

AIV_unknown_spline <- function(data, K = 2, type = "ATE") {
  n <- nrow(data)
  if (length(K) != 1L || !is.numeric(K) || is.na(K) || !is.finite(K) ||
      K != as.integer(K) || K < 2L || K > n) {
    stop("`K` must be an integer between 2 and the sample size.", call. = FALSE)
  }
  K <- as.integer(K)
  
  cross_fit_split_index <- function(n, K) {
    folds <- sample(rep(seq_len(K), length.out = n))
    lapply(seq_len(K), function(k) {
      test_idx <- which(folds == k)
      list(
        train_idx = setdiff(seq_len(n), test_idx),
        test_idx = test_idx
      )
    })
  }
  
  mydat <- data
  if (type == "treat") {
    mydat <- mydat %>% mutate(Y = Y * A)
  } else if (type == "control") {
    mydat <- mydat %>% mutate(Y = Y * (A - 1))
  } else if (type == "ATE"){
    mydat <- mydat
  } else if (is.numeric(type)){
    mydat <- mydat %>% mutate(
      Y = Y * as.numeric(A == type), 
      A = as.numeric(A == type)
    )
  }
  
  split_indices <- cross_fit_split_index(n, K)
  fit_fold <- function(k) {
    train_index <- split_indices[[k]][["train_idx"]]
    test_index <- split_indices[[k]][["test_idx"]]
    train_data <- mydat[train_index, , drop = FALSE]
    test_data <- mydat[test_index, , drop = FALSE]
    
    pi_model <- gam(A ~ s(Z, L, bs = "tp"),
                    data = train_data, family = binomial())
    delta_model <- gam(A ~ s(L, bs = "tp"),
                       data = train_data, family = binomial())
    xi_model <- gam(Y ~ s(Z, L, bs = "tp"),
                    data = train_data, family = gaussian())
    eta_model <- gam(Y ~ s(L, bs = "tp"),
                     data = train_data, family = gaussian())
    
    pi_pred_train <- as.vector(predict(pi_model, newdata = train_data,
                                       type = "response"))
    delta_pred_train <- as.vector(predict(delta_model, newdata = train_data,
                                          type = "response"))
    xi_pred_train <- as.vector(predict(xi_model, newdata = train_data,
                                       type = "response"))
    eta_pred_train <- as.vector(predict(eta_model, newdata = train_data,
                                        type = "response"))
    
    train_data$kappa_score <-
      (train_data$A - delta_pred_train)^2 -
      (train_data$A - pi_pred_train)^2
    train_data$Gamma_score <-
      train_data$Y * (pi_pred_train - delta_pred_train) +
      xi_pred_train * (train_data$A - pi_pred_train) -
      eta_pred_train * (train_data$A - delta_pred_train)
    
    kappa_model <- gam(kappa_score ~ s(L, bs = "tp"),
                       data = train_data, family = gaussian())
    Gamma_model <- gam(Gamma_score ~ s(L, bs = "tp"),
                       data = train_data, family = gaussian())
    
    pi_pred <- as.vector(predict(pi_model, newdata = test_data,
                                 type = "response"))
    delta_pred <- as.vector(predict(delta_model, newdata = test_data,
                                    type = "response"))
    xi_pred <- as.vector(predict(xi_model, newdata = test_data,
                                 type = "response"))
    eta_pred <- as.vector(predict(eta_model, newdata = test_data,
                                  type = "response"))
    kappa_pred <- as.vector(predict(kappa_model, newdata = test_data,
                                    type = "response"))
    Gamma_pred <- as.vector(predict(Gamma_model, newdata = test_data,
                                    type = "response"))
    
    kappa_pred_train <- as.vector(predict(kappa_model, newdata = train_data,
                                          type = "response"))
    # kappa_floor <- max(as.numeric(quantile(abs(kappa_pred_train), 0.1)), 1e-8)
    # den <- pmax(kappa_pred, kappa_floor)
    den <- kappa_pred
    gamma_pred <- Gamma_pred / den
    
    est_aug <- (pi_pred - delta_pred) / den * test_data$Y +
      gamma_pred / den * (
        den + (test_data$A - pi_pred)^2 -
          (test_data$A - delta_pred)^2
      ) +
      1 / den * (
        xi_pred * (test_data$A - pi_pred) -
          eta_pred * (test_data$A - delta_pred)
      )
    
    list(
      test_idx = test_index,
      est_aug = est_aug,
      est_or = gamma_pred,
      est_ipw = test_data$Y * (pi_pred - delta_pred) / den
    )
  }
  
  fold_results <- .aiv_parallel_folds(fit_fold, K)
  
  EST_AUG <- numeric(n)
  EST_OR <- numeric(n)
  EST_IPW <- numeric(n)
  for (fold_result in fold_results) {
    idx <- fold_result$test_idx
    EST_AUG[idx] <- fold_result$est_aug
    EST_OR[idx] <- fold_result$est_or
    EST_IPW[idx] <- fold_result$est_ipw
  }
  
  list(
    est_aug = mean(EST_AUG),
    sd_aug = sd(EST_AUG) / sqrt(n),
    est_or = mean(EST_OR),
    est_ipw = mean(EST_IPW)
  )
}





AIV_known_spline <- function(data, K=2,pi=function(Z,L){return(Z)},type = 'ATE'){
  n <- nrow(data)
  if (length(K) != 1L || !is.numeric(K) || is.na(K) || !is.finite(K) ||
      K != as.integer(K) || K < 2L || K > n) {
    stop("`K` must be an integer between 2 and the sample size.", call. = FALSE)
  }
  K <- as.integer(K)
  cross_fit_split_index <- function(n, K = 5) {
    folds <- sample(rep(1:K, length.out = n))
    splits <- vector("list", K)
    
    for (k in 1:K) {
      test_idx <- which(folds == k)
      train_idx <- setdiff(1:n, test_idx)
      splits[[k]] <- list(train_idx = train_idx, test_idx = test_idx)
    }
    
    return(splits)
  }
  
  mydat <- data
  if (type == "treat") {
    mydat <- mydat %>% mutate(Y = Y * A)
  } else if (type == "control") {
    mydat <- mydat %>% mutate(Y = Y * (A - 1))
  } else if (type == "ATE"){
    mydat <- mydat
  } else if (is.numeric(type)){
    mydat <- mydat %>% mutate(
      Y = Y * as.numeric(A == type), 
      A = as.numeric(A == type)
    )
  }
  mydat <- mydat %>% mutate(Z=pi(Z,L))
  
  
  EST_AUG <- matrix(0,nrow(mydat),1)
  EST_OR <- matrix(0,nrow(mydat),1)
  EST_IPW <- matrix(0,nrow(mydat),1)
  split_indices <- cross_fit_split_index(n = nrow(mydat), K)
  fit_fold <- function(k) {
    train_index <- split_indices[[k]][["train_idx"]]
    test_index <- split_indices[[k]][["test_idx"]]
    
    rho_model <- gam(Z ~ s(L, bs = "tp"),data = mydat[train_index,], family = gaussian())
    delta_model <- gam(A ~ s(L, bs = "tp"),data = mydat[train_index,], family = binomial())
    eta_model <- gam(Y ~ s(L, bs = "tp"),data = mydat[train_index,], family = gaussian())
    
    rho_pred_train   <- as.vector(predict(rho_model, newdata = mydat[train_index,] %>% select(-Y,-A,-Z),
                                          type = 'response'))
    delta_pred_train <- as.vector(predict(delta_model, newdata = mydat[train_index,] %>% select(-Y,-A,-Z),
                                          type = 'response'))
    eta_pred_train   <- as.vector(predict(eta_model, newdata = mydat[train_index,] %>% select(-Y,-A,-Z),
                                          type = 'response'))
    
    kappa_model<- gam((Z-rho_pred_train)*(A-delta_pred_train) ~ s(L, bs = "tp"), data = mydat[train_index,], family = gaussian())
    zeta_model <- gam((Z-rho_pred_train)*(Y-eta_pred_train) ~ s(L, bs = "tp"),data = mydat[train_index,], family = gaussian())
    
    rho_pred  <- as.vector(predict(rho_model, newdata = mydat[test_index,] %>% select(L),type = 'response'))
    delta_pred<- as.vector(predict(delta_model, newdata = mydat[test_index,] %>% select(L),type = 'response'))
    eta_pred  <- as.vector(predict(eta_model, newdata = mydat[test_index,] %>% select(L),type = 'response'))
    kappa_pred<- as.vector(predict(kappa_model, newdata = mydat[test_index,] %>% select(L),type = 'response'))
    zeta_pred <- as.vector(predict(zeta_model, newdata = mydat[test_index,] %>% select(L),type = 'response'))
    
    kappa_pred_train <- as.vector(predict(kappa_model,
                                          newdata = mydat[train_index,], type = 'response'))
    bounds <- quantile(kappa_pred_train, c(0.1, 0.9), names = FALSE)
    den <- pmin(pmax(kappa_pred, bounds[1]), bounds[2])
    den <- ifelse(den < 0, -1, 1) * pmax(abs(den), 1e-8)
    gamma_pred <- zeta_pred/den
    
    est_aug <- 1/den * (mydat[test_index,]$Z - rho_pred) * (mydat[test_index,]$Y - eta_pred) + (
      1 - (mydat[test_index,]$Z - rho_pred) * (mydat[test_index,]$A - delta_pred) / den
    ) * gamma_pred
    
    
    list(test_idx = test_index, est_aug = est_aug, est_or = gamma_pred,
         est_ipw = mydat[test_index,]$Y * (mydat[test_index,]$Z-rho_pred)/den)
  }
  fold_results <- .aiv_parallel_folds(fit_fold, K)
  for (fold_result in fold_results) {
    idx <- fold_result$test_idx
    EST_AUG[idx,] <- fold_result$est_aug
    EST_OR[idx,] <- fold_result$est_or
    EST_IPW[idx,] <- fold_result$est_ipw
  }
  RES = list(est_aug = mean(EST_AUG), sd_aug = sd(EST_AUG)/sqrt(n),
             est_or = mean(EST_OR), est_ipw = mean(EST_IPW))
  return(RES)
}
