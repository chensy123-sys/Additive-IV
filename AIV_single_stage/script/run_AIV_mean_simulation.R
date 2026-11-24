run_AIV_mean_simulation <- function(
    index,
    n,
    types,
    Looptime = getOption("AIV.Looptime", 500L),
    K = getOption("AIV.K", 5L),
    seed = getOption("AIV.seed", 2025L),
    output_dir = getOption("AIV.output_dir", file.path("results", sprintf("data(%d)", n))),
    max_attempts = 3L) {
  check_integer <- function(x, lower, name) {
    if (length(x) != 1L || !is.numeric(x) || is.na(x) ||
        !is.finite(x) || x < lower || x > .Machine$integer.max ||
        x != floor(x)) {
      stop(sprintf("`%s` must be an integer >= %d.", name, lower), call. = FALSE)
    }
  }
  check_integer(index, 1L, "index")
  if (!index %in% 1:6) stop("`index` must be between 1 and 6.", call. = FALSE)
  check_integer(n, 2L, "n")
  check_integer(Looptime, 1L, "Looptime")
  check_integer(K, 2L, "K")
  if (K > n) stop("`K` cannot exceed `n`.", call. = FALSE)
  check_integer(seed, 0L, "seed")
  check_integer(max_attempts, 1L, "max_attempts")
  binary_types <- c("ATE", "treat", "control")
  ternary_types <- c(-1, 0, 1)
  if (!length(types) || anyNA(types) || anyDuplicated(types)) {
    stop("`types` must be non-empty, with no missing or duplicate values.", call. = FALSE)
  }
  if (index <= 4L) {
    if (!is.character(types) || !all(types %in% binary_types)) {
      stop("For index 1--4, `types` must contain only \"ATE\", \"treat\", or \"control\".",
           call. = FALSE)
    }
    type_names <- types
  } else {
    if (!is.numeric(types) || !all(types %in% ternary_types)) {
      stop("For index 5--6, `types` must contain only -1, 0, or 1.", call. = FALSE)
    }
    type_names <- c("mean_m1", "mean_0", "mean_1")[match(types, ternary_types)]
  }

  project_dir <- normalizePath(getwd(), mustWork = TRUE)
  aiv_file <- file.path(project_dir, "AIV.R")
  if (!file.exists(aiv_file)) {
    stop("Run this script from the AIV_single_stage directory.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  if (index <= 4L) {
    output_files <- setNames(
      file.path(output_dir, sprintf("%s%d.RData", type_names, index)),
      type_names
    )
  } else {
    output_files <- setNames(
      file.path(output_dir, sprintf("%s_%d.RData", type_names, index)),
      type_names
    )
  }
  # if (any(file.exists(output_files))) {
  #   stop("Result files already exist. Choose another AIV.output_dir.", call. = FALSE)
  # }

  aiv_env <- new.env(parent = globalenv())
  sys.source(aiv_file, envir = aiv_env)
  metrics <- c("est_aug", "sd_aug", "est_or", "est_ipw")
  records <- setNames(lapply(type_names, function(type_name) {
    list(
      unknown = matrix(NA_real_, Looptime, length(metrics),
                       dimnames = list(NULL, metrics)),
      known = matrix(NA_real_, Looptime, length(metrics),
                     dimnames = list(NULL, metrics))
    )
  }), type_names)

  set.seed(seed)
  data_seeds <- sample.int(100000000L, Looptime, replace = TRUE)
  fold_seeds <- sample.int(100000000L, Looptime, replace = TRUE)

  for (i in seq_len(Looptime)) {
    set.seed(data_seeds[i])
    sim_dat <- aiv_env$sim_AIV(n, index = index, print_state = FALSE)
    success <- FALSE
    for (attempt in seq_len(max_attempts)) {
      fit_results <- tryCatch({
        setNames(lapply(seq_along(types), function(j) {
          type <- types[j]
          set.seed(fold_seeds[i] + attempt - 1L)
          unknown <- aiv_env$AIV_unknown_spline(sim_dat, type = type, K = K)
          set.seed(fold_seeds[i] + attempt - 1L)
          if (index >= 5L && identical(as.numeric(type), 0)) {
            known <- aiv_env$AIV_known_spline(
              sim_dat,
              type = type,
              K = K,
              pi = function(Z, L) {Z^2}
            )
          } else {
            known <- aiv_env$AIV_known_spline(sim_dat, type = type, K = K)
          }
          unknown <- unlist(unknown)[metrics]
          known <- unlist(known)[metrics]
          if (any(!is.finite(c(unknown, known)))) {
            stop(sprintf("Non-finite estimates for type=%s.", type), call. = FALSE)
          }
          list(unknown = unknown, known = known)
        }), type_names)
      }, error = function(e) e)
      if (!inherits(fit_results, "error")) {
        success <- TRUE
        break
      }
      message(sprintf("DGP %d, repetition %d, attempt %d: %s",
                      index, i, attempt, conditionMessage(fit_results)))
    }
    if (!success) {
      stop(sprintf("DGP %d, repetition %d failed after %d attempts.",
                   index, i, max_attempts), call. = FALSE)
    }
    for (type_name in type_names) {
      records[[type_name]]$unknown[i, ] <- fit_results[[type_name]]$unknown
      records[[type_name]]$known[i, ] <- fit_results[[type_name]]$known
      message(sprintf("unknown: type: %s est: %.4f, std: %.4f", type_name, records[[type_name]]$unknown[i, 1], records[[type_name]]$unknown[i, 2]))
      message(sprintf("known: type: %s est: %.4f, std: %.4f", type_name, records[[type_name]]$known[i, 1], records[[type_name]]$known[i, 2]))
    }
    message(sprintf("DGP %d: %d/%d (types: %s completed)",
                    index, i, Looptime, paste(types, collapse = ", ")))
    
  }

  for (j in seq_along(types)) {
    type_name <- type_names[j]
    REC_unknown <- records[[type_name]]$unknown
    REC_known <- records[[type_name]]$known
    pi_known <- if (index >= 5L && identical(as.numeric(types[j]), 0)) {
      "Z^2"
    } else {
      "Z"
    }
    settings <- list(index = index, target = types[j], type = types[j], n = n,
                     pi_known = pi_known, Looptime = Looptime, K = K, seed = seed,
                     data_seeds = data_seeds, fold_seeds = fold_seeds)
    save(REC_unknown, REC_known, settings, file = output_files[[type_name]])
  }
  message("Saved mean potential outcome results to: ", output_dir)
  invisible(output_files)
}

run_AIV_mean_simulation_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste(
    "Usage: Rscript script/run_AIV_mean_simulation.R <index> <n> <type> [<type> ...]",
    "  index: 1--6; n: integer >= 5 (default K = 5)",
    "  type for index 1--4: ATE, treat, control, or all",
    "  type for index 5--6: -1, 0, 1, or all",
    "  Multiple types: separate with spaces or commas; use all by itself.",
    "Run from the AIV_single_stage directory.",
    sep = "\n"
  )
  if (length(args) == 1L && args %in% c("-h", "--help")) {
    cat(usage, "\n", sep = "")
    return(invisible(NULL))
  }
  if (length(args) < 3L) stop(usage, call. = FALSE)

  parse_integer <- function(value, name, lower, upper = .Machine$integer.max) {
    number <- suppressWarnings(as.numeric(value))
    if (!grepl("^[0-9]+$", value) || !is.finite(number) ||
        number < lower || number > upper) {
      stop(sprintf("`%s` must be an integer between %d and %d.",
                   name, lower, upper), call. = FALSE)
    }
    as.integer(number)
  }
  index <- parse_integer(args[1L], "index", 1L, 6L)
  n <- parse_integer(args[2L], "n", 2L)
  allowed_types <- if (index <= 4L) c("ATE", "treat", "control") else c("-1", "0", "1")
  type_args <- args[-c(1L, 2L)]
  # strsplit() drops a trailing empty field, so reject trailing commas explicitly.
  if (any(grepl(",[[:space:]]*$", type_args))) {
    stop("Each comma must be followed by a type value.", call. = FALSE)
  }
  types <- trimws(unlist(strsplit(type_args, ",", fixed = TRUE), use.names = FALSE))
  if (!length(types) || any(!nzchar(trimws(type_args))) || any(!nzchar(types))) {
    stop("`type` values must not be empty.", call. = FALSE)
  }
  if ("all" %in% types) {
    if (length(types) != 1L) stop("Use `all` by itself.", call. = FALSE)
    types <- allowed_types
  }
  if (!all(types %in% allowed_types)) {
    stop(sprintf("For index %d, `type` values must be chosen from: %s, all.",
                 index, paste(allowed_types, collapse = ", ")), call. = FALSE)
  }
  if (index >= 5L) types <- as.numeric(types)

  run_AIV_mean_simulation(index = index, n = n, types = types)
}

# source() only defines the functions; Rscript invokes the command-line entry point.
if (sys.nframe() == 0L) {
  tryCatch(
    run_AIV_mean_simulation_cli(),
    error = function(e) {
      message(conditionMessage(e))
      quit(save = "no", status = 1L)
    }
  )
}
