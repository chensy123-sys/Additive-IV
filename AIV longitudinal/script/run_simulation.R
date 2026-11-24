run_simulation <- function(
    n, action_list,
    Looptime = getOption("LIV.Looptime", 500L),
    seed = getOption("LIV.seed", 2025L),
    project_dir = getwd(),
    output_dir = getOption("LIV.output_dir",
                           file.path(project_dir, "results", sprintf("data(%d)", n)))) {
  check_integer <- function(value, lower, name) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value != floor(value) || value < lower ||
        value > .Machine$integer.max) {
      stop(sprintf("`%s` must be an integer >= %d.", name, lower), call. = FALSE)
    }
  }
  check_integer(n, 5L, "n")
  check_integer(Looptime, 1L, "Looptime")
  check_integer(seed, 0L, "seed")
  if (!is.numeric(action_list) || !length(action_list) ||
      anyNA(action_list) || any(!action_list %in% c(0, 1))) {
    stop("`action_list` must contain only 0 or 1, in time order.", call. = FALSE)
  }

  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  # Rscript may start from script/, the project root, or an unrelated directory.
  # The CLI supplies the script's own directory; resolve its parent when needed.
  candidates <- unique(c(project_dir, dirname(project_dir)))
  valid <- vapply(candidates, function(directory) {
    all(file.exists(file.path(directory, c("script/par.R", "LIV.R"))))
  }, logical(1L))
  if (!any(valid)) {
    stop("Cannot locate the project root containing script/par.R and LIV.R.", call. = FALSE)
  }
  project_dir <- candidates[which(valid)[1L]]
  dependencies <- file.path(project_dir, c("script/par.R", "LIV.R"))
  model <- new.env(parent = globalenv())
  suppressMessages({
    sys.source(dependencies[1L], envir = model)
    sys.source(dependencies[2L], envir = model)
  })
  if (length(action_list) != model$TM + 1L) {
    stop(sprintf("`action_list` must contain %d values for times 0 through %d.",
                 model$TM + 1L, model$TM), call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  action_name <- paste(action_list, collapse = "")
  output_files <- c(
    RES = file.path(output_dir, paste0("RES_", action_name, ".RData")),
    RES_adap = file.path(output_dir, paste0("RES_adap_", action_name, ".RData"))
  )

  set.seed(seed)
  RES <- array(0, dim = c(model$TM + 2L, 5L, Looptime))
  RES_adap <- array(0, dim = c(model$TM + 2L, 5L, Looptime))
  loop <- 0L
  attempts <- 0L
  while (loop < Looptime) {
    attempts <- attempts + 1L
    message(sprintf("Monte Carlo %d/%d (attempt %d)", loop + 1L, Looptime, attempts))
    # Retry the entire replication with fresh data if either estimator fails.
    # Keep the RNG advancing: resetting the seed here would repeat a failed draw.
    result <- tryCatch({
      data_LIV <- model$simLIV(n, model$D)
      fit_LIV <- model$DML_LIV(data_LIV, action_list = action_list, K = 5L)
      fit_LIV_adap <- model$DML_LIV_adap(data_LIV, action_list = action_list, K = 5L)
      result <- list(
        RES = as.matrix(model$predict_LIV(fit_LIV)),
        RES_adap = as.matrix(model$predict_LIV(fit_LIV_adap))
      )
      expected_dim <- c(model$TM + 2L, 5L)
      if (!all(vapply(result, function(x) identical(dim(x), as.integer(expected_dim)),
                      logical(1L)))) {
        stop("Prediction results must have TM+2 rows and 5 columns.", call. = FALSE)
      }
      result
    }, error = function(e) {
      message(sprintf("Attempt %d failed; completed %d/%d. Retrying: %s",
                      attempts, loop, Looptime, conditionMessage(e)))
      NULL
    })
    if (is.null(result)) next

    # Commit both estimates only after the complete replication succeeds.
    loop <- loop + 1L
    RES[, , loop] <- result$RES
    RES_adap[, , loop] <- result$RES_adap
  }
  settings <- list(n = n, action_list = action_list, Looptime = Looptime, K = 5L,
                   seed = seed, attempts = attempts, failures = attempts - loop)
  save(RES, settings, file = output_files[["RES"]])
  save(RES_adap, settings, file = output_files[["RES_adap"]])
  message("Saved results to: ", output_dir)
  invisible(output_files)
}

run_simulation_cli <- function(args = commandArgs(trailingOnly = TRUE), project_dir = NULL) {
  usage <- paste(
    "Usage: Rscript run_simulation.R <n> <action_list>",
    "  n: integer >= 5 (K=5)",
    "  action_list: one 0/1 value per treatment time, in time order.",
    "  Separate actions with commas or spaces; their number must match TM+1 in par.R.",
    "From script/: Rscript run_simulation.R 2000 0,0",
    "From project root: Rscript script/run_simulation.R 2000 0,0",
    "Space-separated actions are also accepted: 2000 0 0",
    "Results are always written under the project root:",
    "  results/data(n)/RES_<actions>.RData and RES_adap_<actions>.RData",
    sep = "\n"
  )
  if (length(args) == 1L && args %in% c("-h", "--help")) {
    cat(usage, "\n", sep = "")
    return(invisible(NULL))
  }
  if (length(args) < 2L) stop(usage, call. = FALSE)
  n <- suppressWarnings(as.numeric(args[1L]))
  if (!grepl("^[0-9]+$", args[1L]) || !is.finite(n) ||
      n < 5 || n > .Machine$integer.max) {
    stop("`n` must be an integer between 5 and .Machine$integer.max.", call. = FALSE)
  }
  action_args <- args[-1L]
  if (any(!nzchar(trimws(action_args))) || any(grepl(",[[:space:]]*$", action_args))) {
    stop("Each action value and each value following a comma must be 0 or 1.", call. = FALSE)
  }
  actions <- trimws(unlist(strsplit(action_args, ",", fixed = TRUE), use.names = FALSE))
  if (!length(actions) || any(!actions %in% c("0", "1"))) {
    stop("`action_list` must contain only 0 or 1, separated by commas or spaces.", call. = FALSE)
  }
  if (is.null(project_dir)) {
    script_arg <- commandArgs(trailingOnly = FALSE)
    script_arg <- script_arg[startsWith(script_arg, "--file=")]
    project_dir <- if (length(script_arg)) {
      dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
    } else getwd()
  }
  run_simulation(n = as.integer(n), action_list = as.numeric(actions), project_dir = project_dir)
}

# source() defines functions only; Rscript starts the requested simulation.
if (sys.nframe() == 0L) {
  tryCatch(run_simulation_cli(), error = function(e) {
    message(conditionMessage(e))
    quit(save = "no", status = 1L)
  })
}
