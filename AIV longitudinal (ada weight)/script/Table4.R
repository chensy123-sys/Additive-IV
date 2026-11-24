summary_source_dir <- local({
  for (frame in rev(sys.frames())) {
    source_file <- frame$ofile
    if (is.character(source_file) && length(source_file) == 1L &&
        !is.na(source_file) && file.exists(source_file)) {
      break
    }
    source_file <- NULL
  }
  if (!is.null(source_file)) dirname(normalizePath(source_file))
})

find_summary_project <- function(directory) {
  candidates <- c(directory, dirname(directory))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "LIV.R")) &&
        file.exists(file.path(candidate, "script", "par.R")) &&
        dir.exists(file.path(candidate, "results"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot find results, LIV.R, and script/par.R under the same project root. ",
       "Run from the project root or supply project_dir.", call. = FALSE)
}

simulation_sizes <- function(results_dir) {
  dirs <- basename(list.dirs(results_dir, full.names = TRUE, recursive = FALSE))
  dirs <- grep("^data[(][0-9]+[)]$", dirs, value = TRUE)
  sort(as.integer(sub("^data[(]([0-9]+)[)]$", "\\1", dirs)))
}

render_simulation_table <- function(sizes = NULL,
                                    project_dir = if (is.null(summary_source_dir)) getwd()
                                                  else summary_source_dir,
                                    accuracy = 1e6) {
  project_dir <- find_summary_project(normalizePath(project_dir, mustWork = TRUE))
  results_dir <- file.path(project_dir, "results")
  if (is.null(sizes)) sizes <- simulation_sizes(results_dir)
  if (!is.numeric(sizes) || !length(sizes) || anyNA(sizes) ||
      any(!is.finite(sizes) | sizes < 1 | sizes != floor(sizes) |
          sizes > .Machine$integer.max)) {
    stop("Supply positive integer sample sizes or add results/data(n) folders.",
         call. = FALSE)
  }
  sizes <- sort(unique(as.integer(sizes)))

  # The second row leaves A_0 natural and only intervenes at time 1.
  interventions <- c("$(0,0)$", "$(0,1)$", "$(1,0)$", "$(1,1)$",
                     "$(A_0,0)$", "$(A_0,1)$")
  files <- c("00", "01", "10", "11")
  result_ids <- c(files, "00", "01")
  result_rows <- c(1L, 1L, 1L, 1L, 2L, 2L)
  metrics <- c("Bias", "SE", "SD", "CR")
  paths <- outer(sizes, files, function(n, id) {
    file.path(results_dir, sprintf("data(%d)", n), paste0("RES_adap_", id, ".RData"))
  })
  available <- matrix(file.exists(paths), nrow = length(sizes))
  if (!any(available)) {
    stop("No adaptive simulation results found for the requested sample sizes.",
         call. = FALSE)
  }
  if (any(!available)) {
    warning("Missing adaptive results (shown as --):\n",
            paste(paths[!available], collapse = "\n"), call. = FALSE)
  }

  model <- new.env(parent = globalenv())
  invisible(capture.output(suppressMessages({
    sys.source(file.path(project_dir, "script", "par.R"), envir = model)
    sys.source(file.path(project_dir, "LIV.R"), envir = model)
  })))
  if (model$TM != 1L) stop("This table requires two treatment times.", call. = FALSE)
  set.seed(2025)
  truth <- setNames(lapply(files, function(id) {
    actions <- as.integer(strsplit(id, "", fixed = TRUE)[[1L]])
    invisible(capture.output(
      values <- suppressMessages(model$ATE(model$D, model$TM, actions, accuracy = accuracy))
    ))
    values
  }), files)

  latex <- c(
    "\\begin{tabular}{ccc|cccccc}",
    "\\toprule",
    "& \\multirow{2}{*}{Size} & \\multirow{2}{*}{Metric} & \\multicolumn{6}{c}{Intervention} \\\\",
    paste("&&&", paste(interventions, collapse = " & "), "\\\\"),
    "\\midrule"
  )
  for (i in seq_along(sizes)) {
    results <- setNames(vector("list", length(files)), files)
    for (j in which(available[i, ])) {
      data <- new.env(parent = emptyenv())
      load(paths[i, j], envir = data)
      res <- data$RES_adap
      if (!is.numeric(res) || length(dim(res)) != 3L ||
          dim(res)[1L] < 2L || dim(res)[2L] < 5L || dim(res)[3L] < 1L) {
        stop("Invalid RES_adap array in ", paths[i, j], call. = FALSE)
      }
      results[[j]] <- res
    }

    values <- vapply(seq_along(interventions), function(j) {
      id <- result_ids[j]
      res <- results[[id]]
      if (is.null(res)) return(rep(NA_real_, length(metrics)))
      row <- result_rows[j]
      target <- truth[[id]][row]
      est <- res[row, 2L, ]
      c(abs(mean(est) - target),
        mean(res[row, 3L, ]),
        stats::sd(est),
        mean(res[row, 4L, ] <= target & target <= res[row, 5L, ]))
    }, numeric(length(metrics)))

    for (j in seq_along(metrics)) {
      size <- if (j == 1L) sprintf("\\multirow{4}{*}{%d}", sizes[i]) else ""
      precision <- if (metrics[j] == "CR") "%.3f" else "%.4f"
      cells <- ifelse(is.finite(values[j, ]),
                      sub("^0", "", sprintf(precision, values[j, ])), "--")
      ending <- if (j == length(metrics) && i < length(sizes)) "\\\\[5pt]" else "\\\\"
      latex <- c(latex, paste(paste(c("", size, metrics[j], cells), collapse = " & "), ending))
    }
  }
  c(latex, "\\bottomrule", "\\end{tabular}")
}

if (sys.nframe() == 0L) {
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  project_dir <- getwd()
  args <- character()
  if (length(script_arg)) {
    script <- sub("^--file=", "", script_arg[1L])
    if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
    if (file.exists(script)) {
      project_dir <- dirname(normalizePath(script, mustWork = TRUE))
      args <- commandArgs(trailingOnly = TRUE)
    }
  }
  if (length(args) == 1L && args %in% c("-h", "--help")) {
    cat("Usage: Rscript summary.R [sample_size ...]\n")
  } else {
    sizes <- if (length(args)) suppressWarnings(as.numeric(args)) else NULL
    writeLines(render_simulation_table(sizes, project_dir))
  }
}
