local({
  alpha <- getOption("AIV.summary.alpha", 0)
  sample_sizes <- c(5000L, 10000L)
  dgp_indices <- 5:6
  dgp_labels <- c("Y1, A3", "Y2, A3")
  dgp_tex <- c("\\ref{Y1}, \\ref{A3}", "\\ref{Y2}, \\ref{A3}")
  type_values <- c(-1, 0, 1)
  type_names <- c("mean_m1", "mean_0", "mean_1")
  metric_labels <- c("Bias", "SE", "SD", "CR")
  results_dir <- file.path(getwd(), "results")

  if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) ||
      !is.finite(alpha) || alpha < 0 || alpha >= 0.5) {
    stop("`AIV.summary.alpha` must be a number in [0, 0.5).", call. = FALSE)
  }

  robust_sd <- function(x) {
    as.numeric((quantile(x, 0.75, names = FALSE) -
                  quantile(x, 0.25, names = FALSE)) /
                 (qnorm(0.75) - qnorm(0.25)))
  }

  summarize_matrix <- function(rec, truth = 0) {
    if (!is.matrix(rec) || !is.numeric(rec) || nrow(rec) < 2L ||
        ncol(rec) < 2L) {
      stop("Each result must be a numeric matrix with at least two rows and columns.",
           call. = FALSE)
    }
    named <- !is.null(colnames(rec)) &&
      all(c("est_aug", "sd_aug") %in% colnames(rec))
    estimate <- rec[, if (named) "est_aug" else 1L]
    se <- rec[, if (named) "sd_aug" else 2L]
    if (any(!is.finite(estimate)) || any(!is.finite(se)) || any(se <= 0)) {
      stop("Estimates must be finite and estimated standard errors positive.",
           call. = FALSE)
    }

    limits <- quantile(estimate, c(alpha, 1 - alpha), names = FALSE)
    keep <- estimate >= limits[1] & estimate <= limits[2]
    c(
      Bias = abs(mean(estimate[keep]) - truth),
      SE = mean(se[keep]),
      SD = sd(estimate[keep]),
      CR = mean(abs((estimate - truth) / se) < qnorm(0.975))
    )
    # c(
    #   Bias = abs(median(estimate[keep]) - truth),
    #   SE = median(se[keep]),
    #   SD = robust_sd(estimate[keep]),
    #   CR = mean(abs((estimate - truth) / se) < qnorm(0.975))
    # )
  }

  format_value <- function(value, metric) {
    digits <- if (metric == "CR") 3L else 4L
    value <- round(value, digits)
    if (value == 0) value <- 0
    sub("^(-?)0\\.", "\\1.", formatC(value, format = "f", digits = digits))
  }

  load_summary <- function(n, index, type_value, type_name, estimator) {
    path <- file.path(results_dir, sprintf("data(%d)", n),
                      sprintf("%s_%d.RData", type_name, index))
    if (!file.exists(path)) stop("Missing result file: ", path, call. = FALSE)

    loaded <- new.env(parent = emptyenv())
    load(path, envir = loaded)
    if (!all(c("REC_unknown", "REC_known", "settings") %in% ls(loaded))) {
      stop("Missing REC_unknown, REC_known, or settings in: ", path,
           call. = FALSE)
    }

    settings <- loaded$settings
    expected_pi <- if (type_value == 0) "Z^2" else "Z"
    if (!identical(as.integer(settings$index), as.integer(index)) ||
        !identical(as.integer(settings$n), as.integer(n)) ||
        !isTRUE(as.numeric(settings$type) == type_value) ||
        !identical(as.character(settings$pi_known), expected_pi)) {
      stop("Result metadata does not match file: ", path, call. = FALSE)
    }

    rec <- if (estimator == "Adaptive") loaded$REC_unknown else loaded$REC_known
    list(summary = summarize_matrix(rec, truth = 0), repetitions = nrow(rec))
  }

  column_keys <- unlist(lapply(sample_sizes, function(n) {
    c(paste(n, "Adaptive", type_names, sep = "_"),
      paste(n, "Prespecified", type_names, sep = "_"))
  }), use.names = FALSE)
  table <- data.frame(
    DGP = rep(dgp_labels, each = length(metric_labels)),
    Metric = rep(metric_labels, length(dgp_indices)),
    stringsAsFactors = FALSE
  )
  for (key in column_keys) table[[key]] <- NA_real_

  repetition_counts <- integer()
  for (d in seq_along(dgp_indices)) {
    index <- dgp_indices[d]
    rows <- (d - 1L) * length(metric_labels) + seq_along(metric_labels)
    for (n in sample_sizes) {
      for (estimator in c("Adaptive", "Prespecified")) {
        for (j in seq_along(type_values)) {
          result <- load_summary(n, index, type_values[j], type_names[j], estimator)
          key <- paste(n, estimator, type_names[j], sep = "_")
          table[rows, key] <- result$summary[metric_labels]
          repetition_counts <- c(repetition_counts, result$repetitions)
        }
      }
    }
  }
  if (length(unique(repetition_counts)) != 1L) {
    stop("The result files have inconsistent repetition counts.", call. = FALSE)
  }
  repetitions <- unique(repetition_counts)

  tex <- c(
    sprintf("%% repetitions = %d; alpha = %.3f.", repetitions, alpha),
    "\\begin{tabular}{cc| ccc |ccc |ccc |ccc}",
    "\\toprule",
    "\\multirow{3}{*}{DGP}&\\multirow{3}{*}{Metric}&",
    "\\multicolumn{6}{c}{$n=5000$}&\\multicolumn{6}{|c}{$n=10000$}\\\\",
    "\\cmidrule(lr){3-8}\\cmidrule(lr){9-14}",
    "&&\\multicolumn{3}{c}{Adaptive weight}&\\multicolumn{3}{|c}{Prespecified}&",
    "\\multicolumn{3}{|c}{Adaptive weight}&\\multicolumn{3}{|c}{Prespecified}\\\\",
    paste0("&& $\\mathbb{E}[Y(-1)]$ & $\\mathbb{E}[Y(0)]$ & $\\mathbb{E}[Y(1)]$ &",
           " $\\mathbb{E}[Y(-1)]$ & $\\mathbb{E}[Y(0)]$ & $\\mathbb{E}[Y(1)]$ &",
           " $\\mathbb{E}[Y(-1)]$ & $\\mathbb{E}[Y(0)]$ & $\\mathbb{E}[Y(1)]$ &",
           " $\\mathbb{E}[Y(-1)]$ & $\\mathbb{E}[Y(0)]$ & $\\mathbb{E}[Y(1)]$\\\\"),
    "\\midrule"
  )

  for (d in seq_along(dgp_indices)) {
    tex <- c(tex, sprintf("\\multirow{4}{*}{%s}", dgp_tex[d]))
    for (metric in metric_labels) {
      row <- (d - 1L) * length(metric_labels) + match(metric, metric_labels)
      values <- vapply(column_keys, function(key) {
        format_value(table[[key]][row], metric)
      }, character(1))
      suffix <- if (metric == "CR" && d < length(dgp_indices)) " \\\\[5pt]" else " \\\\"
      tex <- c(tex, paste0("&", metric, " &", paste(values, collapse = " &"), suffix))
    }
  }
  tex <- c(tex, "\\bottomrule", "\\end{tabular}")

  output_path <- file.path(results_dir, "Table3.tex")
  writeLines(tex, output_path)
  message(sprintf("Robust additional table (%d repetitions) saved to %s",
                  repetitions, output_path))
  print(table, row.names = FALSE, digits = 4)
  invisible(table)
})
