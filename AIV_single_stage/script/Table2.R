local({
  alpha <- 0.000
  sample_sizes <- getOption("AIV.summary.n", c(5000L, 10000L))
  results_dir <- file.path(getwd(), "results")
  targets <- c("treat", "control", "ATE")
  dgp_labels <- c("Y1, A1", "Y2, A1", "Y1, A2", "Y2, A2")
  dgp_tex <- c(
    "\\ref{Y1}, \\ref{A1}", "\\ref{Y2}, \\ref{A1}",
    "\\ref{Y1}, \\ref{A2}", "\\ref{Y2}, \\ref{A2}"
  )
  metric_labels <- c("Bias", "SE", "SD", "CR")
  column_labels <- c(
    "Adaptive_Treated", "Adaptive_Control", "Adaptive_ATE",
    "Prespecified_Treated", "Prespecified_Control", "Prespecified_ATE"
  )

  sd_robust <- function(x) {
    as.numeric((quantile(x, 0.75) - quantile(x, 0.25)) /
                 (qnorm(0.75) - qnorm(0.25)))
  }

  summarize_matrix <- function(rec, truth = 0) {
    if (!is.matrix(rec) || !is.numeric(rec) || ncol(rec) < 2L ||
        nrow(rec) < 2L) {
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
    #   SD = sd_robust(estimate[keep]),
    #   CR = mean(abs((estimate - truth) / se) < qnorm(0.975))
    # )
  }

  format_value <- function(value, metric) {
    digits <- if (metric == "CR") 3L else 4L
    value <- round(value, digits)
    if (value == 0) value <- 0
    sub("^(-?)0\\.", "\\1.", formatC(value, format = "f", digits = digits))
  }

  export_table <- function(n) {
    data_dir <- file.path(results_dir, sprintf("data(%d)", n))
    if (!dir.exists(data_dir)) {
      stop("Missing result directory: ", data_dir, call. = FALSE)
    }
    table <- data.frame(DGP = rep(dgp_labels, each = 4),
                        Metric = rep(metric_labels, 4),
                        stringsAsFactors = FALSE)
    for (column in column_labels) table[[column]] <- NA_real_
    repetition_counts <- integer()
    for (i in 1:4) {
      rows <- (i - 1) * 4 + 1:4
      for (j in seq_along(targets)) {
        path <- file.path(data_dir, sprintf("%s%d.RData", targets[j], i))
        if (!file.exists(path)) stop("Missing result file: ", path, call. = FALSE)
        loaded <- new.env(parent = emptyenv())
        load(path, envir = loaded)
        if (!all(c("REC_unknown", "REC_known") %in% ls(loaded))) {
          stop("Missing REC_unknown or REC_known in: ", path, call. = FALSE)
        }
        if (exists("settings", envir = loaded, inherits = FALSE)) {
          settings <- loaded$settings
          if (!identical(as.integer(settings$index), as.integer(i)) ||
              !identical(as.integer(settings$n), as.integer(n)) ||
              !identical(settings$target, targets[j])) {
            stop("Result metadata does not match file: ", path, call. = FALSE)
          }
        }
        repetition_counts <- c(repetition_counts, nrow(loaded$REC_unknown),
                               nrow(loaded$REC_known))
        # Under Y1 and Y2, E[Y(1)], E[Y(0)] and ATE are all zero.
        table[rows, column_labels[j]] <- summarize_matrix(loaded$REC_unknown)
        table[rows, column_labels[j + 3L]] <- summarize_matrix(loaded$REC_known)
      }
    }
    if (length(unique(repetition_counts)) != 1L) {
      stop("The result files have inconsistent repetition counts.", call. = FALSE)
    }
    repetitions <- repetition_counts[1]
    tex <- c(
      sprintf("%% n = %d; repetitions = %d; alpha = %.3f.", n, repetitions, alpha),
      "% Robust summaries: Bias = abs(median(est_aug) - truth); SE = median(sd_aug).",
      "% SD = IQR(est_aug)/(qnorm(.75)-qnorm(.25)); CR uses all repetitions.",
      "% Truth = 0 for Treated, Control and ATE; absolute Bias is reported.",
      "% Requires booktabs and multirow; Y1/Y2/A1/A2 labels resolve in main.tex.",
      "\\begin{tabular}{cc| ccc |ccc}",
      "\\toprule",
      "\\multirow{2}{*}{DGP}&\\multirow{2}{*}{Metric}&",
      "\\multicolumn{3}{c}{Adaptive weight}&\\multicolumn{3}{|c}{Prespecified ($\\pi(Z,L)=Z$)}\\\\",
      "&& Treated & Control & ATE & Treated & Control & ATE\\\\",
      "\\midrule"
    )
    for (i in 1:4) {
      tex <- c(tex, sprintf("\\multirow{4}{*}{%s}", dgp_tex[i]))
      for (metric in metric_labels) {
        row <- (i - 1) * 4 + match(metric, metric_labels)
        values <- vapply(column_labels, function(column) {
          format_value(table[[column]][row], metric)
        }, character(1))
        tex <- c(tex, paste0("&", metric, " &", paste(values, collapse = " &"), " \\\\"))
      }
      if (i < 4L) tex <- c(tex, "[5pt]")
    }
    tex <- c(tex, "\\bottomrule", "\\end{tabular}")
    prefix <- file.path(results_dir, sprintf("summary_%d", n))
    writeLines(tex, paste0(prefix, ".tex"))
    message(sprintf("n=%d: %d repetitions; robust summaries saved to %s.tex",
                    n, repetitions, prefix))
    print(table, row.names = FALSE, digits = 4)
    invisible(table)
  }

  export_combined_table <- function(tables) {
    sizes <- c(5000L, 10000L)
    tex <- c(
      sprintf("%% n = 5000 and n = 10000; alpha = %.3f.", alpha),
      "\\begin{tabular}{cc| ccc |ccc |ccc |ccc}",
      "\\toprule",
      "\\multirow{3}{*}{DGP}&\\multirow{3}{*}{Metric}&",
      "\\multicolumn{6}{c}{$n=5000$}&\\multicolumn{6}{|c}{$n=10000$}\\\\",
      "\\cmidrule(lr){3-8}\\cmidrule(lr){9-14}",
      "&&\\multicolumn{3}{c}{Adaptive weight}&\\multicolumn{3}{|c}{Prespecified ($\\pi(Z,L)=Z$)}&",
      "\\multicolumn{3}{|c}{Adaptive weight}&\\multicolumn{3}{|c}{Prespecified ($\\pi(Z,L)=Z$)}\\\\",
      paste0("&& Treated & Control & ATE & Treated & Control & ATE &",
             " Treated & Control & ATE & Treated & Control & ATE\\\\"),
      "\\midrule"
    )
    for (i in 1:4) {
      tex <- c(tex, sprintf("\\multirow{4}{*}{%s}", dgp_tex[i]))
      for (metric in metric_labels) {
        row <- (i - 1) * 4 + match(metric, metric_labels)
        values <- unlist(lapply(as.character(sizes), function(size) {
          vapply(column_labels, function(column) {
            format_value(tables[[size]][[column]][row], metric)
          }, character(1))
        }), use.names = FALSE)
        tex <- c(tex, paste0("&", metric, " &", paste(values, collapse = " &"), " \\\\"))
      }
      if (i < 4L) tex <- c(tex, "[5pt]")
    }
    tex <- c(tex, "\\bottomrule", "\\end{tabular}")
    path <- file.path(results_dir, "Table2.tex")
    writeLines(tex, path)
    message("Combined n=5000 and n=10000 table saved to ", path)
  }

  if (!length(sample_sizes) || !is.numeric(sample_sizes) ||
      any(!is.finite(sample_sizes)) || any(sample_sizes != floor(sample_sizes))) {
    stop("AIV.summary.n must contain integer sample sizes.", call. = FALSE)
  }
  sample_sizes <- unique(sample_sizes)
  tables <- setNames(lapply(sample_sizes, export_table), as.character(sample_sizes))
  if (all(c("5000", "10000") %in% names(tables))) export_combined_table(tables)
})
