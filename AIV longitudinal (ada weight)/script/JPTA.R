rm(list=ls())
suppressMessages({source('LIV.R')})
data_JPTA <- read.csv('results/jpta_han.csv')[,2:8] %>%
  mutate(id = row_number()) %>%
  pivot_longer(
    cols = c(A0, Z0, L0, A1, Z1, L1, Y2),
    names_to = c(".value", "time"),
    names_pattern = "([A-Z]+)([0-9]+|Y)?"
  ) %>%
  arrange(id, time) %>% rename(ID=id, t=time) %>%
  type_convert()
print(data_JPTA)

set.seed(2025); smoothness = 8
act_0p <- list(function(fY) {return(0)},function(fY){as.numeric(fY$L<=0.8368623)}) # 0.8 quantile
act_0n <- list(function(fY) {return(0)},function(fY){as.numeric(fY$L>0.8368623)})
act_1p <- list(function(fY) {return(1)},function(fY){as.numeric(fY$L<=0.8368623)})
act_1n <- list(function(fY) {return(1)},function(fY){as.numeric(fY$L>0.8368623)})

fit_LIV_adap_00 <- DML_LIV_adap(data_JPTA, action_list = c(0, 0), K = 1, H = 2,smoothness = smoothness)
fit_LIV_adap_01 <- DML_LIV_adap(data_JPTA, action_list = c(0, 1), K = 1, H = 2,smoothness = smoothness)
fit_LIV_adap_0p <- DML_LIV_adap(data_JPTA, action_list = act_0p,  K = 1, H = 2,smoothness = smoothness)
fit_LIV_adap_0n <- DML_LIV_adap(data_JPTA, action_list = act_0n, K = 1, H = 2,smoothness = smoothness)

fit_LIV_adap_10 <- DML_LIV_adap(data_JPTA, action_list = c(1, 0), K = 1, H = 2,smoothness = smoothness)
fit_LIV_adap_11 <- DML_LIV_adap(data_JPTA, action_list = c(1, 1), K = 1, H = 2,smoothness = smoothness)
fit_LIV_adap_1p <- DML_LIV_adap(data_JPTA, action_list = act_1p, K = 1, H = 2,smoothness = smoothness)
fit_LIV_adap_1n <- DML_LIV_adap(data_JPTA, action_list = act_1n, K = 1, H = 2,smoothness = smoothness)
save.image("results/JPTA.RData")

res <- rbind(
  (predict_LIV(fit_LIV_adap_00) %>% select(est,std))[1,],
  (predict_LIV(fit_LIV_adap_01) %>% select(est,std))[1,],
  (predict_LIV(fit_LIV_adap_0p) %>% select(est,std))[1,],
  (predict_LIV(fit_LIV_adap_0n) %>% select(est,std))[1,],
  
  (predict_LIV(fit_LIV_adap_10) %>% select(est,std))[1,],
  (predict_LIV(fit_LIV_adap_11) %>% select(est,std))[1,],
  (predict_LIV(fit_LIV_adap_1p) %>% select(est,std))[1,],
  (predict_LIV(fit_LIV_adap_1n) %>% select(est,std))[1,],
  
  (predict_LIV(fit_LIV_adap_00) %>% select(est,std))[2,],
  (predict_LIV(fit_LIV_adap_01) %>% select(est,std))[2,],
  (predict_LIV(fit_LIV_adap_0p) %>% select(est,std))[2,],
  (predict_LIV(fit_LIV_adap_0n) %>% select(est,std))[2,]
) %>% round(.,3) %>% data.frame() %>% 
  mutate(DTR=c('00','01','0d+','0d-','10','11','1d+','1d-','x0','x1','xd+','xd-'))
print(res)

# Export the 12 regimes in the same order as the manuscript table.
latex_headers <- c(
  "\\texttt{00}", "\\texttt{01}", "\\texttt{0$d^+$}", "\\texttt{0$d^-$}",
  "\\texttt{10}", "\\texttt{11}", "\\texttt{1$d^+$}", "\\texttt{1$d^-$}",
  "\\texttt{x0}", "\\texttt{x1}", "\\texttt{x$d^+$}", "\\texttt{x$d^-$}"
)
stopifnot(
  identical(res$DTR, c("00", "01", "0d+", "0d-", "10", "11",
                       "1d+", "1d-", "x0", "x1", "xd+", "xd-")),
  all(is.finite(res$est)),
  all(is.finite(res$std)),
  all(res$std >= 0)
)

latex_rows <- vapply(
  list(DTRs = latex_headers,
       EST = sub("^(-?)0\\.", "\\1.", sprintf("%.3f", res$est)),
       SE = sub("^(-?)0\\.", "\\1.", sprintf("%.3f", res$std))),
  function(values) paste(values, collapse = " & "),
  character(1)
)
latex_table <- c(
  "\\begin{tabular}{c| cccc| cccc| cccc}",
  "\\toprule",
  paste0("DTRs & ", latex_rows[["DTRs"]], " \\\\"),
  "\\midrule",
  paste0("EST & ", latex_rows[["EST"]], " \\\\"),
  paste0("SE & ", latex_rows[["SE"]], " \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)
latex_file <- file.path("results", "Table5.tex")
writeLines(latex_table, con = latex_file, useBytes = TRUE)
message("LaTeX table saved to ", normalizePath(latex_file, mustWork = TRUE))
