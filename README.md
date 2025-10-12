# Additive-IV
## AIV_single_stage

In this folder, a simulation study is designed to estimate the mean of potential outcomes using additive instrumental variables (AIVs) in cross-sectional data. 

The implementation of the cross-fitting procedure with **prespecified weights** and **adaptive weights** can be found in the script `AIV.R`. 

One can estimate the mean potential outcome $\mathbb E[Y(0)]$ using `control1.R`, ..., `control4.R`, estimating $\mathbb E[Y(1)]$ using `treat1.R`, ..., `treat4.R`, $\mathbb E[Y(1)-Y(0)]$ using `ATE1.R`, ..., `ATE4.R`. The results can be seen by running `summary_ATE(5000).R`, `summary_control(5000).R`, `summary_treat(5000).R`. Parallel results for a sample size of (n = 2000) are available in `Summary(2000).R`.

## AIV Longitudinal

In this folder, a simulation study is designed to estimate the mean of potential outcomes using additive instrumental variables (AIVs) in longitudinal data. 

The implementation of the backward cross-fitting procedure can be found in the script `LIV_overall.R`. Additionally, one can run the R-scripts `Overall00.R`, `Overall01.R`, ..., `Overall11.R` to conduct the simulations for estimating  $\mathbb E[Y(0,0)],\mathbb E[Y(0,1)],\ldots,\mathbb E[Y(1,1)]$, with the hyperparameters recorded in `par_overall1.R`. 

The results for a sample size of (n = 5000) can be found in `Overall_00.RData`, ..., `Overall_11.RData`. For detailed summaries, you can use the R-script `Summary(5000).R`. Parallel results for a sample size of (n = 2000) are available in `Summary(2000).R`.

## JTPA

In this folder, we apply our method to study returns to schooling and post-school training as sequential treatments and to conduct a policy analysis, using the dataset `jtpa.tab` provided in the online supplementary materials of [Full article: Optimal Dynamic Treatment Regimes and Partial Welfare Ordering](https://www.tandfonline.com/doi/full/10.1080/01621459.2023.2238941#supplemental-material-section). 

The `R` script `data_change.R` preprocesses the `jtpa.tab` to an applicable form `jtpa.csv`.  The implementation of the backward cross-fitting procedure can be found in the script `LIV_overall.R`. The potential terminal income when 8 different types of dynamic treatment rules are designed are estimated by `overall_00.R`, ..., `overall_13.R`. 

The results are reported by running `summary.R`.
