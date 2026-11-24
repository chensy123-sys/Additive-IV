# Additive IV

R code for additive instrumental-variable estimation with prespecified and adaptive weighting functions.

## Setup

Install R and the required packages:

```bash
Rscript -e 'install.packages(c("dplyr", "tibble", "mgcv", "simcausal", "purrr", "readr", "tidyr", "ggplot2"), repos = "https://cloud.r-project.org")'
```

Start from the repository root. Both simulation runners use 500 repetitions and five-fold cross-fitting by default.

## Single-stage simulations

```bash
cd AIV_single_stage
Rscript script/run_AIV_mean_simulation.R 1 5000 all
Rscript script/run_AIV_mean_simulation.R 5 5000 all
```

Arguments are `<DGP> <sample_size> <target>`:

| DGP | Target |
| --- | --- |
| `1`–`4` | `ATE`, `treat` (mean under treatment), `control` (mean under control), or `all` |
| `5`–`6` | `-1`, `0`, `1` (mean under that treatment level), or `all` |

Multiple targets can be separated by spaces or commas. Sample size must be at least 5. Results are saved in `results/data(n)/`.

Generate Tables 2 and 3 from the saved results:

```bash
Rscript script/Table2.R
Rscript script/Table3.R
```

These summaries require results for all corresponding DGPs and targets at both `n = 5000` and `n = 10000`; these result files are included. LaTeX tables are saved in `results/`.

## Longitudinal simulations

From the repository root:

```bash
cd "AIV longitudinal"
for actions in 0,0 0,1 1,0 1,1; do
  Rscript script/run_simulation.R 5000 "$actions"
done
Rscript script/Table4.R 5000 10000 > results/Table4.tex
```

The simulation arguments are `<sample_size> <treatment_sequence>`. For example, `0,1` assigns treatment 0 at time 0 and treatment 1 at time 1. Repeat with `10000` to regenerate results at that sample size. Results are saved in `results/data(n)/` as `RES_*.RData` (prespecified weights) and `RES_adap_*.RData` (adaptive weights).

## Data application

Run from `AIV longitudinal`:

```bash
Rscript script/JPTA.R
```

This uses the included `results/jpta_han.csv` to evaluate static and dynamic treatment regimes, saving `results/JPTA.RData` and `results/Table5.tex`. To rebuild the CSV from `results/jtpa_han.tab`, run `Rscript script/data_change.R` first.

## Other files and background runs

- `AIV_single_stage/AIV.R` and `AIV longitudinal/LIV.R` contain the estimators and simulation functions.
- Each `main.R` provides a single-dataset example, runnable with `Rscript main.R` from its own directory.
- `AIV longitudinal/script/par.R` defines the longitudinal simulation model and is loaded automatically.

For a background run, execute the following from `AIV_single_stage`:

```bash
mkdir -p logs
nohup Rscript script/run_AIV_mean_simulation.R 1 5000 all > logs/dgp1_5000.log 2>&1 &
tail -f logs/dgp1_5000.log
```

Rerunning a simulation with the same settings overwrites its result files. Use `--help` with either simulation runner for command-line usage.

