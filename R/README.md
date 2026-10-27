# R code — data construction, evaluation, tables and figures

Run scripts from this directory. Paths are relative: raw/processed data in
`../Data/`, Ox forecast output in `../Ox/<INDEX>_prediction[...]/`, and
intermediate results in `./prediction_results_2026/` (git-ignored).
Required packages: `MCS`, `sandwich`, `tidyverse` (`readr`, `ggplot2`),
`reshape2`, `readxl`, `ggrepel`, `ggnewscale`.

Suffix convention: `*_2026` = main-sample (revision) scripts;
`*_n225_2020` = COVID-period (N225) scripts.

## 1. Data construction
- `construct_data_djia.R`, `construct_data_n225.R` — build the processed
  return + realized-measure CSVs (`*_rcc_rv5_rkth2_bv_medrv*.csv`).
- `construct_rv_hl_djia.R`, `construct_rv_hl_n225.R`,
  `construct_rv_hl_n225_2020.R` — rolling Hansen-Lunde adjustment of the evaluation proxies
  (`*_hl_roll*`). The main-sample scripts were corrected on 2026-09-28 (same estimation window for
  returns and realized measures; proxy of the forecast target day) and write `*_r3corrected.csv`;
  the versions used before the correction are kept in `previous_versions/` for reference.
- `read_data.R` — shared data-loading helper.

## 2. Merge Ox forecast output → `prediction_results_2026/`
- `merge_prediction_results_2026.R` — main sample.
- `merge_prediction_results_n225_2020.R`,
  `merge_prediction_results_fs_n225_2020.R` — COVID period.

## 3. Forecast evaluation
- `examine_forecast_vol_2026.R` — MSE / QLIKE, Giacomini–White tests, and the
  Model Confidence Set (Hansen, Lunde & Nason 2011).
- `examine_forecast_var_es_2026.R` — VaR/ES: violation rate, Dynamic Quantile
  test, FZ0 loss, GW tests, MCS.
- `gwtest_mcs_vol.R` — standalone GW + MCS for volatility.
- `*_n225_2020.R` — COVID-period counterparts.

The MCS is reported at **both the 75% and 90% levels** (`**` = 75% MCS,
`*` = 90% MCS); the MCS p-value is level-independent so both sets follow from
the same p-value.

## 4. Tables and figures
- Tables: `summarize_forecast_vol_2026.R`, `summarize_prediction_results_2026.R`,
  `output_forecast_vol_2026.R`, `output_forecast_vol_gwtest_2026.R`,
  `output_forecast_var_es_gwtest_2026.R` (+ `_n225_2020` variants).
- Figures: `plot_forecast_vol_2026.R`, `plot_var_es_2026.R` (+ `_n225_2020` variants);
  `fig_dens_analytic_2026.R` (Figure 1 of the article and the density figures of the
  Supplementary Material, analytic densities from `function_dens_skew_t.R`).

## Shared helpers (sourced by the above)
- `function_forecast_evaluation.R` — VaR/ES losses (Clements et al. 2008; FZ0,
  Patton et al. 2019).
- `function_plot_forecast.R`, `function_plot_forecast_2026.R` — plotting helpers.

## Environment variables (IJF R3 pipeline)
The scripts reproduce the tables of the article with the following settings (all optional; unset
reproduces the earlier behaviour and keeps old outputs untouched):
- `RSVST_HL=r3corrected` — use the corrected Hansen-Lunde proxies (`Data/*_r3corrected.csv`) and
  write to `prediction_results_2026_r3corrected/`, `figure2026_r3corrected/`. **Required for the
  tables of the article.**
- `RSVST_AZ=dualfixed` — read the recomputed Azzalini forecasts from the aggregates produced by
  `../Ox/az_recomputation/tools/azdual_aggregate.py` (directory `RSVST_AZ_DIR`, default
  `../Ox/az_recomputation/aggregated/`) and append `_dualfixed` to the output directories.
  **Required for the tables of the article.** If the Azzalini forecasts are regenerated directly with
  the corrected classes in `../Ox/`, leave it unset.
- `RSVST_POINT=median` — evaluate the predictive median instead of the mean (not used in the article).
- `RSVST_MASK_DIR=<dir>` — drop the forecast dates listed there for all models (exclusion
  sensitivity of Section 5.2 of the Supplementary Material).

