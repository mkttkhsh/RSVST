# Recomputation of the Azzalini forecasts (dual-output classes)

Companion to Section 5.2 of the Supplementary Material of the October 2026 revision. See `../README.md`, section "Recomputation of the Azzalini forecasts".

- `RSVAZ.ox`, `SVAZ.ox`, `RSVAZ.h`, `SVAZ.h`, `ReportMCMC.h`: the model classes with the corrected
  one-day-ahead forecasting step (`sqrt(1-delta^2)`), extended to also compute the forecast under
  the original formula (`sqrt(1-c^2)`) from the same draws, and to save the per-draw state.
- `azdual_*.ox`: drivers, one per (sample, model, chunk). Run from this directory, e.g.
  `oxl azdual_djia_main_sv_az_st_c0.ox`; the data path is `../../Data/`.
- `tools/`: aggregation, diagnostics, comparison and multi-chain confirmation scripts (Python 3,
  standard library only). Thresholds: inefficiency factor of h_n > 300 or block ratio of
  mean(exp(h_{n+1})) outside [0.5, 2] flags a window; acceptance requires, for every chain, IF < 300,
  Geweke p > 0.01 and block ratio in [0.5, 2], agreement of the chain means of h_n within half a
  posterior sd, rank-normalized split-Rhat < 1.01 with bulk and tail ESS >= 400, and Monte Carlo
  standard errors of the predictive mean of volatility and of the VaR/ES forecasts below 0.5% and 1%.
