# empirical/ — MCMC Estimation Code for Empirical Analysis

This directory contains the **Ox source code** used for the empirical results in:

> Takahashi, M., Yamauchi, Y., Watanabe, T., & Omori, Y. (2025).  
> *Realized Stochastic Volatility Model with Skew-t Distributions for Improved Volatility and Quantile Forecasting.*  
> [arXiv:2501.13179](https://arxiv.org/abs/2501.13179)

---

## File Overview

### Model Class Libraries

| File | Description |
|------|-------------|
| `RSVAZ.ox` / `RSVAZ.h` | RSV model with Azzalini skew distributions (RSV-AZ-N/SN/T/ST) |
| `RSVGH.ox` / `RSVGH.h` | RSV model with Generalized Hyperbolic distribution (RSV-GH-ST) |
| `SVAZ.ox` / `SVAZ.h`   | Standard SV (no realized variance) with Azzalini distributions |
| `SVGH.ox` / `SVGH.h`   | Standard SV with Generalized Hyperbolic distribution |

Each class implements Bayesian MCMC estimation (Gibbs + Metropolis-Hastings) with 1-step-ahead rolling forecasts for volatility, VaR, and ES.

### Prediction Scripts (rolling window, 1-step-ahead forecasts)

| Script | Model | Dataset |
|--------|-------|---------|
| `rsv_az_st_djia_prediction.ox` | RSV-AZ-ST | DJIA |
| `rsv_az_sn_djia_prediction.ox` | RSV-AZ-SN | DJIA |
| `rsv_az_n_djia_prediction.ox`  | RSV-AZ-N  | DJIA |
| `rsv_gh_st_djia_prediction.ox` | RSV-GH-ST | DJIA |
| `sv_az_st_djia_prediction.ox`  | SV-AZ-ST  | DJIA |
| `sv_az_sn_djia_prediction.ox`  | SV-AZ-SN  | DJIA |
| `rsv_az_st_n225_prediction.ox` | RSV-AZ-ST | Nikkei 225 |
| `rsv_az_sn_n225_prediction.ox` | RSV-AZ-SN | Nikkei 225 |
| `rsv_az_n_n225_prediction.ox`  | RSV-AZ-N  | Nikkei 225 |
| `rsv_gh_st_n225_prediction.ox` | RSV-GH-ST | Nikkei 225 |
| `sv_az_st_n225_prediction.ox`  | SV-AZ-ST  | Nikkei 225 |
| `sv_az_sn_n225_prediction.ox`  | SV-AZ-SN  | Nikkei 225 |
| `*_n225_prediction_2020.ox`    | Various   | Nikkei 225 (extended, 2012–2020) |

### Simulation Scripts (Monte Carlo study)

| Script | Model |
|--------|-------|
| `rsv_az_st_simulation.ox` | RSV-AZ-ST |
| `rsv_az_sn_simulation.ox` | RSV-AZ-SN |
| `rsv_az_n_simulation.ox`  | RSV-AZ-N  |
| `rsv_az_t_simulation.ox`  | RSV-AZ-T  |
| `rsv_gh_st_simulation.ox` | RSV-GH-ST |
| `rsv_gh_n_simulation.ox`  | RSV-GH-N  |

### Post-processing

| Script | Description |
|--------|-------------|
| `merge_prediction_results_az.ox` | Collects rolling forecast outputs for AZ models |
| `merge_prediction_results_gh.ox` | Collects rolling forecast outputs for GH models |
| `merge_prediction_results_az_2020.ox` | Same, for 2020 extended sample |
| `merge_prediction_results_gh_2020.ox` | Same |

---

## Requirements

- **Ox 9.x or later** — verify with `oxl -v`
- macOS / Linux / WSL

---

## Data Setup

Empirical data (Oxford Realized Library) are **not included** due to licensing restrictions.  
Place your data files in `../data/` (i.e., `RSVST/data/`), or update the `sdata_dir` variable at the top of each prediction script to point to your data directory.

### Expected data format

**For DJIA / Nikkei 225 (2009–2019 sample):**  
Filename pattern: `{MARKET}_rcc_rv_200906_201909.csv`

| Column | Description |
|--------|-------------|
| col 0 | Date (YYYYMMDD) |
| col 1 | Close-to-close return (decimal, e.g. 0.012) |
| col 2 | 5-min realized variance (decimal, e.g. 0.0002) |

**For extended N225 (2012–2020 sample):**  
Filename: `N225_rcc_rv5_rkth2_bv_medrv_2012_2020.csv`

| Column | Description |
|--------|-------------|
| col 0 | Date (YYYYMMDD) |
| col 1 | `rcc` — close-to-close return |
| col 2 | `rv5` — 5-min realized variance |
| col 3 | `rk_th2` — realized kernel (Parzen) |
| col 4 | `bv` — bipower variation |
| col 5 | `medrv` — median RV |

---

## Running a Prediction Script

```bash
cd RSVST/empirical

# 1. Edit the data directory in the script (if not using ../data/)
#    Open rsv_az_st_djia_prediction.ox and set:
#      sdata_dir = "/absolute/path/to/your/data/";

# 2. Run
oxl rsv_az_st_djia_prediction.ox
```

Output CSV files are written to a subdirectory created automatically  
(e.g., `DJIA_prediction/` for DJIA scripts).

---

## Running a Simulation Script

```bash
cd RSVST/empirical
oxl rsv_az_st_simulation.ox
```

No data files required — all data are generated synthetically.

---

## MCMC Settings

The following settings are hardcoded in each script and can be adjusted:

| Setting | Default | Variable |
|---------|---------|----------|
| Iterations (first window) | 15,000 | `vburn[0] = 5000` |
| Iterations (subsequent) | 15,000 | `vburn[1:] = 1500` |
| Rolling windows | 603 (DJIA/N225) | `nprd` |
| Multi-move sampler | ON | `flmm = 1` |

---

## Prior Specification

Common priors across all AZ models:

| Parameter | Prior | Hyperparameters |
|-----------|-------|-----------------|
| μ | N(0, 10²) | `dmm0=0, dsm0=10` |
| (φ+1)/2 | Beta(20, 1.5) | `dap0=20, dbp0=1.5` |
| (ρ+1)/2 | Beta(1, 2) | `dar0=1, dbr0=2` |
| σ² | InvGamma(2.5, 0.025) | `dn0=5, dS0=0.05` |
| (δ+1)/2 | Beta(1, 2) | `dad0=1, dbd0=2` |
| ν | Gamma(5, 0.5) | `dnn0=5, dnS0=0.5` |
| ξ | N(0, 1) | `dmx0=0, dsx0=1` |
| σ_u² | InvGamma(2.5, 0.1) | `dsun0=5, dsuS0=0.2` |
