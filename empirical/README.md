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
| `RSV_skewFS.ox` / `RSV_skewFS.h` | RSV/SV model with Fernández-Steel skew distributions (FS-N/SN/T/ST) |
| `ReportMCMC.ox` / `ReportMCMC.h` | MCMC diagnostics and reporting utility (used by RSV_skewFS) |

Each class implements Bayesian MCMC estimation (Gibbs + Metropolis-Hastings) with 1-step-ahead rolling forecasts for volatility, VaR, and ES.

#### RSV_skewFS model variants

`RSV_skewFS` selects the error distribution via three flags set before calling `MCMC`:

| Method | Default | Effect |
|--------|---------|--------|
| `SetIsRSV(1)` | 1 (on) | Include realized variance measurement equation (RSV); set 0 for plain SV |
| `SetIsT(1)`   | 1 (on) | Use Student-t mixing (FS-ST); set 0 for normal (FS-SN) |
| `SetIsSkew(1)`| 1 (on) | Enable skewness (γ ≠ 1); set 0 for symmetric (FS-T or FS-N) |

The skewness parameter is γ (gamma). Values γ < 1 produce left-skewed returns; γ = 1 gives the symmetric t (or normal) special case.

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

### AZ / GH models

The following settings are hardcoded in each script and can be adjusted:

| Setting | Default | Variable |
|---------|---------|----------|
| Iterations (first window) | 15,000 | `vburn[0] = 5000` |
| Iterations (subsequent) | 15,000 | `vburn[1:] = 1500` |
| Rolling windows | 603 (DJIA/N225) | `nprd` |
| Multi-move sampler | ON | `flmm = 1` |

### FS model (`RSV_skewFS`)

Settings are passed directly to `MCMC(nsim, nburn)`. Typical values used in estimation:

| Setting | Typical value |
|---------|---------------|
| `nsim` (kept draws) | 90,000 |
| `nburn` (burn-in) | 60,000 |
| Multi-move sampler | Not implemented (single-move only) |

---

## Prior Specification

### AZ models

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

### FS model

| Parameter | Prior | Hyperparameters |
|-----------|-------|-----------------|
| μ | N(0, σ²) | `dm0_mu=0, dsig0_mu` (set large for flat) |
| (φ+1)/2 | Beta(a, b) | `da0_phi, db0_phi` |
| (ρ+1)/2 | Beta(a, b) | `da0_rho, db0_rho` |
| σ_η² | InvGamma(n/2, S/2) | `dn0_eta, dS0_eta` |
| γ (skewness) | Gamma(n, S) | `dn0_gamma, dS0_gamma` |
| ν (d.o.f.) | Gamma(n, S) | `dn0_nu, dS0_nu` |
| ξ | N(0, σ²) | `dm0_xi, dsig0_xi` |
| σ_u² | InvGamma(n/2, S/2) | `dn0_u, dS0_u` |

All priors are passed via `SetPrior(...)`. See `RSV_skewFS_run1.ox` for reference values.

---

## Smoke Test Results

The following results verify correct operation of the refactored code on simulated data (T=600, nsim=500, nburn=200, `ranseed(1234567)`).

### AZ / GH models

| Model | φ (0.95) | σ (0.20) | ρ (−0.30) | δ/β (−0.70/−0.50) | ν (10) | σ_u (0.30) |
|-------|----------|----------|-----------|-------------------|--------|------------|
| RSV-AZ-N  | 0.949 | 0.199 | −0.227 | — | — | 0.325 |
| RSV-AZ-ST | 0.972 | 0.205 | −0.547 | −0.881 | 7.73 | 0.309 |
| RSV-GH-ST | 0.950 | 0.217 | −0.342 | −0.672 | 10.87 | 0.312 |
| SV-AZ-ST  | 0.902 | 0.230 | −0.230 | −0.333 | 5.62 | — |
| SV-GH-ST  | 0.928 | 0.232 | −0.228 | −0.727 | 9.36 | — |

### FS model

| Parameter | True | Posterior mean | 95% CI |
|-----------|------|---------------|--------|
| φ | 0.95 | 0.953 | [0.926, 0.978] |
| σ_η | 0.20 | 0.212 | [0.182, 0.247] |
| ρ | −0.30 | −0.328 | [−0.475, −0.185] |
| γ | 0.80 | 0.809 | [0.723, 0.899] |
| ν | 10.0 | 11.3 | [6.6, 22.2] |
| ξ | −0.50 | −0.518 | [−0.629, −0.420] |
| σ_u | 0.30 | 0.324 | [0.297, 0.350] |

All true values are within their 95% posterior credible intervals.
