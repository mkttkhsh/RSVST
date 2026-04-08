# RSVST: Realized Stochastic Volatility Model with Skew-t Distributions

**Companion code for:**
> Takahashi, M., Yamauchi, Y., Watanabe, T., & Omori, Y. (2025).  
> *Realized Stochastic Volatility Model with Skew-t Distributions for Improved Volatility and Quantile Forecasting.*  
> [arXiv:2501.13179](https://arxiv.org/abs/2501.13179)

---

## Overview

This repository provides **Ox source code** and a **synthetic-data pipeline** for the RSV models introduced in the paper above.

The **Realized Stochastic Volatility (RSV)** framework augments the standard SV model with a measurement equation for a realized variance proxy:

```
y_t   = mu + exp(h_t / 2) * eps_t
logRV_t = h_t + xi + sigma_u * u_t
h_{t+1} = phi * h_t + sigma_eta * eta_t,   Corr(eps_t, eta_t) = rho
```

where `eps_t` follows a skew-t distribution. Models implemented:

| Model | Distribution of eps_t | Type |
|-------|-----------------------|------|
| RSV-AZ-N  | Normal | RSV |
| RSV-AZ-SN | Azzalini skew-normal | RSV |
| RSV-AZ-T  | Student-*t* | RSV |
| RSV-AZ-ST | Azzalini skew-*t* | RSV |
| RSV-GH-ST | Generalized Hyperbolic skew-*t* | RSV |
| EGARCH-GED | Generalized Error Distribution | GARCH |

Each RSV model is estimated by **Bayesian MCMC** (Gibbs + Metropolis-Hastings) with **1-step-ahead rolling forecasts** for volatility, Value-at-Risk (VaR), and Expected Shortfall (ES).

The empirical high-frequency data (Oxford Realized Library) are **not included** due to licensing restrictions. All estimation and evaluation procedures can be replicated with the provided synthetic-data routines.

---

## Repository Contents

```
RSVST/
├── empirical/                # MCMC estimation code (empirical analysis)
│   ├── RSVAZ.ox / RSVAZ.h   # RSV-Azzalini model class
│   ├── RSVGH.ox / RSVGH.h   # RSV-GH model class
│   ├── SVAZ.ox  / SVAZ.h    # SV-Azzalini model class (benchmark)
│   ├── SVGH.ox  / SVGH.h    # SV-GH model class (benchmark)
│   ├── rsv_az_st_*_prediction.ox   # Rolling forecast scripts
│   ├── rsv_gh_st_*_prediction.ox
│   ├── sv_az_st_*_prediction.ox
│   ├── *_simulation.ox              # Monte Carlo simulation scripts
│   ├── merge_prediction_results_*.ox
│   └── README.md
├── ox/                       # Synthetic-data pipeline (skeleton)
│   ├── run_estimation.ox
│   ├── simulate_data.ox
│   ├── models.oxh / models_impl.ox
│   ├── postproc.ox / utils.ox
├── config/
│   └── config.ini.example
├── scripts/
│   ├── run_synth.sh
│   └── run_real_template.sh
├── data/                     # Place data files here (not included)
├── results/                  # Output directory (empty by default)
├── data_description.txt
├── LICENSE
└── README.md
```

---

## Requirements

- **Ox 9.x or later** — verify with `oxl -v`
- macOS / Linux / WSL

---

## Quick Start

### Empirical analysis (MCMC, requires licensed data)

```bash
# Place data files in data/ (see data_description.txt)
cd RSVST/empirical
oxl rsv_az_st_djia_prediction.ox   # RSV-AZ-ST on DJIA
oxl rsv_gh_st_n225_prediction.ox   # RSV-GH-ST on Nikkei 225
```

See [`empirical/README.md`](empirical/README.md) for data format and MCMC settings.

### Synthetic-data pipeline

```bash
git clone https://github.com/mkttkhsh/RSVST.git
cd RSVST
bash scripts/run_synth.sh
```

Outputs (tables, figures) appear in `results/`.

Alternatively, run directly:

```bash
oxl ox/run_estimation.ox --config=config/config.ini.example --mode=synthetic
```

---

## Configuration

Edit `config/config.ini.example` (or copy it to `config/config.ini`):

```ini
[run]
mode    = synthetic   # synthetic | real
seed    = 1234567
outdir  = results

[models]
use_rsv       = 1
use_egarch_ged = 1
horizon       = 1     # forecast horizon (days)

[simulation]
T = 4000              # sample size for synthetic run
```

---

## Real-Data Replication (licensed users)

Researchers with access to the **Oxford Realized Library** (or equivalent) can replicate the empirical results:

```bash
cp config/config.ini.example config/config.ini
# Edit config.ini: set data_dir to your local data path
oxl ox/run_estimation.ox --config=config/config.ini --mode=real
```

Expected CSV columns (see `data_description.txt`):

| Column | Description |
|--------|-------------|
| `date` | YYYYMMDD |
| `rcc`  | Daily close-to-close log-return |
| `rv`   | 5-min realized variance |
| `rk`   | Realized kernel (Parzen) |

---

## Data and Code Availability

All **Ox source code** and the **synthetic-data pipeline** are released under the MIT License.  
The empirical data are obtained under a paid license (Oxford Realized Library) and cannot be redistributed.

---

## Citation

```bibtex
@article{takahashi2025rsvst,
  title   = {Realized Stochastic Volatility Model with Skew-t Distributions
             for Improved Volatility and Quantile Forecasting},
  author  = {Takahashi, Makoto and Yamauchi, Yuta and
             Watanabe, Toshiaki and Omori, Yasuhiro},
  journal = {arXiv preprint arXiv:2501.13179},
  year    = {2025}
}
```

---

## License

Released under the [MIT License](LICENSE).

---

## Contact

**Makoto Takahashi**  
Faculty of Business Administration, Hosei University  
GitHub: [mkttkhsh](https://github.com/mkttkhsh)
