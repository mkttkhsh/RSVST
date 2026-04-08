---
layout: home
title: RSVST
---

# Realized Stochastic Volatility Model with Skew-t Distributions

**Companion code** for the paper:

> Takahashi, M., Yamauchi, Y., Watanabe, T., & Omori, Y. (2025).  
> *Realized Stochastic Volatility Model with Skew-t Distributions for Improved Volatility and Quantile Forecasting.*  
> [arXiv:2501.13179](https://arxiv.org/abs/2501.13179)

---

## Overview

Realized stochastic volatility (RSV) models augment the standard SV framework with a **measurement equation** for a realized variance proxy (e.g., 5-min RV or a realized kernel):

$$
y_t = \mu + \exp(h_t / 2)\,\varepsilon_t, \qquad
\log\mathrm{RV}_t = h_t + \xi + \sigma_u\,u_t
$$

$$
h_{t+1} = \phi\,h_t + \sigma_\eta\,\eta_t, \qquad \rho = \mathrm{Corr}(\varepsilon_t, \eta_t)
$$

where $\varepsilon_t$ follows one of several skew-t distributions. This repository implements:

| Model | Distribution | Type |
|-------|-------------|------|
| RSV-AZ-N | Normal | RSV |
| RSV-AZ-SN | Azzalini skew-normal | RSV |
| RSV-AZ-T | Student-*t* | RSV |
| RSV-AZ-ST | Azzalini skew-*t* | RSV |
| RSV-GH-ST | Generalized Hyperbolic skew-*t* | RSV |
| EGARCH-GED | Generalized Error Distribution | GARCH |

Each model is estimated by **Bayesian MCMC** (Gibbs + Metropolis-Hastings), with **1-step-ahead rolling forecasts** for volatility, Value-at-Risk, and Expected Shortfall.

---

## Repository Contents

```
RSVST/
├── ox/
│   ├── run_estimation.ox     # Main entry point
│   ├── simulate_data.ox      # Synthetic data generator
│   ├── models_impl.ox        # RSV & EGARCH-GED estimation
│   ├── models.oxh            # Class declarations
│   ├── postproc.ox           # Forecast evaluation and export
│   └── utils.ox              # CLI parser and helpers
├── config/
│   └── config.ini.example    # Example configuration
├── scripts/
│   ├── run_synth.sh          # One-command synthetic run
│   └── run_real_template.sh  # Template for real-data run
├── results/                  # Output directory (empty by default)
├── data_description.txt      # Expected data format
├── LICENSE                   # MIT License
└── README.md
```

---

## Quick Start

### Requirements

- **Ox 9.x or later** — verify with `oxl -v`
- macOS / Linux / WSL

### Synthetic-data run

```bash
git clone https://github.com/mkttkhsh/RSVST.git
cd RSVST
bash scripts/run_synth.sh
```

Outputs (tables, figures) appear in `results/`.

### Real-data run (licensed users)

```bash
cp config/config.ini.example config/config.ini
# Edit config.ini to set data_dir
oxl ox/run_estimation.ox --config=config/config.ini --mode=real
```

Expected CSV columns: `date` (YYYYMMDD), `rcc`, `rv`, `rk`.

---

## Configuration

`config/config.ini.example`:

```ini
[run]
mode    = synthetic   # synthetic | real
seed    = 1234567
outdir  = results

[models]
use_rsv      = 1
use_egarch_ged = 1
horizon      = 1      # forecast horizon (days)

[simulation]
T = 4000              # sample size for synthetic run
```

---

## Data Policy

Empirical high-frequency data (Oxford Realized Library) are **not included** due to licensing restrictions.  
The full synthetic-data pipeline can be run without any proprietary data.  
Researchers with licensed access can reproduce the empirical results by setting `--mode=real`.

See [`data_description.txt`](https://github.com/mkttkhsh/RSVST/blob/main/data_description.txt) for the expected variable format.

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

Released under the [MIT License](https://github.com/mkttkhsh/RSVST/blob/main/LICENSE).
