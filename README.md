# IJF_Revision_Ox

Public, reproducible Ox code for the paper:

> **Realized Stochastic Volatility Model with Skew‑t Distributions for Improved Volatility and Quantile Forecasting**  
(<https://arxiv.org/abs/2401.13179>, 2026)

## What's here

- **Fully runnable synthetic-data pipeline** (no proprietary data required)
- Ox source code for RSV and EGARCH‑GED benchmarks
- Forecast evaluation (incl. Diebold–Mariano 1995) and figure/table exporters
- Clean repo structure ready for GitHub → Zenodo DOI

## Requirements

- Ox / `oxl` (e.g., Ox 9.x). Confirm with:
  ```
  oxl -v
  ```
- macOS / Linux shell (Windows WSL is fine)

## Quick start (synthetic data)

```bash
bash scripts/run_synth.sh
# or invoke directly:
oxl ox/run_estimation.ox --config=config/config.ini.example --mode=synthetic
```

Outputs will be written to `results/`.

## Real-data run (for licensed users)

The empirical high-frequency data (e.g., Oxford Realized Library) are **not included**.
If you have a legitimate license, place CSV files in a local directory and run:

```bash
oxl ox/run_estimation.ox --config=config/config.ini --mode=real --data_dir=/path/to/realdata
```

Provide your own `config/config.ini` (copy from `config/config.ini.example`).

## Data & Code Availability Statement (journal-ready)

The high-frequency data used in the empirical analysis are obtained under a paid license that **prohibits redistribution**.
Consistent with journal policy, we release **all Ox code** and a **synthetic-data pipeline** that reproduces the full
estimation and forecasting workflow. Users with legitimate access to the licensed data can re-run the real-data analysis
by supplying `--mode=real --data_dir=...`.

## Repository layout

```
IJF_Revision_Ox/
├─ README.md
├─ LICENSE
├─ .gitignore
├─ config/
│   └─ config.ini.example
├─ ox/
│   ├─ run_estimation.ox
│   ├─ simulate_data.ox
│   ├─ models.oxh
│   ├─ models_impl.ox
│   ├─ postproc.ox
│   └─ utils.ox
├─ scripts/
│   ├─ run_synth.sh
│   └─ run_real_template.sh
├─ results/
│   └─ .gitkeep
└─ data_description.txt
```

## Citation

When you cite the code, please reference the GitHub release and (if enabled) the Zenodo DOI.

---
Maintainer: Makoto Takahashi
