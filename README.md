# RSVST: Realized Stochastic Volatility Model with Skew-t Distributions

**Companion code for the arXiv preprint:**  
> *Realized Stochastic Volatility Model with Skew-t Distributions for Improved Volatility and Quantile Forecasting*  
> Makoto Takahashi, Wataru Omori, Yusuke Watanabe, and Shinya Yamauchi  
> [arXiv:2501.13179](https://arxiv.org/abs/2401.13179)

---

## Overview

This repository provides **Ox source code** and a **synthetic-data simulation pipeline** for reproducing the results in the above paper.  
The project implements **Realized Stochastic Volatility (RSV)** and **EGARCH-GED** models for volatility and quantile forecasting, with extensions to the **skew-t** distribution.

The empirical high-frequency data (Oxford Realized Library) are **not included** due to licensing restrictions, but all estimation, simulation, and evaluation procedures can be replicated using the provided synthetic-data routines.

---

## Contents

```
RSVST/
├─ ox/
│  ├─ run_estimation.ox         # Main entry point
│  ├─ simulate_data.ox          # Synthetic data generator
│  ├─ models_impl.ox            # RSV & EGARCH-GED estimation (skeleton)
│  ├─ postproc.ox               # Forecast evaluation and export
│  └─ utils.ox                  # CLI parser and helpers
├─ config/config.ini.example    # Example configuration file
├─ scripts/run_synth.sh         # One-command synthetic run
├─ results/                     # Output directory (empty by default)
├─ data_description.txt         # Variables and data policy
└─ README.md                    # This file
```

---

## Requirements

- **Ox / `oxl`** version 9.x or later  
  Check installation with:
  ```bash
  oxl -v
  ```
- macOS / Linux / WSL (Windows Subsystem for Linux)

---

## Quick start (synthetic run)

```bash
bash scripts/run_synth.sh
# or
oxl ox/run_estimation.ox --config=config/config.ini.example --mode=synthetic
```

Outputs (tables, figures) will be generated in `results/`.

---

## Real-data replication (licensed users)

Researchers with legitimate access to the **Oxford Realized Library** or equivalent high-frequency data can replicate the empirical results as follows:

```bash
cp config/config.ini.example config/config.ini
# Edit config.ini to specify your data directory
oxl ox/run_estimation.ox --config=config/config.ini --mode=real --data_dir=/path/to/realdata
```

Expected variables:
- `date` (YYYYMMDD)
- `rcc` – daily return
- `rv` – 5-min realized variance
- `rk` – realized kernel (Parzen)

---

## Data and Code Availability Statement

All **Ox codes** and **synthetic data pipelines** are publicly released under the MIT License.  
The empirical data are obtained under a paid license (Oxford Realized Library) and cannot be redistributed.  
Users with legitimate access can reproduce the real-data analysis by specifying `--mode=real` and providing their local data path.

---

## Citation

If you use this code or derivative work, please cite:

```
Takahashi, M., Omori, W., Watanabe, Y., & Yamauchi, S. (2025).
Realized Stochastic Volatility Model with Skew-t Distributions for Improved Volatility and Quantile Forecasting.
arXiv preprint arXiv:2501.xxxxx
```

BibTeX:
```bibtex
@article{takahashi2025rsvst,
  title={Realized Stochastic Volatility Model with Skew-t Distributions for Improved Volatility and Quantile Forecasting},
  author={Takahashi, Makoto and Omori, Wataru and Watanabe, Yusuke and Yamauchi, Shinya},
  journal={arXiv preprint arXiv:2501.xxxxx},
  year={2025}
}
```

---

## License

This repository is released under the **MIT License**.  
You are free to use, modify, and redistribute the code, provided that appropriate credit is given.

---

## Contact

Maintainer: **Makoto Takahashi**  
Faculty of Business Administration, Hosei University  
GitHub: [mkttkhsh](https://github.com/mkttkhsh)
