# Data

The raw and processed data files are **not redistributed** in this repository
because the underlying index returns and realized measures are subject to the
terms of their original providers. This file documents (i) which files the code
expects, (ii) where to obtain the raw inputs, and (iii) how to rebuild the
processed files.

> Place the processed CSVs listed below in this `Data/` directory (the Ox and R
> scripts load them via the relative path `../Data/...`). The Fernández–Steel
> code additionally expects two FS CSVs in `Data/` (`DJIA_rcc_rv_2009_2019.csv`, `Nikkei_rcc_rv_2009_2019.csv`).

## 1. Processed files expected by the code

All CSVs are daily and (unless noted) have the layout
`Date(yyyymmdd), close-to-close return, realized measure(s)`.

| File (in `Data/`)                                   | Used by                                   | Columns |
|-----------------------------------------------------|-------------------------------------------|---------|
| `DJIA_rcc_rv_200906_201909.csv`                     | `Ox/rsv_*_djia_prediction.ox`         | Date, return, RV(5-min) |
| `N225_rcc_rv_200906_201909.csv`                     | `Ox/rsv_*_n225_prediction.ox`         | Date, return, RV(5-min) |
| `DJIA_rcc_rv5_rkth2_bv_medrv.csv`                   | volatility-proxy evaluation               | Date, return, RV5, RKth2, BV, MedRV |
| `N225_rcc_rv5_rkth2_bv_medrv.csv`                   | volatility-proxy evaluation               | Date, return, RV5, RKth2, BV, MedRV |
| `N225_rcc_rv5_rkth2_bv_medrv_2012_2020.csv`         | COVID-period (`*_n225_prediction_2020.ox`)| Date, return, RV5, RKth2, BV, MedRV |
| `DJIA_rcc_rv5_rkth2_bv_medrv_2012_2020_tick.csv`    | COVID-period (`*_djia_prediction_2020.ox`)| Date, return, RV5, RKth2, BV, MedRV |
| `N225_rcc_rv5_rkth2_bv_medrv_hl_roll_2020.csv`      | COVID-period evaluation (`merge_prediction_results_n225_2020.R`) | Date, return, RV5/RKth2/BV/MedRV (+ HL-corrected) |
| `DJIA_rcc_rv5_rkth2_bv_medrv_2012_2020_hl_roll_tick.csv` | COVID-period evaluation (`merge_prediction_results_djia_2020.R`) | Date, return, RV5/RKth2/BV/MedRV (+ HL-corrected) |
| `Data/DJIA_rcc_rv_2009_2019.csv` (FS)      | `RSV_skewFS_run*.ox`                       | Date, return, RV |
| `Data/Nikkei_rcc_rv_2009_2019.csv` (FS)    | `RSV_skewFS_run*.ox`                       | Date, return, RV |

> **`_tick` suffix (DJIA only, 2026 revision):** the DJIA realized measures were
> recomputed from raw tick data to match the N225 construction convention
> (RV5/BV=5-min, FTRK/MedRV=1-min; see `PROGRESS.md`). Files carrying the
> `_tick` suffix are this recomputed version, kept alongside the original
> Oxford-Man-Institute-based files during the revision for validation; the
> suffix will be dropped in a final consolidation pass once the co-authors
> confirm the MedRV choice.

(Returns are multiplied by 100 and realized variances by 100^2 inside the
scripts, so store them in raw/decimal units as produced by the build scripts.)

## 2. Raw data sources

As stated in Section 4.1 of the article:

- **DJIA.** Daily returns and realized measures (5-minute realized variance, realized kernel,
  bipower variation, median realized variance) are from the Oxford-Man Institute's *Realized
  Library* (Heber, Lunde, Shephard and Sheppard, 2009). The library has been discontinued and its
  website (https://realized.oxford-man.ox.ac.uk/) is no longer accessible; no official archive
  exists. The file is referenced locally as `oxfordmanrealizedvolatilityindices.csv`. The DJIA
  COVID-period files with the `_tick` suffix were recomputed from tick data with the same
  construction as the N225 series.
- **N225.** Daily returns and realized measures are constructed from the Nikkei NEEDS-TICK dataset
  (Nikkei Media Marketing), a commercial dataset; see Ubukata and Watanabe (2014) for the
  construction. The realized measures are referenced locally as an Excel workbook with sheets `RV`,
  `FTRK`, `BPV`, `MedRV` (`Nikkei_RV_1996_202509.xlsx`).

Neither the raw inputs nor the processed CSVs derived from them are redistributed here, in
accordance with the providers' terms. Researchers with access to the same sources can rebuild the
processed files with the scripts in section 3; the Ox and R code runs unchanged on any file with
the column layout of section 1.

## 3. Rebuilding the processed files

The `R/construct_*.R` scripts turn the raw inputs into the processed CSVs:

- `R/construct_data_djia.R`, `R/construct_data_n225.R` — build the per-index
  `*_rcc_rv5_rkth2_bv_medrv*.csv` (and the COVID-period `*_2012_2020.csv`) from
  the raw return + realized-measure inputs.
- `R/construct_rv_hl_djia.R`, `R/construct_rv_hl_n225.R`,
  `R/construct_rv_hl_n225_2020.R` — append the high–low range / Rogers–Satchell
  measures (`*_hl_roll*.csv`).
- The DJIA `_tick` files (main-sample and COVID-period) are currently built by
  scratch scripts (`compute_rv_djia_1min.R`, `construct_rv_hl_djia_tick.R`,
  `construct_rv_hl_djia_2020_tick.R`) that have not yet been migrated into this
  package; see `PROGRESS.md` in the revision working directory for status.

Run them from the `R/` directory (they read/write via `../Data/`). Adjust the
input paths at the top of each script to point at your downloaded raw files.
