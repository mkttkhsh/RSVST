# IJF R2 (2026): merge psi=1 (main spec) vs psi-free (sensitivity) RSV and
# REGARCH forecasts for the N225 COVID period, side by side, for the psi
# sensitivity analysis (Reviewer 1 Comment 2). Companion to
# merge_prediction_results_n225_2020.R (which builds the full
# SV/RSV/EGARCH/REGARCH comparison table for the main non-psi results):
# this one is narrower in model scope (RSV and REGARCH only, since SV/EGARCH
# have no psi) but wider in spec (psi=1 and psi-free per model, for direct
# comparison) than that script.
#
# Inputs: the psi=1 var_/es_/volf_<model>.csv files already existed (copied
# in from Ox2022/ 2026-07-11); the psi-free ones were consolidated from raw
# per-window output by merge_ox_perwindow_psi_n225.R (AZ/GH) and
# merge_prediction_results_fs_psi_n225_2020.R (FS) in this session, since no
# local Ox install exists to run the usual Ox-level merge scripts. REGARCH
# psi=1 (REGARCHn_/REGARCHt_N225_2020.xlsx) and psi-free
# (REGARCHnvarphiF_/REGARCHtvarphiF_N225_2020.xlsx) forecasts were delivered
# by Watanabe-sensei (2026-07-24); both share the same column layout as the
# main (non-psi) REGARCH files, so the extraction mirrors
# merge_prediction_results_n225_2020.R.

rm(list = ls())
library(readr)
library(readxl)

# 7 RSV distributions, each with a psi=1 and a psi-free variant
dists <- c("az_n" = "AZ-N", "az_t" = "AZ-T", "az_sn" = "AZ-SN", "az_st" = "AZ-ST",
           "fs_sn" = "FS-SN", "fs_st" = "FS-ST", "gh_st" = "GH-ST")
dist.n <- length(dists)

data.name <- "N225"
path.data <- "../Data/"
path.in <- paste0("../Ox/", data.name, "_prediction_2020/")
# 2026-10-01 (§86): RSVST_AZ -> RSV-AZ-SN/ST (psi=1 and psi-free) are read from the corrected-recomputation aggregates; AZ-N/T, FS, GH stay as before
path.in.az <- if (nzchar(Sys.getenv("RSVST_AZ", ""))) paste0(Sys.getenv("RSVST_AZ_DIR", "../Ox/az_recomputation/aggregated/"), sub("_excl$", "", Sys.getenv("RSVST_AZ")), "/", data.name, "_prediction_2020/") else path.in
p.in.az <- function(model) if (model %in% c("rsv_az_sn", "rsv_az_st", "rsv_az_sn_psi", "rsv_az_st_psi")) path.in.az else path.in
path.out <- "./prediction_results_n225_2020/"
path.regarch <- "../REGARCH_Bayes/"
true.sets <- "N225_rcc_rv5_rkth2_bv_medrv_hl_roll_2020"

# 2026-09-28 (IJF R3, PROGRESS.md §81): RSVST_HL=r3corrected switches to the corrected Hansen-Lunde proxies
# (Data/*_hl_roll*_r3corrected.csv) and to separate output directories (*_r3corrected/); old results stay untouched.
hl.tag <- Sys.getenv("RSVST_HL", "")
if (nzchar(hl.tag)) {
  if (exists("true.sets")) true.sets <- paste0(true.sets, "_", hl.tag)
  for (.p in intersect(c("path.result", "path.out", "path.figure"), ls())) {
    assign(.p, sub("/?$", paste0("_", hl.tag, "/"), get(.p))); dir.create(get(.p), showWarnings = FALSE, recursive = TRUE)
  }
  cat("RSVST_HL =", hl.tag, "\n")
}

# 2026-10-01 (IJF R3, PROGRESS.md §86): RSVST_AZ=<tag> (dualfixed | dualorig) switches the AZ-family forecasts (SV/RSV-AZ-SN/ST and the
# psi-free RSV-AZ variants) to the corrected-recomputation aggregates (Ox/az_recomputation/aggregated/<tag>/ (override with RSVST_AZ_DIR), read in the
# merge scripts) and appends _<tag> to the result/figure directories (after the RSVST_HL suffix). Old results stay untouched.
az.tag <- Sys.getenv("RSVST_AZ", "")
if (nzchar(az.tag)) {
  for (.p in intersect(c("path.result", "path.out", "path.figure"), ls())) {
    assign(.p, sub("/?$", paste0("_", az.tag, "/"), get(.p))); dir.create(get(.p), showWarnings = FALSE, recursive = TRUE)
  }
  cat("RSVST_AZ =", az.tag, "\n")
}

if (!dir.exists(path.out)) dir.create(path.out, recursive = TRUE)

# true returns and realized measures
data.true <- read_csv(paste0(path.data, true.sets, ".csv"), show_col_types = FALSE)
date <- data.true$yyyymmdd
rcc <- data.true$rcc * 100
rv5.hl <- data.true$rv5_hl * 10000
rk.th2.hl <- data.true$rk_th2_hl * 10000
bv.hl <- data.true$bv_hl * 10000
medrv.hl <- data.true$medrv_hl * 10000

# REGARCH psi=1 (baseline) and psi-free (varphiF) forecasts
f.regarch.n.psi1 <- read_excel(paste0(path.regarch, "REGARCHn_", data.name, "_2020.xlsx"))
f.regarch.t.psi1 <- read_excel(paste0(path.regarch, "REGARCHt_", data.name, "_2020.xlsx"))
f.regarch.n.psifree <- read_excel(paste0(path.regarch, "REGARCHnvarphiF_", data.name, "_2020.xlsx"))
f.regarch.t.psifree <- read_excel(paste0(path.regarch, "REGARCHtvarphiF_", data.name, "_2020.xlsx"))

# build model.names interleaving psi=1 / psi-free per distribution, e.g.
# "RSV-AZ-N (psi=1)", "RSV-AZ-N (psi-free)", "RSV-AZ-T (psi=1)", ...,
# followed by "REGARCH-N (psi=1)", "REGARCH-N (psi-free)", "REGARCH-T (psi=1)", "REGARCH-T (psi-free)"
model.names <- c(as.vector(rbind(paste0("RSV-", dists, " (psi=1)"),
                                  paste0("RSV-", dists, " (psi-free)"))),
                  "REGARCH-N (psi=1)", "REGARCH-N (psi-free)",
                  "REGARCH-T (psi=1)", "REGARCH-T (psi-free)")

# containers
vol.cols <- list()
var1p.cols <- list()
var5p.cols <- list()
es1p.cols <- list()
es5p.cols <- list()

for (i in seq_len(dist.n)) {
  key <- names(dists)[i]
  model.psi1 <- paste0("rsv_", key)
  model.psifree <- paste0("rsv_", key, "_psi")

  f.var.psi1 <- read_csv(paste0(p.in.az(model.psi1), "var_", model.psi1, ".csv"), show_col_types = FALSE)
  f.es.psi1 <- read_csv(paste0(p.in.az(model.psi1), "es_", model.psi1, ".csv"), show_col_types = FALSE)
  f.volf.psi1 <- read_csv(paste0(p.in.az(model.psi1), "volf_", model.psi1, ".csv"), show_col_types = FALSE)

  f.var.psifree <- read_csv(paste0(p.in.az(model.psifree), "var_", model.psifree, ".csv"), show_col_types = FALSE)
  f.es.psifree <- read_csv(paste0(p.in.az(model.psifree), "es_", model.psifree, ".csv"), show_col_types = FALSE)
  f.volf.psifree <- read_csv(paste0(p.in.az(model.psifree), "volf_", model.psifree, ".csv"), show_col_types = FALSE)

  is.fs <- grepl("^fs", key)
  # AZ/GH: var_/es_ columns are Var1..Var4 = 0.5%,1%,5%,10% quantiles; volf_ Var3 = mean(exp(h))
  # FS: var_/es_ columns are VaR_1p/VaR_5p directly; volf_ MeanExp = mean(exp(h))
  if (is.fs) {
    vol.cols[[paste0(model.psi1, "_vol")]] <- f.volf.psi1$MeanExp
    vol.cols[[paste0(model.psifree, "_vol")]] <- f.volf.psifree$MeanExp
    var1p.cols[[paste0(model.psi1, "_var1p")]] <- f.var.psi1$VaR_1p
    var1p.cols[[paste0(model.psifree, "_var1p")]] <- f.var.psifree$VaR_1p
    var5p.cols[[paste0(model.psi1, "_var5p")]] <- f.var.psi1$VaR_5p
    var5p.cols[[paste0(model.psifree, "_var5p")]] <- f.var.psifree$VaR_5p
    es1p.cols[[paste0(model.psi1, "_es1p")]] <- f.es.psi1$ES_1p
    es1p.cols[[paste0(model.psifree, "_es1p")]] <- f.es.psifree$ES_1p
    es5p.cols[[paste0(model.psi1, "_es5p")]] <- f.es.psi1$ES_5p
    es5p.cols[[paste0(model.psifree, "_es5p")]] <- f.es.psifree$ES_5p
  } else {
    vol.cols[[paste0(model.psi1, "_vol")]] <- f.volf.psi1$Var3
    vol.cols[[paste0(model.psifree, "_vol")]] <- f.volf.psifree$Var3
    var1p.cols[[paste0(model.psi1, "_var1p")]] <- f.var.psi1$Var2
    var1p.cols[[paste0(model.psifree, "_var1p")]] <- f.var.psifree$Var2
    var5p.cols[[paste0(model.psi1, "_var5p")]] <- f.var.psi1$Var3
    var5p.cols[[paste0(model.psifree, "_var5p")]] <- f.var.psifree$Var3
    es1p.cols[[paste0(model.psi1, "_es1p")]] <- f.es.psi1$Var2
    es1p.cols[[paste0(model.psifree, "_es1p")]] <- f.es.psifree$Var2
    es5p.cols[[paste0(model.psi1, "_es5p")]] <- f.es.psi1$Var3
    es5p.cols[[paste0(model.psifree, "_es5p")]] <- f.es.psifree$Var3
  }
}

# REGARCH-N/REGARCH-T, psi=1 and psi-free
vol.cols[["regarch_n_psi1_vol"]] <- f.regarch.n.psi1$`volatility forecast (mean)`
vol.cols[["regarch_n_psifree_vol"]] <- f.regarch.n.psifree$`volatility forecast (mean)`
vol.cols[["regarch_t_psi1_vol"]] <- f.regarch.t.psi1$`volatility forecast (mean)`
vol.cols[["regarch_t_psifree_vol"]] <- f.regarch.t.psifree$`volatility forecast (mean)`

var1p.cols[["regarch_n_psi1_var1p"]] <- f.regarch.n.psi1$`1% VaR`
var1p.cols[["regarch_n_psifree_var1p"]] <- f.regarch.n.psifree$`1% VaR`
var1p.cols[["regarch_t_psi1_var1p"]] <- f.regarch.t.psi1$`1% VaR`
var1p.cols[["regarch_t_psifree_var1p"]] <- f.regarch.t.psifree$`1% VaR`

var5p.cols[["regarch_n_psi1_var5p"]] <- f.regarch.n.psi1$`5% VaR`
var5p.cols[["regarch_n_psifree_var5p"]] <- f.regarch.n.psifree$`5% VaR`
var5p.cols[["regarch_t_psi1_var5p"]] <- f.regarch.t.psi1$`5% VaR`
var5p.cols[["regarch_t_psifree_var5p"]] <- f.regarch.t.psifree$`5% VaR`

es1p.cols[["regarch_n_psi1_es1p"]] <- f.regarch.n.psi1$`1% ES`
es1p.cols[["regarch_n_psifree_es1p"]] <- f.regarch.n.psifree$`1% ES`
es1p.cols[["regarch_t_psi1_es1p"]] <- f.regarch.t.psi1$`1% ES`
es1p.cols[["regarch_t_psifree_es1p"]] <- f.regarch.t.psifree$`1% ES`

es5p.cols[["regarch_n_psi1_es5p"]] <- f.regarch.n.psi1$`5% ES`
es5p.cols[["regarch_n_psifree_es5p"]] <- f.regarch.n.psifree$`5% ES`
es5p.cols[["regarch_t_psi1_es5p"]] <- f.regarch.t.psi1$`5% ES`
es5p.cols[["regarch_t_psifree_es5p"]] <- f.regarch.t.psifree$`5% ES`

df.vol <- data.frame(date, rv5.hl, rk.th2.hl, bv.hl, medrv.hl, vol.cols)
names(df.vol) <- c("date", "RV5", "RKth2", "BV", "MedRV", model.names)

df.var.1p <- data.frame(date, rcc, var1p.cols)
names(df.var.1p) <- c("date", "Return", model.names)
df.var.5p <- data.frame(date, rcc, var5p.cols)
names(df.var.5p) <- c("date", "Return", model.names)

df.es.1p <- data.frame(date, rcc, es1p.cols)
names(df.es.1p) <- c("date", "Return", model.names)
df.es.5p <- data.frame(date, rcc, es5p.cols)
names(df.es.5p) <- c("date", "Return", model.names)

# 2026-10-11 (IJF R3, PROGRESS.md §100-101): RSVST_MASK_DIR=<dir> drops the rows whose date is listed in <dir>/<NAME>_exclude_dates.csv
# (NAME = DJIA | N225 | DJIA_2020 | N225_2020) from every merged data frame, so that all models and proxies are evaluated on the same
# reduced sample (exclusion sensitivity for the windows that did not pass the MCMC confirmation pass; rule §7-6). Used with RSVST_AZ=<tag>_excl.
.mask.dir <- Sys.getenv("RSVST_MASK_DIR", "")
if (nzchar(.mask.dir)) {
  .mask.file <- file.path(.mask.dir, paste0(paste0(data.name, "_2020"), "_exclude_dates.csv")); .mask <- as.integer(read.csv(.mask.file)$date); .n0 <- nrow(df.vol)
  for (.d in c("df.vol", "df.vol.med", "df.var.1p", "df.var.5p", "df.es.1p", "df.es.5p")) if (exists(.d)) assign(.d, get(.d)[!(as.integer(get(.d)$date) %in% .mask), ])
  cat(sprintf("RSVST_MASK_DIR: %s, excluded %d of %d rows (%d dates listed)\n", .mask.file, .n0 - nrow(df.vol), .n0, length(.mask)))
}
write_csv(df.vol, paste0(path.out, data.name, "_forecast_vol_psi_sensitivity.csv"))
write_csv(df.var.1p, paste0(path.out, data.name, "_forecast_var1p_psi_sensitivity.csv"))
write_csv(df.var.5p, paste0(path.out, data.name, "_forecast_var5p_psi_sensitivity.csv"))
write_csv(df.es.1p, paste0(path.out, data.name, "_forecast_es1p_psi_sensitivity.csv"))
write_csv(df.es.5p, paste0(path.out, data.name, "_forecast_es5p_psi_sensitivity.csv"))

cat("Done. Wrote 5 files to", path.out, "\n")
