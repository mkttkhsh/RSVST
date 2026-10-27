# IJF R2 (2026): merge psi=1 (main spec) vs psi-free (sensitivity) RSV and
# REGARCH forecasts for the DJIA COVID period, side by side, for the psi
# sensitivity analysis (Reviewer 1 Comment 2). DJIA counterpart of
# merge_prediction_results_n225_2020_psi_sensitivity.R.
#
# Covers AZ x4 and FS x2, each with a psi=1 and a psi-free variant, plus
# REGARCH-N/T. RSV-GH-ST is excluded (see below).
#
# 2026-08-19: FS added. The psi-free FS runs for this sample were executed on
# the Sakura server and the per-window draws sat in
# Ox/DJIA_prediction_2020/RSV_FS_{SN,ST}_psi/ without ever being aggregated;
# merge_prediction_results_fs_psi_djia_2020.R now does that.
#
# Inputs: var_/es_/volf_<model>.csv per-model files produced by the DJIA
# COVID rolling estimation (server_upload_djia_2020[_sakura_reseed]),
# downloaded into reproducibility/Ox/DJIA_prediction_2020/.
#
# 2026-08-17: REGARCH added. The psi=1 (REGARCHn_/REGARCHt_DJIA_2020.xlsx)
# and psi-free (REGARCHnvarphiF_/REGARCHtvarphiF_DJIA_2020.xlsx) forecasts
# were delivered alongside the N225 ones (2026-07-29) and have the same
# column layout and the same 253 dates as the DJIA true-data file, but this
# script never read them, so the supplement stated that DJIA REGARCH results
# under a free psi were unavailable. They were available all along.

rm(list = ls())
library(readr)
library(readxl)

# IJF R2 (2026-07-31): RSV-GH-ST excluded (both psi=1 and psi-free) --
# its volatility forecast is off-scale relative to realized volatility in
# essentially every window, a pervasive structural divergence (see the main
# COVID section), not something a psi-sensitivity comparison can speak to.
# 4 RSV distributions, each with a psi=1 and a psi-free variant.
dists <- c("az_n" = "AZ-N", "az_t" = "AZ-T", "az_sn" = "AZ-SN", "az_st" = "AZ-ST",
           "fs_sn" = "FS-SN", "fs_st" = "FS-ST")
dist.n <- length(dists)

data.name <- "DJIA"
path.data <- "../Data/"
path.in <- paste0("../Ox/", data.name, "_prediction_2020/")
# 2026-10-01 (§86): RSVST_AZ -> RSV-AZ-SN/ST (psi=1 and psi-free) are read from the corrected-recomputation aggregates; AZ-N/T, FS, GH stay as before
path.in.az <- if (nzchar(Sys.getenv("RSVST_AZ", ""))) paste0(Sys.getenv("RSVST_AZ_DIR", "../Ox/az_recomputation/aggregated/"), sub("_excl$", "", Sys.getenv("RSVST_AZ")), "/", data.name, "_prediction_2020/") else path.in
p.in.az <- function(model) if (model %in% c("rsv_az_sn", "rsv_az_st", "rsv_az_sn_psi", "rsv_az_st_psi")) path.in.az else path.in
path.out <- "./prediction_results_djia_2020/"
path.regarch <- "../REGARCH_Bayes/"
true.sets <- "DJIA_rcc_rv5_rkth2_bv_medrv_2012_2020_hl_roll_tick"

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

# The REGARCH files carry their own date column. Line them up with the
# evaluation dates rather than trusting row order: a silent off-by-one would
# shift every loss and DQ statistic without any visible error.
for (nm in c("f.regarch.n.psi1", "f.regarch.t.psi1",
             "f.regarch.n.psifree", "f.regarch.t.psifree")) {
  stopifnot(identical(as.numeric(get(nm)[[1]]), as.numeric(date)))
}

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

  # FS aggregates carry only the 1% and 5% quantiles, whereas AZ carries four
  # with 1% in Var2 and 5% in Var3 (PROGRESS.md 36). Rename so the rest of the
  # loop can address both families the same way; volf already matches.
  rd <- function(kind, model) {
    x <- read_csv(paste0(p.in.az(model), kind, "_", model, ".csv"), show_col_types = FALSE)
    if (kind != "volf" && ncol(x) == 2) names(x) <- c("Var2", "Var3")
    else names(x) <- paste0("Var", seq_len(ncol(x)))
    x
  }
  f.var.psi1 <- rd("var", model.psi1)
  f.es.psi1 <- rd("es", model.psi1)
  f.volf.psi1 <- rd("volf", model.psi1)

  f.var.psifree <- rd("var", model.psifree)
  f.es.psifree <- rd("es", model.psifree)
  f.volf.psifree <- rd("volf", model.psifree)

  # AZ/GH: var_/es_ columns are Var1..Var4 = 0.5%,1%,5%,10% quantiles; volf_ Var3 = mean(exp(h))
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

# Guard against the model-count drift that produced the missing DJIA GW
# figures (see PROGRESS.md 44): 4 RSV distributions x 2 psi variants + 4
# REGARCH columns = 12.
stopifnot(dist.n == 6, length(model.names) == 2 * dist.n + 4,
          length(vol.cols) == length(model.names),
          length(var1p.cols) == length(model.names),
          length(es5p.cols) == length(model.names))

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
