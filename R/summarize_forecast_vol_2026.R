rm(list = ls())

# packages
library(tidyverse)

# set path
path.data <- "../Data/"
path.result <- "./prediction_results_2026/"
path.figure <- "./figure2026/"

# 2026-09-28 (October 2026 revision, PROGRESS.md §81): RSVST_HL=r3corrected switches to the corrected Hansen-Lunde proxies
# (Data/*_hl_roll*_r3corrected.csv) and to separate output directories (*_r3corrected/); old results stay untouched.
hl.tag <- Sys.getenv("RSVST_HL", "")
if (nzchar(hl.tag)) {
  if (exists("true.sets")) true.sets <- paste0(true.sets, "_", hl.tag)
  for (.p in intersect(c("path.result", "path.out", "path.figure"), ls())) {
    assign(.p, sub("/?$", paste0("_", hl.tag, "/"), get(.p))); dir.create(get(.p), showWarnings = FALSE, recursive = TRUE)
  }
  cat("RSVST_HL =", hl.tag, "\n")
}

# 2026-10-01 (October 2026 revision, PROGRESS.md §86): RSVST_AZ=<tag> (dualfixed | dualorig) switches the AZ-family forecasts (SV/RSV-AZ-SN/ST and the
# psi-free RSV-AZ variants) to the corrected-recomputation aggregates (Ox/az_recomputation/aggregated/<tag>/ (override with RSVST_AZ_DIR), read in the
# merge scripts) and appends _<tag> to the result/figure directories (after the RSVST_HL suffix). Old results stay untouched.
az.tag <- Sys.getenv("RSVST_AZ", "")
if (nzchar(az.tag)) {
  for (.p in intersect(c("path.result", "path.out", "path.figure"), ls())) {
    assign(.p, sub("/?$", paste0("_", az.tag, "/"), get(.p))); dir.create(get(.p), showWarnings = FALSE, recursive = TRUE)
  }
  cat("RSVST_AZ =", az.tag, "\n")
}
path.ox <- "../Ox/"

# define sets of data
data.sets <- c("DJIA", "N225")
data.n <- length(data.sets)

# define volatility proxies
rv.sets <- c("RV5", "RKth2", "BV", "MedRV")
rv.n <- length(rv.sets)

# set model names
model.names <- c(paste0("SV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("RSV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("EGARCH-", c("N", "T")),
                 paste0("REGARCH-", c("N", "T")))
model.n <- length(model.names)

# set loss names
loss.names <- c("MSE", "QLIKE")
loss.n <- length(loss.names)

###### select data #####
data.i <- 2
.args <- commandArgs(trailingOnly = TRUE)  # 2026-09-24: Rscript summarize_forecast_vol_2026.R 1|2 (DJIA|N225); default N225 as before
if (length(.args) > 0) data.i <- as.integer(.args[1])
data.name <- data.sets[data.i]

# extract prediction losses for all data sets
df.mse <- data.frame(Model = model.names)
df.qlike <- data.frame(Model = model.names)
for (rv.i in 1:rv.n) {
  rv.name <- rv.sets[rv.i]
  
  # read prediction losses
  loss <- read_csv(paste0(path.result, data.name, "_", rv.name, "_forecast_vol_mse_qlike_mcs_gw06test.csv"))
  mse <- loss[, 2:3]
  qlike <- loss[, 4:5]

  # label
  names(mse) <- c(rv.name, paste0(rv.name, "_pval"))
  names(qlike) <- c(rv.name, paste0(rv.name, "_pval"))
  
  # merge
  df.mse <- data.frame(df.mse, mse)
  df.qlike <- data.frame(df.qlike, qlike)
}

# loss rankings
mse.rank <- data.frame(Model = model.names, apply(df.mse[, rv.sets], 2, min_rank))
qlike.rank <- data.frame(Model = model.names, apply(df.qlike[, rv.sets], 2, min_rank))

# average ranking
mse.rank$Average <- rowMeans(mse.rank[, -1])
qlike.rank$Average <- rowMeans(qlike.rank[, -1])

# label
names(mse.rank)[-1] <- paste0("Rank_", names(mse.rank)[-1])
names(qlike.rank)[-1] <- paste0("Rank_", names(qlike.rank)[-1])

# save results
write.csv(cbind(df.mse, mse.rank[,-1]), paste0(path.result, data.name, "_forecast_vol_mse.csv"))
write.csv(cbind(df.qlike, qlike.rank[,-1]), paste0(path.result, data.name, "_forecast_vol_qlike.csv"))

