rm(list = ls())

# packages
library(tidyverse)

# set path
path.data <- "../Data/"
path.result <- "./prediction_results_2026/"

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

# read QLIKE losses
qlike.djia <- read_csv(paste0(path.result, "DJIA_forecast_vol_qlike.csv"))
qlike.n225 <- read_csv(paste0(path.result, "N225_forecast_vol_qlike.csv"))

# set model names
model.names <- qlike.djia$Model
model.n <- length(model.names)

# MCS membership marker following Hansen, Lunde and Nason (2011):
# report both the 75% MCS (p-value > 0.25) and the 90% MCS (p-value > 0.10).
# Note: the MCS p-value is level-independent; alpha only sets the membership
# threshold, so both sets are read off the same p-value (verified empirically).
#   ** : model belongs to the 75% MCS (hence also to the 90% MCS)
#   *  : model belongs to the 90% MCS only
mcs.mark <- function(p) if (p > 0.25) "^{**}" else if (p > 0.1) "^{*}" else ""

# set volatility proxies
rv.sets <- names(qlike.djia[, c(3, 5, 7, 9)])
rv.n <- length(rv.sets)

# extract QLIKE losses as matrix
mat.qlike.djia <- unname(as.matrix(qlike.djia[, c(-1,-2)]))
mat.qlike.n225 <- unname(as.matrix(qlike.n225[, c(-1,-2)]))

# average rank
average.score <- rowMeans(cbind(mat.qlike.djia[, c(1,3,5,7)], mat.qlike.n225[, c(1,3,5,7)]))
average.rank <- rowMeans(cbind(mat.qlike.djia[, 9:12], mat.qlike.n225[, 9:12]))

# output for Latex: RV5, RKth2, BV, MedRV
cat("\n ### output for Latex: RV5, RKth2, BV, MedRV ### \n")
for (i in 1:model.n) {
  # Models
  xxx <- model.names[i]
  # RV5 for DJIA
  xxx <- paste0(xxx, " & $", format(round(mat.qlike.djia[i, 1], 3), nsmall = 3))
  xxx <- paste0(xxx, mcs.mark(mat.qlike.djia[i, 2]))
  xxx <- paste0(xxx, "$")
  # RKth2 for DJIA
  xxx <- paste0(xxx, " & $", format(round(mat.qlike.djia[i, 3], 3), nsmall = 3))
  xxx <- paste0(xxx, mcs.mark(mat.qlike.djia[i, 4]))
  xxx <- paste0(xxx, "$")
  # BV for DJIA
  xxx <- paste0(xxx, " & $", format(round(mat.qlike.djia[i, 5], 3), nsmall = 3))
  xxx <- paste0(xxx, mcs.mark(mat.qlike.djia[i, 6]))
  xxx <- paste0(xxx, "$")
  # MedRV for DJIA
  xxx <- paste0(xxx, " & $", format(round(mat.qlike.djia[i, 7], 3), nsmall = 3))
  xxx <- paste0(xxx, mcs.mark(mat.qlike.djia[i, 8]))
  xxx <- paste0(xxx, "$")
  #
  xxx <- paste0(xxx, " & ")
  # RV5 for N225
  xxx <- paste0(xxx, " & $", format(round(mat.qlike.n225[i, 1], 3), nsmall = 3))
  xxx <- paste0(xxx, mcs.mark(mat.qlike.n225[i, 2]))
  xxx <- paste0(xxx, "$")
  # RKth2 for N225
  xxx <- paste0(xxx, " & $", format(round(mat.qlike.n225[i, 3], 3), nsmall = 3))
  xxx <- paste0(xxx, mcs.mark(mat.qlike.n225[i, 4]))
  xxx <- paste0(xxx, "$")
  # BV for N225
  xxx <- paste0(xxx, " & $", format(round(mat.qlike.n225[i, 5], 3), nsmall = 3))
  xxx <- paste0(xxx, mcs.mark(mat.qlike.n225[i, 6]))
  xxx <- paste0(xxx, "$")
  # MedRV for N225
  xxx <- paste0(xxx, " & $", format(round(mat.qlike.n225[i, 7], 3), nsmall = 3))
  xxx <- paste0(xxx, mcs.mark(mat.qlike.n225[i, 8]))
  xxx <- paste0(xxx, "$")
  # Average score
  xxx <- paste0(xxx, " & & $", format(round(average.score[i], 3), nsmall = 3), "$")
  # Average rank
  xxx <- paste0(xxx, " & $", format(round(average.rank[i], 1), nsmall = 1), "$")
  #
  xxx <- paste0(xxx, " \\\\")
  #
  cat(xxx, "\n")
}

# output for Latex: RV5 only
cat("\n ### output for Latex: RV5 only ### \n")
for (i in 1:model.n) {
  # Models
  xxx <- model.names[i]
  # RV5 for DJIA
  xxx <- paste0(xxx, " & $", format(round(mat.qlike.djia[i, 1], 4), nsmall = 4))
  xxx <- paste0(xxx, "$ & $(", format(round(mat.qlike.djia[i, 2], 2), nsmall = 2))
  xxx <- paste0(xxx, ")$")
  #
  xxx <- paste0(xxx, " & ")
  # RV5 for N225
  xxx <- paste0(xxx, " & $", format(round(mat.qlike.n225[i, 1], 4), nsmall = 4))
  xxx <- paste0(xxx, "$ & $(", format(round(mat.qlike.n225[i, 2], 2), nsmall = 2))
  xxx <- paste0(xxx, ")$")
  #
  xxx <- paste0(xxx, " \\\\")
  #
  cat(xxx, "\n")
}
