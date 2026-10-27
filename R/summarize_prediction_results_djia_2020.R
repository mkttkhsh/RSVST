rm(list = ls())

# packages
library(tidyverse)

# set path
path.data <- "../Data/"
path.result <- "./prediction_results_djia_2020/"

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

# define data
data.name <- "DJIA"

# set model names
# IJF R3 (2026-09-24): SV-GH-ST restored (forecast-step bug fixed, re-estimated; PROGRESS.md section 74).
model.names <- c(paste0("SV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("RSV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST")),
                 paste0("EGARCH-", c("N", "T")),
                 paste0("REGARCH-", c("N", "T")))
model.n <- length(model.names)

# set loss names
loss.names <- c("QL(RV5)", "QL(RK)", "QL(BV)", "QL(Med)", "FZ0(1%)", "FZ0(5%)")
loss.n <- length(loss.names)

# read prediction losses
loss.qlike.rv5 <- read_csv(paste0(path.result, data.name, "_RV5_forecast_vol_mse_qlike_mcs_gw06test.csv"))[, 4]
loss.qlike.rk <- read_csv(paste0(path.result, data.name, "_RKth2_forecast_vol_mse_qlike_mcs_gw06test.csv"))[, 4]
loss.qlike.bv <- read_csv(paste0(path.result, data.name, "_BV_forecast_vol_mse_qlike_mcs_gw06test.csv"))[, 4]
loss.qlike.med <- read_csv(paste0(path.result, data.name, "_MedRV_forecast_vol_mse_qlike_mcs_gw06test.csv"))[, 4]
loss.fz0.1p <- read_csv(paste0(path.result, data.name, "_forecast_var_es_1p_rate_fz0_mcs_dq_gw06test.csv"))[, 3]
loss.fz0.5p <- read_csv(paste0(path.result, data.name, "_forecast_var_es_5p_rate_fz0_mcs_dq_gw06test.csv"))[, 3]

# merge
loss <- data.frame(Model = model.names, loss.qlike.rv5, loss.qlike.rk, loss.qlike.bv, loss.qlike.med, loss.fz0.1p, loss.fz0.5p)

# loss rankings
loss.rank <- data.frame(Model = model.names, apply(loss[, -1], 2, min_rank))

# average ranking
loss.rank.avg <- data.frame(Model = model.names, average.rank = rowMeans(loss.rank[, -1]))

# output for latex
mat.loss <- as.matrix(loss[, -1])
mat.rank <- as.matrix(loss.rank[, -1])
mat.avg <- as.matrix(loss.rank.avg[, -1])
for (i in 1:model.n) {
  xxx <- model.names[i]
  for (j in 1:loss.n) {
    xxx <- paste0(xxx, " & $", format(round(mat.loss[i, j], 2), nsmall = 2), "$")
    xxx <- paste0(xxx, " [", mat.rank[i, j], "]")
  }
  xxx <- paste0(xxx, " & $", format(round(mat.avg[i, 1], 2), nsmall = 2), "$")
  xxx <- paste0(xxx, " {\\smallskip} \\\\")
  cat(xxx, "\n")
}
# ranking only
cat("\n### Ranking only ###")
for (i in 1:model.n) {
  xxx <- model.names[i]
  for (j in 1:loss.n) {
    xxx <- paste0(xxx, " & \\multicolumn{1}{c}{")
    if(mat.rank[i,j] < 10)
      xxx <- paste0(xxx, "\\hspace{0.4em}")
    xxx <- paste0(xxx, mat.rank[i, j], "}")
    if(j == 4){
      xxx <- paste0(xxx, " & ")
    }
  }
  xxx <- paste0(xxx, " & $", format(round(mat.avg[i, 1], 2), nsmall = 2), "$")
  xxx <- paste0(xxx, " \\\\")
  cat(xxx, "\n")
}

# save rankings
df.rank.all <- merge(loss.rank, loss.rank.avg)
names(df.rank.all) <- c("Model",
                        "NK-QLIKE-RV5", "NK-QLIKE-RK", "NK-QLIKE-BV", "NK-QLIKE-Med",
                        "NK-FZ0(1%)", "NK-FZ0(5%)",
                        "Average")
# Model列をmodel.namesの順序を持つ因子型に変換し、その順序で並べ替える
df.rank.all <- df.rank.all %>%
  mutate(Model = factor(Model, levels = model.names)) %>%
  arrange(Model)
write.csv(df.rank.all, paste0(path.result, "rank_loss.csv"))
