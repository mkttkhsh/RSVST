rm(list = ls())

# packages
library(tidyverse)
library(ggrepel)
library(readxl)
library(reshape2)
library(MCS)
library(sandwich)

# load function
source("function_forecast_evaluation.R")

# set path
path.data <- "../Data/"
path.result <- "./prediction_results_djia_2020/"
path.figure <- "./figure_djia_2020/"

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

# define data
data.name <- "DJIA"

# set model names
# revision (2026-08-01): SV-AZ-ST restored after re-estimation (see
# merge_prediction_results_djia_2020.R for the full note). GH-ST is
# permanently excluded for both RSV and SV on this sample; FS models are
# FS models added 2026-08-06 after all four finished and passed qc_relative.py.
# revision (2026-09-24): SV-GH-ST restored (forecast-step bug fixed, 253 windows re-estimated; PROGRESS.md section 74).
model.names <- c(paste0("SV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("RSV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST")),
                 paste0("EGARCH-", c("N", "T")),
                 paste0("REGARCH-", c("N", "T")))
model.n <- length(model.names)

# set sets of target rates
alpha.sets <- c("1p", "5p")

# read VaR forecasts
f.var1p <- read_csv(paste0(path.result, data.name, "_forecast_var1p.csv"))
f.var5p <- read_csv(paste0(path.result, data.name, "_forecast_var5p.csv"))

# read ES forecasts
f.es1p <- read_csv(paste0(path.result, data.name, "_forecast_es1p.csv"))
f.es5p <- read_csv(paste0(path.result, data.name, "_forecast_es5p.csv"))

# set dates
f.date <- as.Date(as.character(f.var1p[["date"]]), format = '%Y%m%d')

# number of observations
n <- length(f.date)

#### Computations ####

##### compute violations and various lossfunctions #####
# loss function of Clements et al. (2008), 
# standardized residuals for VaR and ES, and FZ0 loss function of Patton et al. (2019) 
loss1p <- array(NA, c(model.n, n, 5))
loss5p <- array(NA, c(model.n, n, 5))
for (i in 1:model.n) {
  loss1p[i, , ] <- loss_var_es_fz0(f.var1p[[2]], f.var1p[[i + 2]], f.es1p[[i + 2]], 0.01)
  loss5p[i, , ] <- loss_var_es_fz0(f.var5p[[2]], f.var5p[[i + 2]], f.es5p[[i + 2]], 0.05)
}

##### Dynamic quantile test of Engle and Manganelli (2004) #####
DQtest1p <- matrix(NA, model.n, 2)
DQtest5p <- matrix(NA, model.n, 2)
for (i in 1:model.n) {
  DQtest1p[i, 1] <- DQtest(loss1p[i, , 3], f.var1p[[i + 2]])
  DQtest1p[i, 2] <- DQtest(loss1p[i, , 4], f.es1p[[i + 2]])
  DQtest5p[i, 1] <- DQtest(loss5p[i, , 3], f.var5p[[i + 2]])
  DQtest5p[i, 2] <- DQtest(loss5p[i, , 4], f.es5p[[i + 2]])
}

# compute violation rates and average losses
loss1p.mean <- apply(loss1p, c(1,3), mean)
loss5p.mean <- apply(loss5p, c(1,3), mean)

# extract FZ0 loss
loss.fz0.1p <- t(loss1p[, , 5])
loss.fz0.5p <- t(loss5p[, , 5])

##### Model confidence set of Hansen et al. (2011) using FZ0 loss #####
# 1%
colnames(loss.fz0.1p) <- model.names
set.seed(20260801)  # reproducible MCS bootstrap (2026 revision)
loss.fz0.1p.mcs <- MCSprocedure(loss.fz0.1p, alpha = 0.1, B=1000, statistic = "Tmax", k = 10)
loss.fz0.1p.mcs.pval <- loss.fz0.1p.mcs@show[, c(1,3)]
loss.fz0.1p.mcs.pval.names <- row.names(loss.fz0.1p.mcs.pval)
loss.fz0.1p.mcs.n <- length(loss.fz0.1p.mcs.pval.names)
# 5%
colnames(loss.fz0.5p) <- model.names
set.seed(20260801)  # reproducible MCS bootstrap (2026 revision)
loss.fz0.5p.mcs <- MCSprocedure(loss.fz0.5p, alpha = 0.1, B=1000, statistic = "Tmax", k = 10)
loss.fz0.5p.mcs.pval <- loss.fz0.5p.mcs@show[, c(1,3)]
loss.fz0.5p.mcs.pval.names <- row.names(loss.fz0.5p.mcs.pval)
loss.fz0.5p.mcs.n <- length(loss.fz0.5p.mcs.pval.names)

##### unconditional predictive ability test of Giacomini and White (2006) #####
loss1p.d <- matrix(NA, model.n, model.n)
loss5p.d <- matrix(NA, model.n, model.n)
loss1p.tstat <- matrix(NA, model.n, model.n)
loss5p.tstat <- matrix(NA, model.n, model.n)
loss1p.tstat.pval <- matrix(NA, model.n, model.n)
loss5p.tstat.pval <- matrix(NA, model.n, model.n)
loss1p.tstat.prop <- matrix(NA, model.n, model.n)
loss5p.tstat.prop <- matrix(NA, model.n, model.n)
m <- floor(0.75 * n^(1/3))
for(i in 1:(model.n-1)){
  for(j in (i+1):model.n){
    # 1%
    loss1p.d[i, j] <- loss1p.mean[i, 5] - loss1p.mean[j, 5]
    diff1p <- loss.fz0.1p[, i] - loss.fz0.1p[, j]
    fm <- lm(diff1p ~ 1)
    loss1p.tstat.prop[i, j] <- mean(fm$fitted.values > 0)
    v <- as.numeric(NeweyWest(fm, lag = m - 1, prewhite = F, adjust = T))
    loss1p.tstat[i, j] <- loss1p.d[i, j] / sqrt(v)
    loss1p.tstat.pval[i, j] <- 2 * (1 - pnorm(abs(loss1p.tstat[i, j])))
    # 5%
    loss5p.d[i, j] <- loss5p.mean[i, 5] - loss5p.mean[j, 5]
    diff5p <- loss.fz0.5p[, i] - loss.fz0.5p[, j]
    fm <- lm(diff5p ~ 1)
    loss5p.tstat.prop[i, j] <- mean(fm$fitted.values > 0)
    v <- as.numeric(NeweyWest(fm, lag = m - 1, prewhite = F, adjust = T))
    loss5p.tstat[i, j] <- loss5p.d[i, j] / sqrt(v)
    loss5p.tstat.pval[i, j] <- 2 * (1 - pnorm(abs(loss5p.tstat[i, j])))
  }
}

##### predictive ability test based on Giacomini and White (2006) using FZ0 loss #####
loss1p.Qstat <- matrix(NA, model.n, model.n)
loss5p.Qstat <- matrix(NA, model.n, model.n)
loss1p.Qstat.pval <- matrix(NA, model.n, model.n)
loss5p.Qstat.pval <- matrix(NA, model.n, model.n)
loss1p.Qstat.prop <- matrix(NA, model.n, model.n)
loss5p.Qstat.prop <- matrix(NA, model.n, model.n)
for(i in 1:(model.n-1)){
  for(j in (i+1):model.n){
    # 1%
    diff1p <- loss.fz0.1p[, i] - loss.fz0.1p[, j]
    fm <- lm(diff1p ~ c(0, diff1p[1:(n-1)]))
    loss1p.Qstat.prop[i, j] <- mean(fm$fitted.values > 0)
    h <- cbind(rep(1,n), c(0, diff1p[1:(n-1)]))
    Z <- t(h) * diff1p
    Z.bar <- as.matrix(apply(Z, 1, mean))
    Omega <- Z %*% t(Z) / n
    loss1p.Qstat[i, j] <- n * t(Z.bar) %*% solve(Omega) %*% Z.bar
    loss1p.Qstat.pval[i, j] <- 1 - pchisq(loss1p.Qstat[i,j], 2)
    # 5%
    diff5p <- loss.fz0.5p[, i] - loss.fz0.5p[, j]
    fm <- lm(diff5p ~ c(0, diff5p[1:(n-1)]))
    loss5p.Qstat.prop[i, j] <- mean(fm$fitted.values > 0)
    h <- cbind(rep(1,n), c(0, diff5p[1:(n-1)]))
    Z <- t(h) * diff5p
    Z.bar <- as.matrix(apply(Z, 1, mean))
    Omega <- Z %*% t(Z) / n
    loss5p.Qstat[i, j] <- n * t(Z.bar) %*% solve(Omega) %*% Z.bar
    loss5p.Qstat.pval[i, j] <- 1 - pchisq(loss5p.Qstat[i,j], 2)
  }
}

#### Output ####

##### Merge the results #####
# 1%
stat1p <- -t(loss1p.tstat)
stat1p[upper.tri(stat1p)] <- loss1p.Qstat[upper.tri(loss1p.Qstat)]
pval1p <- t(loss1p.tstat.pval)
pval1p[upper.tri(pval1p)] <- loss1p.Qstat.pval[upper.tri(loss1p.Qstat.pval)]
prop1p <- t(loss1p.tstat.prop)
prop1p[upper.tri(prop1p)] <- loss1p.Qstat.prop[upper.tri(loss1p.Qstat.prop)]
mcs.pval1p <- numeric(model.n)
for (i in 1:model.n) {
  for (j in 1:loss.fz0.1p.mcs.n) {
    if(model.names[i] == loss.fz0.1p.mcs.pval.names[j]){
      mcs.pval1p[i] <- loss.fz0.1p.mcs.pval[j, 2]
    }
  }
}
# 5%
stat5p <- -t(loss5p.tstat)
stat5p[upper.tri(stat5p)] <- loss5p.Qstat[upper.tri(loss5p.Qstat)]
pval5p <- t(loss5p.tstat.pval)
pval5p[upper.tri(pval5p)] <- loss5p.Qstat.pval[upper.tri(loss5p.Qstat.pval)]
prop5p <- t(loss5p.tstat.prop)
prop5p[upper.tri(prop5p)] <- loss5p.Qstat.prop[upper.tri(loss5p.Qstat.prop)]
mcs.pval5p <- numeric(model.n)
for (i in 1:model.n) {
  for (j in 1:loss.fz0.5p.mcs.n) {
    if(model.names[i] == loss.fz0.5p.mcs.pval.names[j]){
      mcs.pval5p[i] <- loss.fz0.5p.mcs.pval[j, 2]
    }
  }
}

##### save violation rates, average losses, and DQ test results #####
# 1%
df.loss1p <- data.frame(loss1p.mean[, c(1,5)], mcs.pval1p, DQtest1p, stat1p, prop1p, pval1p)
names(df.loss1p) <- c("1%", "FZ0", "MCS p-value", "DQ", "DES", 
                    paste0("Test Stat vs ", sub("RSV-", "", model.names)),
                    paste0("Prop vs ", sub("RSV-", "", model.names)),
                    paste0("P-value vs ", sub("RSV-", "", model.names)))
write.csv(df.loss1p, paste0(path.result, data.name, "_forecast_var_es_1p_rate_fz0_mcs_dq_gw06test.csv"))
# 5%
df.loss5p <- data.frame(loss5p.mean[, c(1,5)], mcs.pval5p, DQtest5p, stat5p, prop5p, pval5p)
names(df.loss5p) <- c("1%", "FZ0", "MCS p-value", "DQ", "DES", 
                      paste0("Test Stat vs ", sub("RSV-", "", model.names)),
                      paste0("Prop vs ", sub("RSV-", "", model.names)),
                      paste0("P-value vs ", sub("RSV-", "", model.names)))
write.csv(df.loss5p, paste0(path.result, data.name, "_forecast_var_es_5p_rate_fz0_mcs_dq_gw06test.csv"))

##### output for latex ######
# violation rate, FZ0 loss
mat.loss1p <- as.matrix(df.loss1p)
mat.loss5p <- as.matrix(df.loss5p)
cat("\n### Violation rates and FZ0 loss ###\n")
for (i in 1:model.n) {
  xxx <- model.names[i]
  # 1%
  xxx <- paste0(xxx, " & $", format(round(mat.loss1p[i, 1]*100, 2), nsmall = 2)) # violation rate
  xxx <- paste0(xxx, "$ & $", format(round(mat.loss1p[i, 4], 2), nsmall = 2)) # p-value of DQ test
  xxx <- paste0(xxx, "$ & $", format(round(mat.loss1p[i, 2], 4), nsmall = 4)) # FZ0 loss
  xxx <- paste0(xxx, "$ & $", format(round(mat.loss1p[i, 3], 2), nsmall = 2)) # p-value of MCS
  # MCS membership following Hansen, Lunde and Nason (2011): report both levels.
  #   ** : 75% MCS (p-value > 0.25)   * : 90% MCS only (0.10 < p-value <= 0.25)
  if(mat.loss1p[i, 3] > 0.25){
    xxx <- paste0(xxx, "^{**}")
  }
  else if(mat.loss1p[i, 3] > 0.1){
    xxx <- paste0(xxx, "^{*\\hspace{0.4em}}")
  }
  else{
    xxx <- paste0(xxx, "^{\\hspace{0.8em}}")
  }
  xxx <- paste0(xxx, "$")
  xxx <- paste0(xxx, " & ")
  # 5%
  xxx <- paste0(xxx, " & $", format(round(mat.loss5p[i, 1]*100, 2), nsmall = 2)) # violation rate
  xxx <- paste0(xxx, "$ & $", format(round(mat.loss5p[i, 4], 2), nsmall = 2)) # p-value of DQ test
  xxx <- paste0(xxx, "$ & $", format(round(mat.loss5p[i, 2], 4), nsmall = 4)) # FZ0 loss
  xxx <- paste0(xxx, "$ & $", format(round(mat.loss5p[i, 3], 2), nsmall = 2)) # p-value of MCS
  # MCS membership following Hansen, Lunde and Nason (2011): report both levels.
  #   ** : 75% MCS (p-value > 0.25)   * : 90% MCS only (0.10 < p-value <= 0.25)
  if(mat.loss5p[i, 3] > 0.25){
    xxx <- paste0(xxx, "^{**}")
  }
  else if(mat.loss5p[i, 3] > 0.1){
    xxx <- paste0(xxx, "^{*\\hspace{0.4em}}")
  }
  else{
    xxx <- paste0(xxx, "^{\\hspace{0.8em}}")
  }
  xxx <- paste0(xxx, "$")
  xxx <- paste0(xxx, " \\\\")
  cat(xxx, "\n")
}
