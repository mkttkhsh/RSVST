rm(list = ls())

#### Settings ####

# packages
library(readr)
library(ggplot2)
library(readxl)
library(reshape2)
library(MCS)
library(sandwich)

# function
source("function_plot_forecast_2026.R")

# set path
path.data <- "../Data/"
path.result <- "./prediction_results_n225_2020/"
path.figure <- "./figure_n225_2020/"

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
data.name <- "N225"

# define volatility proxies
rv.sets <- c("RV5", "RKth2", "BV", "MedRV")
rv.n <- length(rv.sets)

# set model names
model.names <- c(paste0("SV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("RSV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("EGARCH-", c("N", "T")),
                 paste0("REGARCH-", c("N", "T")))
model.n <- length(model.names)

###### select volatility proxies #####
# This script processes ONE proxy per invocation (unlike the DJIA COVID script,
# which loops over all four). The index was hardcoded to 4, so re-running it
# refreshed only the MedRV column and silently left the other three columns of
# the published table untouched. Keep the default at 4 but allow an override:
#   Rscript examine_forecast_vol_n225_2020.R 1   # RV5
rv.i <- 4
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) > 0) rv.i <- as.integer(.args[1])
stopifnot(rv.i %in% seq_along(rv.sets))
rv.name <- rv.sets[rv.i]

#### Set Prediction Results ####

# read volatility forecasts
f.vol <- read_csv(paste0(path.result, data.name, "_forecast_vol.csv"))
f.vol.med <- read_csv(paste0(path.result, data.name, "_forecast_vol_med.csv"))
# 2026-09-24 (R3, AE comment 2): the point forecast is the predictive mean of exp(h) for ALL models.
# Set RSVST_POINT=median for the median variant (robustness, supplement): models with a column in
# *_forecast_vol_med.csv are switched, the others stay at the mean; output gets the suffix "_median".
point.type <- Sys.getenv("RSVST_POINT", "mean")
out.suffix <- ""
if (point.type == "median") {
  med.models <- intersect(model.names, names(f.vol.med))
  med.models <- med.models[sapply(med.models, function(m) all(is.finite(f.vol.med[[m]]) & f.vol.med[[m]] > 0))]  # skip placeholder columns (e.g. RSV-GH-ST N225 COVID = 0)
  f.vol[med.models] <- f.vol.med[med.models]
  out.suffix <- "_median"
  cat("median variant: switched", paste(med.models, collapse = ", "), "\n")
}
f.vol <- f.vol[c("date", rv.name, model.names)] # select volatility proxy

# convert yyyymmdd to dates
yyyymmdd <- f.vol[["date"]]
f.vol[["date"]] <- as.Date(as.character(f.vol[["date"]]),format = '%Y%m%d')
# names(f.vol) <- c("date", "RK", model.names)

# revision (2026-09-24): SV-GH-ST is no longer excluded (forecast-step bug fixed, all 242 windows
# re-estimated; PROGRESS.md section 74). All 18 models are evaluated.

#### Evaluation ####

##### compute MSE and QLIKE #####
n <- nrow(f.vol)
k <- ncol(f.vol) - 2
loss.mse <- matrix(0, n, k)
loss.qlike <- matrix(0, n, k)
for(i in 1:k){
  loss.mse[1:n,i] <- (f.vol[[i+2]] - f.vol[[2]])^2 / 2
  loss.qlike[1:n,i] <- f.vol[[2]]/f.vol[[i+2]] - log(f.vol[[2]]/f.vol[[i+2]]) - 1
}
mse <- colMeans(loss.mse)
qlike <- colMeans(loss.qlike)
dt.mse.qlike <- data.frame(model.names, mse, qlike)

###### unconditional predictive ability test of Giacomini and White (2006) #####
qlike.d <- matrix(NA, k, k)
qlike.tstat <- matrix(NA, k, k)
qlike.tstat.prop <- matrix(NA, k, k)
qlike.tstat.pval <- matrix(NA, k, k)
m <- floor(0.75 * n^(1/3))
for(i in 1:(k-1)){
  if(sum(is.na(loss.qlike[, i])) == 0){
    for(j in (i+1):k){
      if(sum(is.na(loss.qlike[, j])) == 0){
        loss.qlike.d <- loss.qlike[, i] - loss.qlike[, j]
        qlike.d[i, j] <- mean(loss.qlike.d)
        fm <- lm(loss.qlike.d ~ 1)
        qlike.tstat.prop[i, j] <- mean(fm$fitted.values > 0)
        v <- as.numeric(NeweyWest(fm, lag = m - 1, prewhite = F, adjust = T))
        qlike.tstat[i, j] <- qlike.d[i, j] / sqrt(v)
        qlike.tstat.pval[i, j] <- 2 * (1 - pnorm(abs(qlike.tstat[i, j])))
      }
    }
  }
}

###### conditional predictive ability test of Giacomini and White (2006) #####
qlike.Qstat <- matrix(NA, k, k)
qlike.Qstat.prop <- matrix(NA, k, k)
qlike.Qstat.pval <- matrix(NA, k, k)
for(i in 1:(k-1)){
  if(sum(is.na(loss.qlike[, i])) == 0){
    for(j in (i+1):k){
      if(sum(is.na(loss.qlike[, j])) == 0){
        loss.qlike.d <- loss.qlike[, i] - loss.qlike[, j]
        fm <- lm(loss.qlike.d ~ c(0, loss.qlike.d[1:(n-1)]))
        qlike.Qstat.prop[i, j] <- mean(fm$fitted.values > 0)
        h <- cbind(rep(1,n), c(0, loss.qlike.d[1:(n-1)]))
        Z <- t(h) * loss.qlike.d
        Z.bar <- as.matrix(apply(Z, 1, mean))
        Omega <- Z %*% t(Z) / n
        qlike.Qstat[i,j] <- n * t(Z.bar) %*% solve(Omega) %*% Z.bar
        qlike.Qstat.pval[i,j] <- 1 - pchisq(qlike.Qstat[i,j], 2)
      }
    }
  }
}

###### Model confidence set of Hansen et al. (2011) for MSE #####
colnames(loss.mse) <- model.names
set.seed(20260801)  # reproducible MCS bootstrap (2026 revision-08-17)
mse.mcs <- MCSprocedure(loss.mse, alpha = 0.1, B=1000, statistic = "Tmax", k = 10)
mse.mcs.pval <- mse.mcs@show[, c(1,3)]
mse.mcs.pval.names <- row.names(mse.mcs.pval)
mse.mcs.n <- length(mse.mcs.pval.names)

###### Model confidence set of Hansen et al. (2011) for QLIKE #####
colnames(loss.qlike) <- model.names
set.seed(20260801)  # reproducible MCS bootstrap (2026 revision-08-17)
qlike.mcs <- MCSprocedure(loss.qlike, alpha = 0.1, B=1000, statistic = "Tmax", k = 10)
qlike.mcs.pval <- qlike.mcs@show[, c(1,3)]
# BUGFIX (2026-08-01): this used qlike.mcs@Info$model.names, which is the
# ORIGINAL input order, whereas qlike.mcs@show is ordered by MCS elimination.
# Matching one against the other by position assigned each model another
# model's MCS p-value whenever the two orders differed, which is almost
# always. The MSE block two lines above already used row.names correctly.
qlike.mcs.pval.names <- row.names(qlike.mcs.pval)
qlike.mcs.n <- length(qlike.mcs.pval.names)

###### Merge the results #####
stat <- -t(qlike.tstat)
stat[upper.tri(stat)] <- qlike.Qstat[upper.tri(qlike.Qstat)]
prop <- t(qlike.tstat.prop)
prop[upper.tri(prop)] <- qlike.Qstat.prop[upper.tri(qlike.Qstat.prop)]
pval <- t(qlike.tstat.pval)
pval[upper.tri(pval)] <- qlike.Qstat.pval[upper.tri(qlike.Qstat.pval)]
mcs.pval <- matrix(0, model.n, 2)
for (i in 1:model.n) {
  for (j in 1:mse.mcs.n) {
    if(model.names[i] == mse.mcs.pval.names[j]){
      mcs.pval[i, 1] <- mse.mcs.pval[j, 2]
    }
  }
  for (j in 1:qlike.mcs.n) {
    if(qlike.mcs.n > 1){
      if(model.names[i] == qlike.mcs.pval.names[j]){
        mcs.pval[i, 2] <- qlike.mcs.pval[j, 2]
      }
    }
    if(qlike.mcs.n == 1){
      if(model.names[i] == qlike.mcs.pval.names[j]){
        mcs.pval[i, 2] <- qlike.mcs.pval[2]
      }
    }
  }
}
df.qlike <- data.frame(mse, mcs.pval[,1], qlike, mcs.pval[,2], stat, prop, pval, row.names = model.names)
names(df.qlike) <- c("MSE", "MCS p-value", "QLIKE", "MCS p-value", paste0("Test Stat vs ", sub("RSV-", "", model.names)),
                     paste0("Prop vs ", sub("RSV-", "", model.names)), 
                     paste0("P-value vs ", sub("RSV-", "", model.names)))
write.csv(df.qlike, paste0(path.result, data.name, "_", rv.name, "_forecast_vol_mse_qlike_mcs_gw06test", out.suffix, ".csv"))
