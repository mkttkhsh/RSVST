# IJF R2 (2026): volatility forecast evaluation (QLIKE/MSE, GW06 test, MCS)
# for the psi sensitivity analysis -- RSV psi=1 vs psi-free, DJIA COVID
# period. DJIA counterpart of examine_forecast_vol_n225_2020_psi_sensitivity.R.
# Unlike N225 (7 RSV distributions), the DJIA COVID batch has no FS models,
# so this covers the 4 AZ-family RSV distributions x {psi=1, psi-free} plus
# REGARCH-N/T x {psi=1, psi-free} = 12 models. RSV-GH-ST is excluded (both
# psi variants): its volatility forecast is off-scale in essentially every
# window (see merge_prediction_results_djia_2020_psi_sensitivity.R).
#
# 2026-08-17: REGARCH added. Its forecasts were delivered with the N225 ones
# on 2026-07-29 but the DJIA merge script never read them.
#
# Run separately for each rv.i (1..4 = RV5/RKth2/BV/MedRV).

rm(list = ls())

library(readr)
library(MCS)
library(sandwich)

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
data.name <- "DJIA"

rv.sets <- c("RV5", "RKth2", "BV", "MedRV")
rv.n <- length(rv.sets)

dists <- c("AZ-N", "AZ-T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST")
model.names <- c(as.vector(rbind(paste0("RSV-", dists, " (psi=1)"),
                                  paste0("RSV-", dists, " (psi-free)"))),
                  "REGARCH-N (psi=1)", "REGARCH-N (psi-free)",
                  "REGARCH-T (psi=1)", "REGARCH-T (psi-free)")
model.n <- length(model.names)
# Guard against the model-count drift that produced the missing DJIA GW
# figures (see PROGRESS.md 44). 6 RSV distributions x 2 psi variants
# plus 4 REGARCH columns.
stopifnot(model.n == 16)

###### run one volatility proxy #####
run_examine_vol_psi <- function(rv.i) {
rv.name <- rv.sets[rv.i]
cat(sprintf("\n=== examine_forecast_vol_djia_2020_psi_sensitivity: %s ===\n", rv.name))

#### Set Prediction Results ####
f.vol <- read_csv(paste0(path.result, data.name, "_forecast_vol_psi_sensitivity.csv"), show_col_types = FALSE)
f.vol[["date"]] <- as.Date(as.character(f.vol[["date"]]), format = '%Y%m%d')
f.vol <- f.vol[, c("date", rv.name, model.names)]

#### Evaluation ####

##### compute MSE and QLIKE #####
n <- nrow(f.vol)
k <- ncol(f.vol) - 2
loss.mse <- matrix(0, n, k)
loss.qlike <- matrix(0, n, k)
for (i in 1:k) {
  loss.mse[1:n, i] <- (f.vol[[i + 2]] - f.vol[[2]])^2 / 2
  loss.qlike[1:n, i] <- f.vol[[2]] / f.vol[[i + 2]] - log(f.vol[[2]] / f.vol[[i + 2]]) - 1
}
mse <- colMeans(loss.mse)
qlike <- colMeans(loss.qlike)

###### unconditional predictive ability test of Giacomini and White (2006) #####
qlike.d <- matrix(NA, k, k)
qlike.tstat <- matrix(NA, k, k)
qlike.tstat.prop <- matrix(NA, k, k)
qlike.tstat.pval <- matrix(NA, k, k)
m <- floor(0.75 * n^(1 / 3))
for (i in 1:(k - 1)) {
  for (j in (i + 1):k) {
    loss.qlike.d <- loss.qlike[, i] - loss.qlike[, j]
    qlike.d[i, j] <- mean(loss.qlike.d)
    fm <- lm(loss.qlike.d ~ 1)
    qlike.tstat.prop[i, j] <- mean(fm$fitted.values > 0)
    v <- as.numeric(NeweyWest(fm, lag = m - 1, prewhite = F, adjust = T))
    qlike.tstat[i, j] <- qlike.d[i, j] / sqrt(v)
    qlike.tstat.pval[i, j] <- 2 * (1 - pnorm(abs(qlike.tstat[i, j])))
  }
}

###### conditional predictive ability test of Giacomini and White (2006) #####
qlike.Qstat <- matrix(NA, k, k)
qlike.Qstat.prop <- matrix(NA, k, k)
qlike.Qstat.pval <- matrix(NA, k, k)
for (i in 1:(k - 1)) {
  for (j in (i + 1):k) {
    loss.qlike.d <- loss.qlike[, i] - loss.qlike[, j]
    fm <- lm(loss.qlike.d ~ c(0, loss.qlike.d[1:(n - 1)]))
    qlike.Qstat.prop[i, j] <- mean(fm$fitted.values > 0)
    h <- cbind(rep(1, n), c(0, loss.qlike.d[1:(n - 1)]))
    Z <- t(h) * loss.qlike.d
    Z.bar <- as.matrix(apply(Z, 1, mean))
    Omega <- Z %*% t(Z) / n
    qlike.Qstat[i, j] <- n * t(Z.bar) %*% solve(Omega) %*% Z.bar
    qlike.Qstat.pval[i, j] <- 1 - pchisq(qlike.Qstat[i, j], 2)
  }
}

###### Model confidence set of Hansen et al. (2011) #####
colnames(loss.mse) <- model.names
set.seed(20260801)  # reproducible MCS bootstrap (IJF R2, 2026-08-01)
mse.mcs <- MCSprocedure(loss.mse, alpha = 0.1, B = 1000, statistic = "Tmax", k = 10)
mse.mcs.pval <- mse.mcs@show[, c(1, 3)]
mse.mcs.pval.names <- row.names(mse.mcs.pval)
mse.mcs.n <- length(mse.mcs.pval.names)

colnames(loss.qlike) <- model.names
set.seed(20260801)  # reproducible MCS bootstrap (IJF R2, 2026-08-01)
qlike.mcs <- MCSprocedure(loss.qlike, alpha = 0.1, B = 1000, statistic = "Tmax", k = 10)
qlike.mcs.pval <- qlike.mcs@show[, c(1, 3)]
# BUGFIX (2026-08-01): was qlike.mcs@Info$model.names (original input order),
# matched by position against qlike.mcs@show (MCS elimination order). Each
# model therefore received another model's MCS p-value. The MSE block just
# above already used row.names correctly.
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
    if (model.names[i] == mse.mcs.pval.names[j]) mcs.pval[i, 1] <- mse.mcs.pval[j, 2]
  }
  for (j in 1:qlike.mcs.n) {
    if (qlike.mcs.n > 1) {
      if (model.names[i] == qlike.mcs.pval.names[j]) mcs.pval[i, 2] <- qlike.mcs.pval[j, 2]
    } else {
      if (model.names[i] == qlike.mcs.pval.names[j]) mcs.pval[i, 2] <- qlike.mcs.pval[2]
    }
  }
}
df.qlike <- data.frame(mse, mcs.pval[, 1], qlike, mcs.pval[, 2], stat, prop, pval, row.names = model.names)
names(df.qlike) <- c("MSE", "MCS p-value", "QLIKE", "MCS p-value",
                     paste0("Test Stat vs ", model.names),
                     paste0("Prop vs ", model.names),
                     paste0("P-value vs ", model.names))
write.csv(df.qlike, paste0(path.result, data.name, "_", rv.name, "_forecast_vol_mse_qlike_mcs_gw06test_psi_sensitivity.csv"))

cat("Done. rv.name =", rv.name, "\n")
}

#### Run for all volatility proxies ####
for (rv.i in 1:rv.n) {
  run_examine_vol_psi(rv.i)
}
