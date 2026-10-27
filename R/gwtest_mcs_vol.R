rm(list = ls())

# packages
library(readr)
library(ggplot2)
library(readxl)
library(reshape2)
library(MCS)
library(sandwich)

# read volatility forecasts
f.vol <- read_csv("forecast_vol.csv")

# set model names
model.names <- c(paste0("SV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("RSV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("EGARCH-", c("N", "T")),
                 paste0("REGARCH-", c("N", "T")))
model.n <- length(model.names)

# convert yyyymmdd to dates
yyyymmdd <- f.vol[["date"]]
f.vol[["date"]] <- as.Date(as.character(f.vol[["date"]]),format = '%Y%m%d')
names(f.vol) <- c("date", "RK", model.names)

# compute MSE and QLIKE
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

# unconditional predictive ability test of Giacomini and White (2006)
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

# conditional predictive ability test of Giacomini and White (2006)
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

# Model confidence set of Hansen et al. (2011) for MSE
colnames(loss.mse) <- model.names
mse.mcs <- MCSprocedure(loss.mse, alpha = 0.1, B=1000, statistic = "Tmax", k = 10)
mse.mcs.pval <- mse.mcs@show[, c(1,3)]
mse.mcs.pval.names <- row.names(mse.mcs.pval)
mse.mcs.n <- length(mse.mcs.pval.names)

# Model confidence set of Hansen et al. (2011) for QLIKE
colnames(loss.qlike) <- model.names
qlike.mcs <- MCSprocedure(loss.qlike, alpha = 0.1, B=1000, statistic = "Tmax", k = 10)
qlike.mcs.pval <- qlike.mcs@show[, c(1,3)]
qlike.mcs.pval.names <- row.names(qlike.mcs.pval)
qlike.mcs.n <- length(qlike.mcs.pval.names)

# Merge the results
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
    if(model.names[i] == qlike.mcs.pval.names[j]){
      mcs.pval[i, 2] <- qlike.mcs.pval[j, 2]
    }
  }
}
df.qlike <- data.frame(mse, mcs.pval[,1], qlike, mcs.pval[,2], stat, prop, pval, row.names = model.names)
names(df.qlike) <- c("MSE", "MCS p-value", "QLIKE", "MCS p-value", paste0("Test Stat vs ", sub("RSV-", "", model.names)),paste0("Prop vs ", sub("RSV-", "", model.names)), 
                     paste0("P-value vs ", sub("RSV-", "", model.names)))
# write.csv(df.qlike, "forecast_vol_mse_qlike_mcs_gw06test.csv")