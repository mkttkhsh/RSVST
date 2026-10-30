# function to compute VaR violations, loss function of Clements et al. (2008), and 
# standardized residuals and FZ0 loss function of Patton et al. (2019)
loss_var_es_fz0 <- function(y, v, e, a){
  l <- matrix(0, length(y), 5)
  l[, 1] <- (y <= v) # VaR violations
  l[, 2] <- (a - (y <= v)) * (y - v) # loss function of Clements et al. (2008)
  l[, 3] <- (y <= v) - a
  l[, 4] <- (1 / a) * (y <= v) * (y / e) - 1 # standardized residuals
  l[, 5] <- -(1 / (a * e)) * (y <= v) * (v - y) +  (v / e) + log(-e) - 1 # FZ0
  return(l)
}

# Dynamic quantile test of Engle and Manganelli (2004).
# 2026-10-30 (IJF R4): corrected. Hit_t = I(y_t <= VaR_t) - a is regressed on a constant, q lagged hits, and VaR_t;
# the statistic is b'X'Xb / (a(1-a)), chi-square with q+2 degrees of freedom under the null (EM 2004, eq. 9).
# The previous version divided by the residual variance times the number of regressors (s2 * k) and used one lag,
# which shrank the statistic by roughly a factor of eight at the 1% level; its p-values were reported in R1-R3.
# Arguments: y1 = centred hits (Hit_t - a), y2 = VaR forecasts, a = nominal level, q = number of hit lags (default 4).
DQtest <- function(y1, y2, a, q = 4){
  n <- length(y1)
  idx <- (q + 1):n
  x <- cbind(1, sapply(1:q, function(j) y1[idx - j]), y2[idx])
  y <- y1[idx]
  # b'X'Xb equals the sum of squared fitted values, so lm.fit is used: it handles the rank-deficient X'X that arises
  # when a model has no violation in the lagged range (the lagged-hit columns then equal the constant times -a).
  fit <- lm.fit(x, y)
  k <- fit$rank
  dq <- sum(fit$fitted.values^2) / (a * (1 - a))
  pval <- 1 - pchisq(dq, k)
  return(pval)
}

# Regression-based check on the standardized ES residuals (one lag, residual-variance denominator).
# This is the statistic used for both columns before 2026-10-30; it is kept only for the "DES" column, which is not reported.
DQtest_resid <- function(y1, y2){
  n <- length(y1)
  y <- y1[2:n]
  x <- cbind(1, y1[1:(n-1)], y2[2:n])
  k <- ncol(x)
  fit <- lm(y ~ x - 1)
  b <- matrix(fit$coefficients, k, 1)
  e <- unname(fit$residuals)
  s2 <- sum(e^2) / (n - k)
  dq <- t(b) %*% (t(x) %*% x) %*% b / (s2 * k)
  pval <- 1 - pchisq(dq, k)
  return(pval)
}

# Unconditional coverage test of Kupiec (1995): likelihood ratio of the violation frequency against the nominal level,
# chi-square with one degree of freedom under the null. Added 2026-10-30 (IJF R4). hit = violation indicators, a = nominal level.
UCtest <- function(hit, a){
  n <- length(hit); x <- sum(hit); ph <- x / n
  if (x == 0 || x == n) {
    lr <- -2 * (x * log(a) + (n - x) * log(1 - a))
  } else {
    lr <- -2 * (x * log(a) + (n - x) * log(1 - a)) + 2 * (x * log(ph) + (n - x) * log(1 - ph))
  }
  pval <- 1 - pchisq(lr, 1)
  return(pval)
}