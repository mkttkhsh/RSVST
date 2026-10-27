# function to implement Ljung-Box test adjusted for heteroskedasticity following Diebold (1988)
LB_Diebold <- function(x, m){
  n <- length(x)
  rho <- acf(x, lag.max = m, plot = FALSE)$acf[2:(m+1)]
  x2 <- (x - mean(x)) ^ 2
  gam <- acf(x2, lag.max = m, type = "covariance", plot = FALSE)$acf[2:(m+1)]
  sig4 <- var(x) ^ 2
  LBm <- 0
  for (k in 1:m) {
    LBm <- LBm + (sig4 / (sig4 + gam[k])) * rho[k]^2 / (n - k)
  }
  LBm <- LBm * n * (n + 2)
  pval <- 1 - pchisq(LBm, m)
  return(pval)
}