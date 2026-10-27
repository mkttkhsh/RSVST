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

# function to implement dynamic quantile test of Engle and Manganelli (2004)
DQtest <- function(y1, y2){
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