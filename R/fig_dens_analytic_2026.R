# Figure 1 of the main text (panels (a) GH, (b) Azzalini, (c) FS) and the lower panels of the
# supplementary density figures, all drawn from ANALYTIC densities (no simulation). 2026-10-14.
#
# - GH skew-t: dens_gh() in function_dens_skew_t.R (standardized), as in fig_dens_ghst_2026.R.
# - FS skew-t: dens_fs() in function_dens_skew_t.R (standardized), as in fig_dens_fs_2026.R.
# - Azzalini skew-t: NEW dens_az_eq6() below, the density of eps_t in equation (6) of the main text,
#       eps = (delta (z0 - c) + sqrt(1 - delta^2) z) / sqrt(1 - c^2 delta^2) * sqrt(lambda / mu_lambda),
#       z0 ~ TN(0,inf)(0,1), z ~ N(0,1), lambda ~ IG(nu/2, nu/2), mu_lambda = nu/(nu-2), c = sqrt(2/pi).
#   This differs from dens_az() (Azzalini's ST(0,1,alpha,nu) standardized after mixing): in (6) the
#   centering by c is applied BEFORE the scale mixture, which dilutes the skewness (e.g. delta = -0.6
#   gives skewness -0.08, delta = -0.9 gives -0.52). The earlier simulation-based figure
#   (fig_dens_azst_2026.R) drew (6) correctly; this script reproduces it analytically.
#
# Run from R/: Rscript fig_dens_analytic_2026.R
# Output: figure2026/fig/fig_dens_{ghst_beta,azst_delta,fsst_gamma}.png        (main text, Figure 1, (a)(b)(c))
#         figure2026/fig/fig_dens_{ghst_beta,azst_delta,fsst_gamma}_suppl.png  (supplement, upper panels, (i))
#         figure2026/fig/fig_dens_{ghst_nu,azst_nu,fsst_nu}.png                (supplement, lower panels, (ii))
rm(list = ls(all = TRUE))
source("function_dens_skew_t.R")

# ---- density of eq. (6): integrate over q = nu / lambda ~ chi^2_nu --------------------------------
dens_az_eq6 <- function(x, nu, delta) {
  cc <- sqrt(2 / pi)
  s <- sqrt(1 - cc^2 * delta^2)
  alpha <- delta / sqrt(1 - delta^2)
  f_w <- function(w) { u <- s * w + delta * cc; s * 2 * dnorm(u) * pnorm(alpha * u) }
  one <- function(xx) {
    integrate(function(q) { r <- sqrt((nu - 2) / q); f_w(xx / r) / r * dchisq(q, nu) },
              lower = 0, upper = Inf, rel.tol = 1e-8, subdivisions = 400L)$value
  }
  vapply(x, one, numeric(1))
}

# ---- common plotting ------------------------------------------------------------------------------
x.max <- 6; x.min <- -6; y.max <- 0.7; y.min <- 0
dx <- 0.005
x <- seq(x.min, x.max, by = dx)
dir.create("figure2026/fig", showWarnings = FALSE, recursive = TRUE)

check <- function(tag, f) {
  m0 <- sum(f) * dx; m1 <- sum(x * f) * dx; m2 <- sum(x^2 * f) * dx
  v <- m2 - m1^2; m3 <- sum((x - m1)^3 * f) * dx
  cat(sprintf("%-28s mass %.4f mean %7.4f var %.4f skew %7.4f max %.4f at x = %6.3f\n",
              tag, m0, m1, v, m3 / v^1.5, max(f), x[which.max(f)]))
}

draw <- function(file, f, labels, title) {
  png(file, width = 1600, height = 900)
  plot(x, f[, 1], type = "l", col = "gray", lwd = 6, ylab = "", xlab = "", main = "",
       xlim = c(x.min, x.max), ylim = c(y.min, y.max), cex.axis = 2.6)
  lines(x, f[, 2], lty = 2, lwd = 5)
  lines(x, f[, 3], lwd = 3)
  axis(side = 1, at = seq(x.min, x.max, 1), tck = 1.0, lty = "dotted", lwd = 1, col = "gray", labels = NA)
  axis(side = 2, at = seq(y.min, y.max, 0.1), tck = 1.0, lty = "dotted", lwd = 1, col = "gray", labels = NA)
  legend(-6.2, 0.55, legend = labels, lty = c(1, 2, 1), col = c("gray", "black", "black"),
         lwd = c(6, 5, 3), cex = 4.5, bg = "white")
  mtext(side = 3, text = title, line = -3.8, adj = 0.03, cex = 4.5)
  dev.off()
}
lab <- function(sym, vals, fmt) lapply(vals, function(v) bquote(.(sym) == .(sprintf(fmt, v))))
as_expr <- function(l) do.call(expression, l)

# ---- (a) GH skew-t: beta = 0, -1, -2, nu = 10; supplement: nu = 24, 16, 8, beta = -1 --------------
nu <- 10; beta <- c(0, -1, -2)
f <- sapply(beta, function(b) dens_gh(x, nu, b, standardize = TRUE))
for (i in seq_along(beta)) check(sprintf("GH beta = %4.1f nu = %2d", beta[i], nu), f[, i])
draw("figure2026/fig/fig_dens_ghst_beta.png", f,
     as_expr(lab(quote(beta), beta, "%2.0f")),
     expression(paste("(a) GH skew-", italic(t), ", ", nu, " = 10")))
draw("figure2026/fig/fig_dens_ghst_beta_suppl.png", f,
     as_expr(lab(quote(beta), beta, "%2.0f")),
     expression(paste("(i) GH skew-", italic(t), ", ", nu, " = 10")))

nus <- c(24, 16, 8); beta <- -1
f <- sapply(nus, function(n) dens_gh(x, n, beta, standardize = TRUE))
for (i in seq_along(nus)) check(sprintf("GH beta = %4.1f nu = %2d", beta, nus[i]), f[, i])
draw("figure2026/fig/fig_dens_ghst_nu.png", f,
     as_expr(lab(quote(nu), nus, "%2.0f")),
     expression(paste("(ii) GH skew-", italic(t), ", ", beta, " = -1")))

# ---- (b) Azzalini skew-t, eq. (6): delta = 0, -0.8, -0.95, nu = 10 (2026-10-14: was 0, -0.6, -0.9; delta = -0.6 is visually indistinguishable from 0); supplement: nu = 15, 10, 5 -----
nu <- 10; delta <- c(0, -0.8, -0.95)
f <- sapply(delta, function(d) dens_az_eq6(x, nu, d))
for (i in seq_along(delta)) check(sprintf("AZ delta = %5.2f nu = %2d", delta[i], nu), f[, i])
draw("figure2026/fig/fig_dens_azst_delta.png", f,
     as_expr(lab(quote(delta), delta, "%5.2f")),
     expression(paste("(b) Azzalini skew-", italic(t), ", ", nu, " = 10")))
draw("figure2026/fig/fig_dens_azst_delta_suppl.png", f,
     as_expr(lab(quote(delta), delta, "%5.2f")),
     expression(paste("(i) Azzalini skew-", italic(t), ", ", nu, " = 10")))

nus <- c(15, 10, 5); delta <- -0.9
f <- sapply(nus, function(n) dens_az_eq6(x, n, delta))
for (i in seq_along(nus)) check(sprintf("AZ delta = %5.2f nu = %2d", delta, nus[i]), f[, i])
draw("figure2026/fig/fig_dens_azst_nu.png", f,
     as_expr(lab(quote(nu), nus, "%2.0f")),
     expression(paste("(ii) Azzalini skew-", italic(t), ", ", delta, " = -0.9")))

# ---- (c) FS skew-t: gamma = 1, 0.8, 0.5, nu = 10; supplement: nu = 20, 15, 10, gamma = 0.5 --------
nu <- 10; gam <- c(1, 0.8, 0.5)
f <- sapply(gam, function(g) dens_fs(x, nu, g, standardize = TRUE))
for (i in seq_along(gam)) check(sprintf("FS gamma = %4.1f nu = %2d", gam[i], nu), f[, i])
draw("figure2026/fig/fig_dens_fsst_gamma.png", f,
     as_expr(lab(quote(gamma), gam, "%3.1f")),
     expression(paste("(c) FS skew-", italic(t), ", ", nu, " = 10")))
draw("figure2026/fig/fig_dens_fsst_gamma_suppl.png", f,
     as_expr(lab(quote(gamma), gam, "%3.1f")),
     expression(paste("(i) FS skew-", italic(t), ", ", nu, " = 10")))

nus <- c(20, 15, 10); gam <- 0.5
f <- sapply(nus, function(n) dens_fs(x, n, gam, standardize = TRUE))
for (i in seq_along(nus)) check(sprintf("FS gamma = %4.1f nu = %2d", gam, nus[i]), f[, i])
draw("figure2026/fig/fig_dens_fsst_nu.png", f,
     as_expr(lab(quote(nu), nus, "%2.0f")),
     expression(paste("(ii) FS skew-", italic(t), ", ", gamma, " = 0.5")))
