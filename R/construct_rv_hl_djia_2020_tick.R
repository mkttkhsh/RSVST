# IJF R3 (2026-09-28, PROGRESS.md §80-81): Hansen–Lunde (2005) rolling correction for the evaluation proxies, DJIA COVID sample (2012-2020, 253 forecasts in 2020).
# CORRECTED IMPLEMENTATION. The previous version of this script (output *_hl_roll*.csv, kept for reproducibility) had
#  (1) [main samples] returns taken from the UNFILTERED file (data$rcc, from 2000/1996) instead of the filtered sample,
#  (2) [all samples] the adjusted proxy computed from the proxy of the last day of the estimation window (row t+nest-1)
#      instead of the forecast target day (row t+nest), i.e. the evaluation proxy was lagged by one trading day,
#  (3) [all samples] the return window shifted by one day relative to the proxy window.
# Now, for each forecast origin t (estimation window = rows t..t+nest-1 of the sample, target day = row t+nest):
#   c_HL,t = sum_{s in window} (y_s - ybar)^2 / sum_{s in window} x_s   (same rows for returns and proxies)
#   adjusted proxy_{t+nest} = c_HL,t * x_{t+nest}                          (raw proxy of the forecast target day)
# as in the manuscript (Sections 2.2 and 4.3). Output file: DJIA_rcc_rv5_rkth2_bv_medrv_2012_2020_hl_roll_tick_r3corrected.csv
suppressPackageStartupMessages(library(dplyr))
data <- read.csv("../Data/DJIA_rcc_rv5_rkth2_bv_medrv_2012_2020_tick.csv")
filtered_data <- data  # the file already spans exactly the sample (no filtering)
nobs <- nrow(filtered_data); nprd <- 253; nest <- nobs - nprd
nrv <- ncol(filtered_data) - 2
px <- as.matrix(filtered_data[, 3:(2 + nrv)])
c_hl <- matrix(0, nprd, nrv); rv_hl <- matrix(0, nprd, nrv)
for (t in 1:nprd) {
  idx <- t:(t + nest - 1)                       # estimation window (same rows for returns and proxies)
  r <- filtered_data$rcc[idx]
  c_hl[t, ] <- sum((r - mean(r))^2) / colSums(px[idx, , drop = FALSE])
  rv_hl[t, ] <- c_hl[t, ] * px[t + nest, ]      # raw proxy of the forecast target day
}
colnames(c_hl) <- paste0("c_hl_", colnames(px))
colnames(rv_hl) <- paste0(colnames(px), "_hl")
merged_data <- cbind(filtered_data[(nest + 1):nobs, ], c_hl, rv_hl)
stopifnot(nrow(merged_data) == nprd, all(is.finite(as.matrix(merged_data))))
stopifnot(all(abs(merged_data[, paste0(colnames(px), "_hl")] - merged_data[, paste0("c_hl_", colnames(px))] * merged_data[, colnames(px)]) <= 1e-12 * merged_data[, paste0(colnames(px), "_hl")]))
write.csv(merged_data, "../Data/DJIA_rcc_rv5_rkth2_bv_medrv_2012_2020_hl_roll_tick_r3corrected.csv", row.names = FALSE)
cat(sprintf("%s: %d rows (%d-%d), nest=%d nprd=%d, mean c_hl: %s\n", "DJIA_rcc_rv5_rkth2_bv_medrv_2012_2020_hl_roll_tick_r3corrected.csv", nrow(merged_data), min(merged_data$yyyymmdd), max(merged_data$yyyymmdd), nest, nprd, paste(sprintf("%s=%.3f", colnames(px), colMeans(c_hl)), collapse = " ")))
