# IJF R2 (2026-08-19): consolidate the psi-free DJIA COVID-period RSV-FS-{SN,ST}
# raw per-window draws (h_pred_s_<t>.csv / y_pred_s_<t>.csv) into per-model
# var_/es_/volf_ summary files. DJIA counterpart of
# merge_prediction_results_fs_psi_n225_2020.R.
#
# The FS psi-free runs for the DJIA COVID sample were executed on the Sakura
# server (RSV_FS_{SN,ST}_psi_run<t>.ox) and the per-window draws were retrieved
# into Ox/DJIA_prediction_2020/RSV_FS_{SN,ST}_psi/. They were never aggregated,
# so the psi sensitivity analysis covered only the AZ family.
#
# Column layout follows merge_prediction_results_fs_djia_2020.R: volf carries
# mean(h), median(h), mean(exp h), median(exp h); var/es carry only the 1% and
# 5% quantiles, whereas the AZ files carry four quantiles with 1% in Var2 and
# 5% in Var3 (see PROGRESS.md 36). The consumer script shifts the FS names
# accordingly, so we keep the FS-native two-column form here.

rm(list = ls())

models <- c("rsv_fs_sn_psi" = "RSV_FS_SN_psi",
            "rsv_fs_st_psi" = "RSV_FS_ST_psi")
n_pred <- 253
probs <- c(0.01, 0.05)
base_dir <- "../Ox/DJIA_prediction_2020/"
output_dir <- "../Ox/DJIA_prediction_2020/"

for (model in names(models)) {
  input_dir <- file.path(base_dir, models[[model]])
  cat("\nProcessing model:", model, "(", input_dir, ")\n")

  var_results <- matrix(NA_real_, nrow = n_pred, ncol = length(probs))
  es_results <- matrix(NA_real_, nrow = n_pred, ncol = length(probs))
  volf_results <- matrix(NA_real_, nrow = n_pred, ncol = 4)

  for (x in 0:(n_pred - 1)) {
    h_file <- file.path(input_dir, paste0("h_pred_s_", x, ".csv"))
    y_file <- file.path(input_dir, paste0("y_pred_s_", x, ".csv"))
    if (!file.exists(h_file) || !file.exists(y_file)) {
      stop(paste("missing draws for window", x, "in", model))
    }
    h_samples <- read.csv(h_file, header = TRUE)[, 2]
    y_samples <- read.csv(y_file, header = TRUE)[, 2]

    vvar <- quantile(y_samples, probs = probs)
    var_results[x + 1, ] <- vvar
    for (i in seq_along(probs)) {
      es_results[x + 1, i] <- mean(y_samples[y_samples < vvar[i]])
    }
    volf_results[x + 1, 1] <- mean(h_samples)
    volf_results[x + 1, 2] <- median(h_samples)
    volf_results[x + 1, 3] <- mean(exp(h_samples))
    volf_results[x + 1, 4] <- median(exp(h_samples))
  }

  stopifnot(!anyNA(var_results), !anyNA(es_results), !anyNA(volf_results))

  colnames(var_results) <- paste0("VaR_", probs * 100, "p")
  colnames(es_results) <- paste0("ES_", probs * 100, "p")
  colnames(volf_results) <- c("Mean", "Median", "MeanExp", "MedianExp")

  write.csv(var_results, file.path(output_dir, paste0("var_", model, ".csv")), row.names = FALSE)
  write.csv(es_results, file.path(output_dir, paste0("es_", model, ".csv")), row.names = FALSE)
  write.csv(volf_results, file.path(output_dir, paste0("volf_", model, ".csv")), row.names = FALSE)
  cat("  saved var_/es_/volf_", model, ".csv (", n_pred, "rows )\n")
}
cat("\nAll processes completed.\n")
