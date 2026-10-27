# 2026 revision: consolidate the psi-free N225 COVID-period RSV-FS-{SN,ST}
# raw per-window draws (h_pred_s_<t>.csv / y_pred_s_<t>.csv) into per-model
# var_/es_/volf_ summary files, matching the layout
# merge_prediction_results_n225_2020.R expects. Adapted from
# merge_prediction_results_fs_n225_2020.R (which reads from a per-model
# subdirectory under N225_prediction_2020/); the psi-free FS output instead
# lives directly under output_pred_N225_2020_psi_{st,sn}/ (one directory per
# distribution, no per-model subfolder, since each package only ever ran one
# model).

rm(list = ls())

models <- c("rsv_fs_st_psi" = "output_pred_N225_2020_psi_st",
            "rsv_fs_sn_psi" = "output_pred_N225_2020_psi_sn")
n_pred <- 242
probs <- c(0.01, 0.05)
base_dir <- "../Ox/"
output_dir <- "../Ox/N225_prediction_2020/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (model in names(models)) {
  input_dir <- file.path(base_dir, models[[model]])
  cat("\nProcessing model:", model, "(", input_dir, ")\n")

  var_results <- matrix(NA_real_, nrow = n_pred, ncol = length(probs))
  es_results <- matrix(NA_real_, nrow = n_pred, ncol = length(probs))
  volf_results <- matrix(NA_real_, nrow = n_pred, ncol = 4)

  for (x in 0:(n_pred - 1)) {
    h_file <- file.path(input_dir, paste0("h_pred_s_", x, ".csv"))
    y_file <- file.path(input_dir, paste0("y_pred_s_", x, ".csv"))

    if (file.exists(h_file) && file.exists(y_file)) {
      h_samples <- read.csv(h_file, header = TRUE)[, 2]
      y_samples <- read.csv(y_file, header = TRUE)[, 2]

      vvar <- quantile(y_samples, probs = probs)
      var_results[x + 1, ] <- vvar

      for (i in seq_along(probs)) {
        es_samples <- y_samples[y_samples < vvar[i]]
        es_results[x + 1, i] <- mean(es_samples)
      }

      volf_results[x + 1, 1] <- mean(h_samples)
      volf_results[x + 1, 2] <- median(h_samples)
      volf_results[x + 1, 3] <- mean(exp(h_samples))
      volf_results[x + 1, 4] <- median(exp(h_samples))
    } else {
      warning(paste("File not found for x =", x, "in model", model))
    }
  }

  colnames(var_results) <- paste0("VaR_", probs * 100, "p")
  colnames(es_results) <- paste0("ES_", probs * 100, "p")
  colnames(volf_results) <- c("Mean", "Median", "MeanExp", "MedianExp")

  write.csv(var_results, file.path(output_dir, paste0("var_", model, ".csv")), row.names = FALSE)
  write.csv(es_results, file.path(output_dir, paste0("es_", model, ".csv")), row.names = FALSE)
  write.csv(volf_results, file.path(output_dir, paste0("volf_", model, ".csv")), row.names = FALSE)

  cat("  saved var_/es_/volf_", model, ".csv (", n_pred, "rows)\n")
}

cat("\nAll processes completed.\n")
