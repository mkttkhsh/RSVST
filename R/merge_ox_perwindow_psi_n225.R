# 2026 revision: per-window -> per-model consolidation for the psi-free
# N225 COVID-period RSV forecasts, done in R instead of Ox (no local Ox
# install) -- replicates merge_prediction_results_az_2020.ox / _gh_2020.ox
# for the "_psi" models only (the non-psi models were already consolidated
# when copied in from Ox2022/). Reads the 242 per-window files (var_/es_/
# volf_/h1_y_<model>_<t>.csv) and writes one consolidated var_/es_/volf_
# <model>.csv per model, matching the column layout the R evaluation
# scripts (merge_prediction_results_n225_2020.R) expect.

rm(list = ls())
library(readr)

n_pred <- 242
in_dir <- "../Ox/N225_prediction_2020/"

models <- c("rsv_az_n_psi", "rsv_az_sn_psi", "rsv_az_st_psi", "rsv_az_t_psi", "rsv_gh_st_psi")

for (model in models) {
  cat("Processing model:", model, "\n")

  mvar <- matrix(NA_real_, nrow = n_pred, ncol = 4)
  mes  <- matrix(NA_real_, nrow = n_pred, ncol = 4)
  mvolf <- matrix(NA_real_, nrow = n_pred, ncol = 4) # mean(h), median(h), mean(exp(h)), median(exp(h))

  for (t in 0:(n_pred - 1)) {
    stemp <- paste0(model, "_", t)

    f.var <- read_csv(paste0(in_dir, "var_", stemp, ".csv"), col_types = cols(), show_col_types = FALSE)
    f.es  <- read_csv(paste0(in_dir, "es_",  stemp, ".csv"), col_types = cols(), show_col_types = FALSE)
    f.volf <- read_csv(paste0(in_dir, "volf_", stemp, ".csv"), col_types = cols(), show_col_types = FALSE)
    f.h1y <- read_csv(paste0(in_dir, "h1_y_", stemp, ".csv"), col_types = cols(), show_col_types = FALSE)

    mvar[t + 1, ] <- as.numeric(f.var[1, c("Var1", "Var2", "Var3", "Var4")])
    mes[t + 1, ]  <- as.numeric(f.es[1,  c("Var1", "Var2", "Var3", "Var4")])

    vh1 <- f.h1y$Var1 # simulated one-day-ahead log-volatility draws
    mvolf[t + 1, 1] <- as.numeric(f.volf[1, "Var1"])      # mean(h), from the class's own savemat
    mvolf[t + 1, 2] <- median(vh1)                        # median(h), recomputed from draws
    mvolf[t + 1, 3] <- mean(exp(vh1))                     # mean(exp(h))
    mvolf[t + 1, 4] <- median(exp(vh1))                   # median(exp(h))
  }

  colnames(mvar) <- paste0("Var", 1:4)
  colnames(mes) <- paste0("Var", 1:4)
  colnames(mvolf) <- paste0("Var", 1:4)

  write_csv(as.data.frame(mvar), paste0(in_dir, "var_", model, ".csv"))
  write_csv(as.data.frame(mes), paste0(in_dir, "es_", model, ".csv"))
  write_csv(as.data.frame(mvolf), paste0(in_dir, "volf_", model, ".csv"))

  cat("  saved var_/es_/volf_", model, ".csv (", n_pred, "rows)\n")
}

cat("\nDone.\n")
