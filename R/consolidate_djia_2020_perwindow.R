# IJF R2 (2026): consolidate per-window var_/es_/volf_<model>_<t>.csv output
# (t=0..252) from the DJIA COVID-period (2012-2020) rolling estimation into
# single per-model var_/es_/volf_<model>.csv files, replacing the Ox-level
# merge_prediction_results_{az,gh}_2020.ox step (no local Ox install).
#
# Verified against the already-consolidated N225 files (same pipeline, run
# 2026-07-11 via the real Ox merge script) before writing this:
#   AZ family: raw per-window volf_<model>_<t>.csv has 2 columns
#     (Var1 = mean(h), Var2 = median(h), both computed by the Ox model class
#     itself). Var2 numerically matches median(h1_y draws) exactly, confirming
#     these are already final. mean(exp(h)) / median(exp(h)) are NOT in the
#     raw file and must be computed from the h1_y_<model>_<t>.csv draws
#     (column 0 = h draws). Consolidated column order: mean(h), median(h),
#     mean(exp(h)), median(exp(h)).
#   GH family (rsv_gh_st, rsv_gh_st_psi, sv_gh_st): raw per-window
#     volf_<model>_<t>.csv has 4 columns, of which only the first two are
#     used: Var1 = mean(h), Var2 = mean(exp(h)) (cross-checked against
#     mean(exp(h1_y draws)) -- exact match). The original Ox merge script has
#     an indexing quirk (a loop variable meant to flag "not the first model"
#     is always 0 when the model list has a single entry) that leaves the
#     median(h) and median(exp(h)) columns as placeholders (-9999 and 0);
#     this exact placeholder pattern is present in the already-published
#     N225 consolidated files, so it is reproduced here for consistency
#     rather than "fixed" -- downstream evaluation only uses the mean(exp(h))
#     column (Var3).
#
# var_/es_ files: each per-window file is already a single row in the target
# format; consolidation is a plain row bind across t=0..252, no computation.

library(readr)
library(dplyr)

path.in <- "../Ox/DJIA_prediction_2020/"
nprd <- 253

gh_models <- c("rsv_gh_st", "rsv_gh_st_psi", "sv_gh_st")
az_models <- c(
  "rsv_az_n", "rsv_az_sn", "rsv_az_st", "rsv_az_t",
  "rsv_az_n_psi", "rsv_az_sn_psi", "rsv_az_st_psi", "rsv_az_t_psi",
  "sv_az_n", "sv_az_sn", "sv_az_st", "sv_az_t"
)
all_models <- c(az_models, gh_models)
stopifnot(length(all_models) == 15)

# Optional model filter, so that a repaired model can be re-consolidated on its
# own without touching the others (and without tripping over models whose
# per-window output is incomplete):
#   Rscript consolidate_djia_2020_perwindow.R sv_az_n sv_az_st
# With no arguments every model is processed, as before.
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) > 0) {
  .unknown <- setdiff(.args, all_models)
  if (length(.unknown) > 0) {
    stop("unknown model(s): ", paste(.unknown, collapse = ", "))
  }
  all_models <- .args
}

read_raw <- function(prefix, model, t) {
  path <- paste0(path.in, prefix, "_", model, "_", t, ".csv")
  as.numeric(read_csv(path, col_names = FALSE, skip = 1, show_col_types = FALSE)[1, -1])
}

for (model in all_models) {
  cat(sprintf("=== %s ===\n", model))

  mvar <- NULL
  mes <- NULL
  mvolf <- matrix(NA_real_, nprd, 4)

  is_gh <- model %in% gh_models

  for (t in 0:(nprd - 1)) {
    v <- read_raw("var", model, t)
    e <- read_raw("es", model, t)
    if (t == 0) {
      mvar <- matrix(NA_real_, nprd, length(v))
      mes <- matrix(NA_real_, nprd, length(e))
    }
    mvar[t + 1, ] <- v
    mes[t + 1, ] <- e

    vf <- read_raw("volf", model, t)
    if (is_gh) {
      mean_h <- vf[1]
      mean_exp_h <- vf[2]
      median_h <- -9999
      median_exp_h <- 0
    } else {
      mean_h <- vf[1]
      median_h <- vf[2]
      h1 <- read_csv(paste0(path.in, "h1_y_", model, "_", t, ".csv"), show_col_types = FALSE)[[2]]
      mean_exp_h <- mean(exp(h1))
      median_exp_h <- median(exp(h1))
    }
    mvolf[t + 1, ] <- c(mean_h, median_h, mean_exp_h, median_exp_h)
  }

  colnames(mvar) <- paste0("Var", 1:ncol(mvar))
  colnames(mes) <- paste0("Var", 1:ncol(mes))
  colnames(mvolf) <- paste0("Var", 1:4)

  write_csv(as.data.frame(mvar), paste0(path.in, "var_", model, ".csv"))
  write_csv(as.data.frame(mes), paste0(path.in, "es_", model, ".csv"))
  write_csv(as.data.frame(mvolf), paste0(path.in, "volf_", model, ".csv"))

  cat(sprintf("  wrote var(%dx%d) es(%dx%d) volf(%dx%d)\n",
              nrow(mvar), ncol(mvar), nrow(mes), ncol(mes), nrow(mvolf), ncol(mvolf)))
}

cat("Done.\n")
