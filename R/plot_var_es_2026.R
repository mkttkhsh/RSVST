rm(list = ls())

# packages
library(tidyverse)
library(ggrepel)
library(readxl)
library(reshape2)
library(MCS)
library(sandwich)

# load function
source("function_forecast_evaluation.R")
source("function_plot_forecast_2026.R")

# set path
path.data <- "../Data/"
path.result <- "./prediction_results_2026/"
path.figure <- "./figure2026/"

# 2026-09-28 (IJF R3, PROGRESS.md §81): RSVST_HL=r3corrected switches to the corrected Hansen-Lunde proxies
# (Data/*_hl_roll*_r3corrected.csv) and to separate output directories (*_r3corrected/); old results stay untouched.
hl.tag <- Sys.getenv("RSVST_HL", "")
if (nzchar(hl.tag)) {
  if (exists("true.sets")) true.sets <- paste0(true.sets, "_", hl.tag)
  for (.p in intersect(c("path.result", "path.out", "path.figure"), ls())) {
    assign(.p, sub("/?$", paste0("_", hl.tag, "/"), get(.p))); dir.create(get(.p), showWarnings = FALSE, recursive = TRUE)
  }
  cat("RSVST_HL =", hl.tag, "\n")
}

# 2026-10-01 (IJF R3, PROGRESS.md §86): RSVST_AZ=<tag> (dualfixed | dualorig) switches the AZ-family forecasts (SV/RSV-AZ-SN/ST and the
# psi-free RSV-AZ variants) to the corrected-recomputation aggregates (Ox/az_recomputation/aggregated/<tag>/ (override with RSVST_AZ_DIR), read in the
# merge scripts) and appends _<tag> to the result/figure directories (after the RSVST_HL suffix). Old results stay untouched.
az.tag <- Sys.getenv("RSVST_AZ", "")
if (nzchar(az.tag)) {
  for (.p in intersect(c("path.result", "path.out", "path.figure"), ls())) {
    assign(.p, sub("/?$", paste0("_", az.tag, "/"), get(.p))); dir.create(get(.p), showWarnings = FALSE, recursive = TRUE)
  }
  cat("RSVST_AZ =", az.tag, "\n")
}

# define sets of data
data.sets <- c("DJIA", "N225")
data.n <- length(data.sets)

# set model names
model.names <- c(paste0("SV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("RSV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("EGARCH-", c("N", "T")),
                 paste0("REGARCH-", c("N", "T")))
model.n <- length(model.names)

# set colors and linetypes
cols <- ggColorHue(model.n + 1)
cols.n <- length(cols)
ltys <- c("dotted", rep("dashed", 7), rep("solid", 7), rep("longdash", 2), rep("twodash", 2))
ltys.n <- length(ltys)

# set range of data
range.var.sets <- rbind(c(-6, 6), c(-9, 9))

# set sets of target rates
alpha.sets <- c("1p", "5p")

# select data
data.i <- 2
data.name <- data.sets[data.i]
range.var <- c(range.var.sets[data.i, ])

# read VaR forecasts
f.var1p <- read_csv(paste0(path.result, data.name, "_forecast_var1p.csv"))
f.var5p <- read_csv(paste0(path.result, data.name, "_forecast_var5p.csv"))

# read ES forecasts
f.es1p <- read_csv(paste0(path.result, data.name, "_forecast_es1p.csv"))
f.es5p <- read_csv(paste0(path.result, data.name, "_forecast_es5p.csv"))

# set dates
f.date <- as.Date(as.character(f.var1p[["date"]]), format = '%Y%m%d')

# number of observations
n <- length(f.date)

#### Computations ####

##### compute violations and various lossfunctions #####
# loss function of Clements et al. (2008), 
# standardized residuals for VaR and ES, and FZ0 loss function of Patton et al. (2019) 
loss1p <- array(NA, c(model.n, n, 5))
loss5p <- array(NA, c(model.n, n, 5))
for (i in 1:model.n) {
  loss1p[i, , ] <- loss_var_es_fz0(f.var1p[[2]], f.var1p[[i + 2]], f.es1p[[i + 2]], 0.01)
  loss5p[i, , ] <- loss_var_es_fz0(f.var5p[[2]], f.var5p[[i + 2]], f.es5p[[i + 2]], 0.05)
}

# extract FZ0 loss
loss.fz0.1p <- t(loss1p[, , 5])
loss.fz0.5p <- t(loss5p[, , 5])


#### Plots ####

##### plot VaR and ES (1% and 5%) forecasts of each model #####
# data.label <- c("Return", "", "ES 1%", "VaR Violation 1%", "", "ES 5%", "VaR Violation 5%")
# for (i in 1:model.n) {
#   # 1%
#   violation1p <- f.var1p[[2]]
#   violation1p[!loss1p[i, ,1]] <- NA
#   df.var.es.1p <- data.frame(f.date, f.var1p[[2]], f.var1p[[i + 2]], f.es1p[[i + 2]], violation1p)
#   names(df.var.es.1p) <- c("date", "Return", "VaR (1%)", "ES (1%)", "Violation (1%)")
#   plot_var_es(
#     f.var.es = df.var.es.1p,
#     out.name = paste0(path.figure, data.name, "_forecast_var_es_1p_", gsub("-", "_", tolower(model.names[i])), ".png"),
#     size = 40,
#     i.legend = TRUE,
#     # legend.xy = ...  <-- 削除します（デフォルトの左上が適用されます）
#     range.y = range.var
#   )
#   # 5%
#   violation5p <- f.var5p[[2]]
#   violation5p[!loss5p[i, ,1]] <- NA
#   df.var.es.5p <- data.frame(f.date, f.var5p[[2]], f.var5p[[i + 2]], f.es5p[[i + 2]], violation5p)
#   names(df.var.es.5p) <- c("date", "Return", "VaR (5%)", "ES (5%)", "Violation (5%)")
#   plot_var_es(
#     f.var.es = df.var.es.5p,
#     out.name = paste0(path.figure, data.name, "_forecast_var_es_5p_", gsub("-", "_", tolower(model.names[i])), ".png"),
#     size = 40,
#     i.legend = TRUE,
#     # legend.xy = ...  <-- 削除します
#     range.y = range.var
#   )
# }

##### plot cumulative loss difference for FZ0 (1%) #####

# データフレーム作成
loss.fz0.cld.1p <- data.frame(f.date, apply(-loss.fz0.1p[, 2:model.n] + loss.fz0.1p[, 1], 2, cumsum))
names(loss.fz0.cld.1p) <- c("date", model.names[2:model.n])

# 新しい関数でプロット (Faceted)
plot_fz0_cld_faceted(
  f.cld = loss.fz0.cld.1p, 
  out.name = paste0(path.figure, data.name, "_forecast_fz0_cld_1p_faceted.png"), 
  size = 20, 
  cols = cols, 
  i.legend = TRUE
)

plot_fz0_cld_faceted_outer_legend(
  f.cld = loss.fz0.cld.1p, 
  out.name = paste0(path.figure, data.name, "_forecast_fz0_cld_1p_faceted_outer.png"), 
  size = 20, 
  cols = cols, 
  i.legend = TRUE
)  

##### plot cumulative loss difference for FZ0 (5%) #####

# データフレーム作成
loss.fz0.cld.5p <- data.frame(f.date, apply(-loss.fz0.5p[, 2:model.n] + loss.fz0.5p[, 1], 2, cumsum))
names(loss.fz0.cld.5p) <- c("date", model.names[2:model.n])

# 新しい関数でプロット (Faceted)
plot_fz0_cld_faceted(
  f.cld = loss.fz0.cld.5p, 
  out.name = paste0(path.figure, data.name, "_forecast_fz0_cld_5p_faceted.png"), 
  size = 20, 
  cols = cols, 
  i.legend = TRUE
)

plot_fz0_cld_faceted_outer_legend(
  f.cld = loss.fz0.cld.5p, 
  out.name = paste0(path.figure, data.name, "_forecast_fz0_cld_5p_faceted_outer.png"), 
  size = 20, 
  cols = cols, 
  i.legend = TRUE
)  
