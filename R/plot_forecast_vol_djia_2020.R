rm(list = ls())

#### Settings ####

# packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(reshape2)
library(forcats)

# function
source("function_plot_forecast_2026.R")

# set path
path.data <- "../Data/"
path.result <- "./prediction_results_djia_2020/"
path.figure <- "./figure_djia_2020/"

# 2026-09-28 (October 2026 revision, PROGRESS.md §81): RSVST_HL=r3corrected switches to the corrected Hansen-Lunde proxies
# (Data/*_hl_roll*_r3corrected.csv) and to separate output directories (*_r3corrected/); old results stay untouched.
hl.tag <- Sys.getenv("RSVST_HL", "")
if (nzchar(hl.tag)) {
  if (exists("true.sets")) true.sets <- paste0(true.sets, "_", hl.tag)
  for (.p in intersect(c("path.result", "path.out", "path.figure"), ls())) {
    assign(.p, sub("/?$", paste0("_", hl.tag, "/"), get(.p))); dir.create(get(.p), showWarnings = FALSE, recursive = TRUE)
  }
  cat("RSVST_HL =", hl.tag, "\n")
}

# 2026-10-01 (October 2026 revision, PROGRESS.md §86): RSVST_AZ=<tag> (dualfixed | dualorig) switches the AZ-family forecasts (SV/RSV-AZ-SN/ST and the
# psi-free RSV-AZ variants) to the corrected-recomputation aggregates (Ox/az_recomputation/aggregated/<tag>/ (override with RSVST_AZ_DIR), read in the
# merge scripts) and appends _<tag> to the result/figure directories (after the RSVST_HL suffix). Old results stay untouched.
az.tag <- Sys.getenv("RSVST_AZ", "")
if (nzchar(az.tag)) {
  for (.p in intersect(c("path.result", "path.out", "path.figure"), ls())) {
    assign(.p, sub("/?$", paste0("_", az.tag, "/"), get(.p))); dir.create(get(.p), showWarnings = FALSE, recursive = TRUE)
  }
  cat("RSVST_AZ =", az.tag, "\n")
}

# define data
data.name <- "DJIA"

# define volatility proxies
rv.sets <- c("RV5", "RKth2", "BV", "MedRV")
rv.n <- length(rv.sets)

# set model names
# revision (2026-09-24): SV-GH-ST restored (forecast-step bug fixed, re-estimated; PROGRESS.md section 74).
model.names <- c(paste0("SV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("RSV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST")),
                 paste0("EGARCH-", c("N", "T")),
                 paste0("REGARCH-", c("N", "T")))
model.n <- length(model.names)

#### Set Prediction Results ####

# read volatility forecasts
f.vol <- read_csv(paste0(path.result, data.name, "_forecast_vol.csv"))
f.vol.med <- read_csv(paste0(path.result, data.name, "_forecast_vol_med.csv"))
# 2026-09-24 (R3): predictive mean for all models; the median substitution for SV models was removed (AE comment 2).

# convert yyyymmdd to dates
yyyymmdd <- f.vol[["date"]]
f.vol[["date"]] <- as.Date(as.character(f.vol[["date"]]),format = '%Y%m%d')

# log-scale
f.vol.log <- cbind(f.vol[,1], log(f.vol[,-1]))

# ---- 必要な列だけ抽出して長形式に変換 ----
vol_ts_long <- f.vol.log %>%
  select(date, all_of(rv.sets)) %>%
  pivot_longer(cols = all_of(rv.sets), names_to = "proxy_type", values_to = "volatility") 

vol_ts_long$proxy_type <- factor(vol_ts_long$proxy_type, levels = rv.sets)

vol_ts_long <- vol_ts_long %>%
  mutate(proxy_type = fct_recode(proxy_type,
                                 "RV5" = "RV5",
                                 "RK" = "RKth2",
                                 "BV" = "BV",
                                 "Med" = "MedRV"))

# (2026-09-24: the 17-model re-declaration that excluded SV-GH-ST was removed)

# ---- select volatility proxies ----
for (rv.i in 1:rv.n) {
  rv.name <- rv.sets[rv.i]
  f.vol.i <- f.vol[c("date", rv.name, model.names)] # select volatility proxy

  # compute QLIKE
  n <- nrow(f.vol.i)
  k <- ncol(f.vol.i) - 2
  loss.mse <- matrix(0, n, k)
  loss.qlike <- matrix(0, n, k)
  for(i in 1:k){
    loss.mse[1:n,i] <- (f.vol.i[[i+2]] - f.vol.i[[2]])^2 / 2
    loss.qlike[1:n,i] <- f.vol.i[[2]]/f.vol.i[[i+2]] - log(f.vol.i[[2]]/f.vol.i[[i+2]]) - 1
  }
  mse <- colMeans(loss.mse)
  qlike <- colMeans(loss.qlike)
  dt.mse.qlike <- data.frame(model.names, mse, qlike)

  ##### time-series plot #####

  # log-scale
  f.vol.i.log <- cbind(f.vol.i[,1], log(f.vol.i[,-1]))

  # ---- ループ内 (for (rv.i in 1:rv.n) { ... }) の記述 ----
  
  # 1. 色の定義
  # モデル数(model.n) + 実績値1つ分 の色を用意
  cols <- ggColorHue(model.n + 1)
  cols[1] <- "gray40" # 1番目（実績値）をグレーに上書き
  
  # 2. 線種の定義
  # 関数内で "solid" に固定しましたが、引数として渡す必要があるため定義します
  ltys <- rep("solid", model.n + 1) 
  
  # 3. 関数呼び出し
  # legend.xy = "inside" と指定する必要はありませんが（関数内で座標指定済み）、
  # 念のためデフォルト引数あるいは明示的に渡します。

  # --- 実行ファイルの修正箇所 ---
  
  plot_vol_faceted(
    f.vol = f.vol.i,                                     # データ
    out.name = paste0(path.figure, data.name, "_", rv.name, "_forecast_vol_faceted.png"), 
    size = 20,                                           # フォントサイズ
    cols = cols,                                         # 色リスト
    # ltys = ltys,                                       # 関数内で自動割当するため削除または無視
    legend.pos.type = "topright",                      # ★ここを修正（位置のタイプを指定）
    i.legend = TRUE,                                     # 凡例を表示するか
    i.title = 0,                                         # タイトル有無
    s.title = ""                                         # タイトル文字列
  )  

  # ログスケール版も同様に
  plot_vol_faceted(
    f.vol = f.vol.i.log, 
    out.name = paste0(path.figure, data.name, "_", rv.name, "_forecast_logvol_faceted.png"), 
    size = 20,                                           # フォントサイズ
    cols = cols,                                         # 色リスト
    # ltys = ltys,                                       # 関数内で自動割当するため削除または無視
    legend.pos.type = "topright",                      # ★ここを修正（位置のタイプを指定）
    i.legend = TRUE,                                     # 凡例を表示するか
    i.title = 0,                                         # タイトル有無
    s.title = ""                                         # タイトル文字列
  )

  plot_vol_faceted_outer_legend(
    f.vol = f.vol.i.log,
    out.name = paste0(path.figure, data.name, "_", rv.name, "_forecast_logvol_faceted_outer.png"), 
    size = 20,                                           # フォントサイズ
    cols = cols,                                         # 色リスト
    i.legend = TRUE
  )
  
  ##### plot cumulative QLIKE #####
  loss.qlike.cld <- data.frame(apply(-loss.qlike[, 2:model.n] + loss.qlike[, 1], 2, cumsum))
  loss.qlike.cld <- data.frame(f.vol.i[1], loss.qlike.cld)
  names(loss.qlike.cld) <- c("date", model.names[2:model.n])

  # 関数呼び出し
  plot_cld_faceted(
    f.cld = loss.qlike.cld, 
    out.name = paste0(path.figure, data.name, "_", rv.name, "_forecast_vol_qlike_cld_faceted.png"),
    size = 20, 
    cols = cols, 
    legend.pos.type = "topleft"
  )  

  plot_cld_faceted_outer_legend(
    f.cld = loss.qlike.cld,
    out.name = paste0(path.figure, data.name, "_", rv.name, "_forecast_vol_qlike_cld_faceted_outer.png"),
    size = 20,
    cols = cols, 
    i.legend = TRUE
  )
  
  # 累積MSE (Fig 9/11) については、Nikkei225の場合はファイル名に _supp を付けるなどして区別
  # suffix <- ""
  # if(data.name == "DJIA") suffix <- "_Supplementary"
  # 
  # plot_vol(loss.qlike.cld, paste0(path.figure, data.name, "_", rv.name, "_forecast_vol_mse_cld", suffix, ".png"),
  #          20, cols[3:cols.n], ltys[3:ltys.n], TRUE, "right")
}
