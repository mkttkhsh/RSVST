rm(list = ls())

#### Settings ####

# packages
library(readr)
library(ggplot2)
library(readxl)
library(reshape2)
library(MCS)
library(sandwich)
library(ggnewscale)
library(tidyr)   # 2026-09-28: pivot_longer is used below (was only available when tidyverse happened to be attached)
library(tibble)

# function
source("function_plot_forecast_2026.R")

# set path
path.data <- "../Data/"
path.result <- "./prediction_results_2026/"
path.figure <- "./figure2026/"

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

# define sets of data
data.sets <- c("DJIA", "N225")
data.n <- length(data.sets)

# define volatility proxies
rv.sets <- c("RV5", "RKth2", "BV", "MedRV")
rv.n <- length(rv.sets)

# set model names (EG: EGARCH; REG: REGARH)
model.names <- c(paste0("SV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("RSV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("EGARCH-", c("N", "T")),
                 paste0("REGARCH-", c("N", "T")))
model.n <- length(model.names)

###### run one (data, volatility proxy) combination #####
run_output_vol <- function(data.i, rv.i) {
data.name <- data.sets[data.i]
rv.name <- rv.sets[rv.i]
cat(sprintf("\n=== output_forecast_vol_gwtest_2026: %s / %s ===\n", data.name, rv.name))

# read GW test resutls
loss <- read_csv(paste0(path.result, data.name, "_", rv.name, "_forecast_vol_mse_qlike_mcs_gw06test.csv"))
stat <- as.matrix(loss[, 6:(5+model.n)])
prop <- as.matrix(loss[, (6+model.n):(5+2*model.n)])
pval <- as.matrix(loss[, (6+2*model.n):(5+3*model.n)])

##### Output for Latex #####

for (i in 1:model.n) {
  xxx <- model.names[i]
  # xxx <- paste0(xxx, " & ")
  for (j in 1:model.n) {
    if(i == j){
      xxx <- paste0(xxx, " & ")
    }
    if(i != j){
      xxx <- paste0(xxx, " & ")
      if(pval[i,j] <= 0.01){
        xxx <- paste0(xxx, "\\cellcolor{gray!30}{")
      }
      else if(pval[i,j] > 0.01 & pval[i,j] <= 0.05){
        xxx <- paste0(xxx, "\\cellcolor{gray!15}{")
      }
      xxx <- paste0(xxx, "$")
      if(i > j){
        xxx <- paste0(xxx, format(round(stat[i,j], 1), nsmall = 1))
      }
      if(i < j){
        xxx <- paste0(xxx, format(round(prop[i,j]*100, 0), nsmall = 0))
      }
      xxx <- paste0(xxx, "$")
      if(pval[i,j] <= 0.05){
        xxx <- paste0(xxx, "}")
      }
    }
  }
  xxx <- paste0(xxx, " \\\\")
  cat(xxx, "\n")
}

##### 論文改訂用 (2026/01) #####

# --- 1. データの準備（インデックス比較による確実な抽出） ---

rn <- model.names

# --- 関数: マトリックスをLong形式に変換し、インデックスを付与する ---
make_long_data <- function(mat_val, mat_pval, type) {
  # 値データの変換
  df_val <- as.data.frame(mat_val)
  rownames(df_val) <- rn
  colnames(df_val) <- rn
  
  df_val <- df_val |>
    rownames_to_column("row") |>
    pivot_longer(-row, names_to = "col", values_to = "val")
  
  # p値データの変換
  df_pval <- as.data.frame(mat_pval)
  rownames(df_pval) <- rn
  colnames(df_pval) <- rn
  
  df_pval <- df_pval |>
    rownames_to_column("row") |>
    pivot_longer(-row, names_to = "col", values_to = "pval")
  
  # 結合
  df <- left_join(df_val, df_pval, by = c("row", "col"))
  
  # インデックス（番号）を追加してフィルタリング
  df <- df |>
    mutate(
      row_idx = match(row, rn),
      col_idx = match(col, rn)
    )
  
  # タイプに応じた抽出
  if (type == "lower") {
    # 下三角: 行番号 > 列番号
    df <- df |> filter(row_idx > col_idx)
  } else if (type == "upper") {
    # 上三角: 行番号 < 列番号
    df <- df |> filter(row_idx < col_idx)
  }
  
  return(df)
}

# --- スター判定関数 ---
get_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ ""
  )
}

# 色設定
col_good <- "#00BA38" # 緑
col_bad  <- "#F8766D" # オレンジ
val_alpha <- 0.6      # 透明度 (文字が見やすいように少し濃いめに調整)

# --- データの作成 ---
# 下三角（t統計量）
df_lower <- make_long_data(stat, pval, "lower")

# 上三角（勝率 %）
df_upper <- make_long_data(prop * 100, pval, "upper")

# 対角成分（背景用）
df_diag <- data.frame(row = rn, col = rn, val = 1)

# --- Factor順序の統一（プロットの配置決定） ---
all_levels <- rn

# Y軸は逆順（上から下）、X軸は正順
df_lower$row <- factor(df_lower$row, levels = rev(all_levels))
df_lower$col <- factor(df_lower$col, levels = all_levels)

df_upper$row <- factor(df_upper$row, levels = rev(all_levels))
df_upper$col <- factor(df_upper$col, levels = all_levels)

df_diag$row <- factor(df_diag$row, levels = rev(all_levels))
df_diag$col <- factor(df_diag$col, levels = all_levels)


# --- 2. プロット作成 ---

g <- ggplot() +
  
  # --- 1. 対角成分 ---
  geom_tile(data = df_diag, aes(x = col, y = row), 
            fill = "grey95", color = "white") +
  
  # --- 2. 下三角 (t値) ---
  geom_tile(data = df_lower, aes(x = col, y = row, fill = val), 
            color = "white", linewidth = 0.5, alpha = val_alpha) +
  geom_text(data = df_lower, 
            aes(x = col, y = row, label = paste0(sprintf("%.1f", val), get_stars(pval))),
            size = 5, fontface = "bold", color = "black") +
  
  # 下三角スケール
  scale_fill_gradient2(name = "t-stat\n(Lower)", 
                       low = col_good, mid = "white", high = col_bad, midpoint = 0) +
  
  # --- 新しいスケールへ ---
  new_scale_fill() + 
  
  # --- 3. 上三角 (勝率) ---
  geom_tile(data = df_upper, aes(x = col, y = row, fill = val), 
            color = "white", linewidth = 0.5, alpha = val_alpha) +
  geom_text(data = df_upper, 
            aes(x = col, y = row, label = paste0(sprintf("%.0f", val), get_stars(pval))),
            size = 5, fontface = "bold", color = "black") +
  
  # 上三角スケール
  scale_fill_gradient2(name = "Prop%\n(Upper)", 
                       low = col_good, mid = "white", high = col_bad, 
                       midpoint = 50, limits = c(0, 100)) +
  
  # --- 4. デザイン調整 ---
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    # X軸ラベル: 斜め45度、右揃え
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm"),
    legend.title = element_text(size = 11, face = "bold")
  )

print(g)
ggsave(paste0(path.figure, data.name, "_heatmap_gw_vol_", tolower(rv.name), ".png"),
       plot = g, width = 12, height = 10)
}

#### Run for all (data, volatility proxy) combinations ####
for (data.i in 1:data.n) {
  for (rv.i in 1:rv.n) {
    run_output_vol(data.i, rv.i)
  }
}
