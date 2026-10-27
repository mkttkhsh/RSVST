# 必要なライブラリを読み込む
library(dplyr)

# データを読み込む
file_path <- "../Data/N225_rcc_rv5_rkth2_bv_medrv.csv"
data <- read.csv(file_path)

# yyyymmddが20090601以上20190930以下の行を抽出
filtered_data <- data %>%
  filter(yyyymmdd >= 20090601  & yyyymmdd <= 20190930)


# サンプルサイズの設定
nobs <- nrow(filtered_data)
nprd <- 590  # 予測サンプルサイズ
nest <- nobs - nprd

# RV推定量の数
nrv <- ncol(filtered_data) - 2

# 結果の格納用
c_hl <- matrix(0, nprd, nrv)
rv_hl <- matrix(0, nprd, nrv)

# 修正RVの計算
for (t in 1:nprd) {
  r <- data$rcc[(t + 1):(t + nest)]
  rv <- as.matrix(filtered_data[t:(t + nest - 1), 3:(2 + nrv)])  # 明示的に行列に変換
  m <- mean(r)
  
  # Hansen and Lunde (2005)の修正項と修正RVを計算
  c_hl[t,] <- sum((r - m)^2) / colSums(rv)
  rv_hl[t,] <- c_hl[t,] * rv[nest, ]
}

# c_hlとrv_hlに列名を設定
colnames(c_hl) <- paste0("c_hl_", colnames(rv))
colnames(rv_hl) <- paste0(colnames(rv), "_hl")

# filtered_data, c_hl, rv_hlを列方向に結合
merged_data <- cbind(filtered_data[(nest + 1):nobs, ], c_hl, rv_hl)

# 統合されたデータを確認
head(merged_data)

# 結果を保存
write.csv(merged_data, "../Data/N225_rcc_rv5_rkth2_bv_medrv_hl_roll.csv", row.names = FALSE)
