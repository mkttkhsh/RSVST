# 必要なライブラリを読み込む
library(dplyr)

# データを読み込む
file_path <- "../Data/oxfordmanrealizedvolatilityindices.csv"
data <- read.csv(file_path)

# `.DJI`シンボルのデータをフィルタリング
filtered_data <- data %>%
  filter(Symbol == ".DJI") %>%
  select(
    date = X, # 'Unnamed: 0' を 'date' にリネーム
    Symbol,
    rk_th2,
    medrv,
    bv,
    close_price,
    rv5
  )

# 日付列を確認
head(filtered_data$date)

# yyyymmdd, rccを生成
filtered_data <- filtered_data %>%
  mutate(
    yyyymmdd = as.numeric(gsub("-", "", substr(date, 1, 10))), # dateの最初の10文字を抽出し、ハイフンを削除して数値型に変換
    rcc = c(0, diff(log(close_price))) # 対数リターンを計算
  )

# 結果を表示
head(filtered_data)

# 指定した列の順序に並べ替え
ordered_data <- filtered_data %>%
  select(yyyymmdd, rcc, rv5, rk_th2, bv, medrv)

# CSVファイルとして保存
write.csv(ordered_data, "../Data/DJIA_rcc_rv5_rkth2_bv_medrv.csv", row.names = FALSE)

# 確認
head(ordered_data)
