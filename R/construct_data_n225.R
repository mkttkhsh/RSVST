# 必要なライブラリを読み込む
library(dplyr)
library(readxl)

# インデックスデータを読み込む
file_path <- "../Data/Before2022/NK225_19960104-20190930.xlsx"
data_index <- read_excel(file_path)

# yyyymmdd, rccを生成
data_rcc_1 <- data_index %>%
  mutate(
    yyyymmdd = as.numeric(gsub("-", "", substr(date, 1, 10))), # dateの最初の10文字を抽出し、ハイフンを削除して数値型に変換
    rcc = c(0, diff(log(close_price))) # 対数リターンを計算
  ) %>%
  mutate(
    yyyymmdd = ifelse(yyyymmdd > 20190939, date_value, yyyymmdd) # 日付を修正
  ) %>%
  select(yyyymmdd, rcc)

# インデックスデータを読み込む
file_path <- "../Data/nikkei/nikkei_20190919_20251219.csv"
data_index <- read.csv(file_path)

# yyyymmdd, rccを生成
data_rcc_2 <- data_index %>%
  mutate(
    yyyy = as.numeric(substr(date, 1, 4)),
    mm = as.numeric(gsub("/", "", substr(date, 6, 7))),
    dd = as.numeric(gsub("/", "", substr(date, 8, 10))),
    yyyymmdd = 10000*yyyy + 100*mm + dd,
    ymd = as.numeric(gsub("/", "", substr(date, 1, 10))), # dateの最初の10文字を抽出し、ハイフンを削除して数値型に変換
    rcc = c(0, diff(log(as.numeric(gsub(",", "", close))))) # 対数リターンを計算
  ) %>%
  select(yyyymmdd, rcc)

# リターンデータを統合する
data_rcc <- rbind(data_rcc_1, data_rcc_2[data_rcc_2$yyyymmdd >= 20191001, ]) %>%
  filter(yyyymmdd >= 19960311, yyyymmdd <= 20241231)

# RVデータを読み込む
file_path <- "../Data/Nikkei_RV_1996_202509.xlsx"
data_rv <- read_excel(file_path, sheet = "RV")
data_ftrk <- read_excel(file_path, sheet = "FTRK")
data_ftrk <- read_excel(file_path, sheet = "FTRK", skip = 1)
data_bpv <- read_excel(file_path, sheet = "BPV")
data_medrv <- read_excel(file_path, sheet = "MedRV")

# 必要なデータを選択
data_rv5_rkth2_bv_medrv <- data.frame(yyyymmdd = data_rv[[1]],
                                      rv5 = data_rv[[6]],
                                      rk_th2 = data_ftrk[[6]],
                                      bv = data_bpv[[6]],
                                      medrv = data_medrv[[2]]) %>%
  filter(yyyymmdd <= 20241231)

# リターンとRVのデータを統合する
data_rcc_rv5_rkth2_bv_medrv <- full_join(data_rcc, data_rv5_rkth2_bv_medrv, by = "yyyymmdd")

# マッチしなかった日を確認
anti_join(data_rcc, data_rv5_rkth2_bv_medrv, by = "yyyymmdd")

# 確認
head(data_rcc_rv5_rkth2_bv_medrv)

# CSVファイルとして保存
write.csv(data_rcc_rv5_rkth2_bv_medrv, "../Data/N225_rcc_rv5_rkth2_bv_medrv_199603_202412.csv", row.names = FALSE)

##### COVID期予測用サンプル #####

# 該当期間を抽出
data_rcc_rv5_rkth2_bv_medrv_2012_2020 <- data_rcc_rv5_rkth2_bv_medrv %>%
  filter(yyyymmdd >= 20120101, yyyymmdd <= 20201231)

# 確認
head(data_rcc_rv5_rkth2_bv_medrv_2012_2020)

# CSVファイルとして保存
write.csv(data_rcc_rv5_rkth2_bv_medrv_2012_2020,
          "../Data/N225_rcc_rv5_rkth2_bv_medrv_2012_2020.csv", row.names = FALSE)
