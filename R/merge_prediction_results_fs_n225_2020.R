rm(list = ls())

# --- 設定項目 ---
models <- c("SV_FS_SN", "SV_FS_ST", "RSV_FS_SN", "RSV_FS_ST")
n_pred <- 242  # 予測期間数（実際の数に合わせて変更してください）
probs <- c(0.01, 0.05) # 1%と5%
input_base_dir <- "../Ox/N225_prediction_2020"
output_dir <- "../Ox/N225_prediction_2020"

# 出力フォルダが存在しない場合は作成
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- メイン処理 ---
for (model in models) {
  cat(paste0("\nProcessing model: ", model, "\n"))
  
  # 結果を格納するリスト（後で行列に変換）
  var_results <- matrix(NA, nrow = n_pred, ncol = length(probs))
  es_results <- matrix(NA, nrow = n_pred, ncol = length(probs))
  volf_results <- matrix(NA, nrow = n_pred, ncol = 4) # Mean and Median for h and exp(h)
  
  for (x in 0:(n_pred-1)) {
    # ファイルパスの構築
    # ※拡張子は.csvと仮定しています。適宜変更してください。
    h_file <- file.path(input_base_dir, model, paste0("h_pred_s_", x, ".csv"))
    y_file <- file.path(input_base_dir, model, paste0("y_pred_s_", x, ".csv"))
    
    # ファイルが存在するか確認して読み込み
    if (file.exists(h_file) && file.exists(y_file)) {
      # データの読み込み（1列のデータを想定）
      h_samples <- read.csv(h_file, header = TRUE)[, 2]
      y_samples <- read.csv(y_file, header = TRUE)[, 2]
      
      # 1. VaR の計算 (Quantile)
      vvar <- quantile(y_samples, probs = probs)
      var_results[x + 1, ] <- vvar
      
      # 2. ES (Expected Shortfall) の計算
      # 各確率(1%, 5%)に対して、VaR以下のサンプルの平均をとる
      for (i in seq_along(probs)) {
        es_samples <- y_samples[y_samples < vvar[i]]
        es_results[x + 1, i] <- mean(es_samples)
      }
      
      # 3. ボラティリティ予測 (Mean and Median)
      # Oxのコード「meanc(vsamph1) ~ quantilec(vsamph1)」に準拠
      volf_results[x + 1, 1] <- mean(h_samples)
      volf_results[x + 1, 2] <- median(h_samples)
      volf_results[x + 1, 3] <- mean(exp(h_samples))
      volf_results[x + 1, 4] <- median(exp(h_samples))
      
    } else {
      warning(paste("File not found for x =", x, "in model", model))
    }
  }
  
  # --- モデルごとの結果保存 ---
  # 列名の設定
  colnames(var_results) <- paste0("VaR_", probs * 100, "p")
  colnames(es_results) <- paste0("ES_", probs * 100, "p")
  colnames(volf_results) <- c("Mean", "Median", "MeanExp", "MedianExp")
  
  # CSV書き出し
  write.csv(var_results, file.path(output_dir, paste0("var_", tolower(model), ".csv")), row.names = FALSE)
  write.csv(es_results, file.path(output_dir, paste0("es_", tolower(model), ".csv")), row.names = FALSE)
  write.csv(volf_results, file.path(output_dir, paste0("volf_", tolower(model), ".csv")), row.names = FALSE)
  
  cat(paste("Successfully saved results for", model, "\n"))
}

cat("\nAll processes completed.\n")