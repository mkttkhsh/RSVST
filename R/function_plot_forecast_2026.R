library(patchwork)
library(dplyr)
library(ggplot2)
library(reshape2)

# set colors
ggColorHue <- function(n, l=65) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=100)[1:n]
}

plot_var_es <- function(f.var.es, out.name, size, i.legend, 
                        legend.xy = c(0.02, 0.98), # デフォルトを左上に変更
                        range.y){
  
  # melt data
  f.var.es.melt <- gather(f.var.es, key = "key", value = "value", c(-1), factor_key = TRUE)
  
  # time-series plot
  g <- ggplot(f.var.es.melt, aes(x = date, y = value, linetype = key, color = key, shape = key, size = key))
  g <- g + theme_bw(base_size = size)
  g <- g + geom_line()
  g <- g + geom_point(aes(stroke = 2))
  g <- g + scale_linetype_manual(values = c("solid", "solid", "solid", "blank"))
  g <- g + scale_shape_manual(values = c(32, 32, 32, 4))
  g <- g + scale_size_manual(values = c(0.3, 0.6, 0.8, 3))
  g <- g + xlab("") + ylab("")
  
  if(i.legend){
    g <- g + theme(
      legend.position = legend.xy,            # 指定座標（デフォルトは左上）
      legend.justification = c("left", "top"),# 左上を基点にする
      legend.title = element_blank(),
      legend.background = element_rect(fill = scales::alpha("white", 0.6), color = NA), # 半透明・枠なし
      legend.key = element_rect(fill = NA),   # キー背景透明
      legend.text = element_text(size = size * 0.7), # 文字サイズ調整
      legend.margin = margin(2, 2, 2, 2)
    )
    # 凡例を1列に揃える（お好みで）
    # g <- g + guides(color = guide_legend(ncol = 1)) 
  }
  else {
    g <- g + theme(legend.position = "none")
  }
  
  g <- g + theme(axis.text = element_text(size = size))
  g <- g + scale_x_date(date_labels = "%Y-%m")
  g <- g + ylim(range.y)
  
  ggsave(file = out.name, plot = g, dpi = 100, width = 16, height = 10)
}

# モデル名からグループを抽出するヘルパー関数
get_model_group <- function(model_names) {
  groups <- case_when(
    grepl("^SV", model_names) ~ "SV Family",
    grepl("^RSV", model_names) ~ "RSV Family",
    grepl("^EGARCH", model_names) ~ "EGARCH Family",
    grepl("^REGARCH", model_names) ~ "REGARCH Family",
    TRUE ~ "Other" # Proxy (RV5 etc) の場合
  )
  # Proxy自体はすべてのFacetに薄く表示するか、あるいは別扱いにする必要がありますが、
  # ここでは予測モデルの比較に主眼を置くため、予測モデルのグルーピングを優先します。
  return(groups)
}

##### Facet対応版のプロット関数 #####

plot_vol_faceted <- function(f.vol, out.name, size, cols, 
                             legend.pos.type = "topright", # 選べるように追加
                             i.legend = TRUE, i.title = 0, s.title = ""){
  
  f.vol.melt <- melt(f.vol, id.vars = "date")
  proxy_name <- names(f.vol)[2] 
  
  f.vol.melt <- f.vol.melt %>%
    mutate(group = case_when(
      variable == proxy_name ~ "Proxy",
      grepl("^SV", variable) ~ "SV",
      grepl("^RSV", variable) ~ "RSV",
      grepl("^EGARCH", variable) ~ "EGARCH",
      grepl("^REGARCH", variable) ~ "REGARCH",
      TRUE ~ "Other"
    ))
  
  y_range <- range(f.vol.melt$value, na.rm = TRUE)
  
  # --- 線種と形状の定義 ---
  # 各ファミリー内で最大7モデル程度と想定し、区別しやすいスタイルを用意
  custom_ltys <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "11")
  custom_shapes <- c(16, 17, 15, 18, 1, 2, 0) # 塗りつぶしと白抜き
  
  model_cols <- cols[-1]
  names(model_cols) <- model.names
  
  plot_list <- list()
  family_levels <- c("SV", "RSV", "EGARCH", "REGARCH")
  
  # 凡例位置の座標変換
  leg_settings <- switch(legend.pos.type,
                         "topleft"     = list(pos = c(0.02, 0.98), just = c("left", "top")),
                         "topright"    = list(pos = c(0.98, 0.98), just = c("right", "top")),
                         "bottomleft"  = list(pos = c(0.02, 0.02), just = c("left", "bottom")),
                         "bottomright" = list(pos = c(0.98, 0.02), just = c("right", "bottom"))
  )
  
  for(i in seq_along(family_levels)){
    fam <- family_levels[i]
    df_fam <- f.vol.melt %>% filter(group == fam | group == "Proxy")
    df_fam$group <- fam 
    
    # ファミリー内のモデル一覧を取得して線種等を割り当て
    models_in_fam <- setdiff(unique(df_fam$variable), proxy_name)
    
    p <- ggplot(df_fam, aes(x = date, y = value)) + theme_bw(base_size = size)
    
    # (1) 背景（Proxy）
    p <- p + geom_line(data = subset(df_fam, variable == proxy_name),
                       color = "gray70", linetype = "solid", size = 0.5)
    
    # (2) モデル（線種を分ける）
    # aes(linetype = variable) を追加
    p <- p + geom_line(data = subset(df_fam, variable != proxy_name),
                       aes(color = variable, linetype = variable), size = 0.8)
    
    # (3) モデル（形状を分ける：データが多いので間引いて表示）
    # 30日ごとに1点プロットする例
    df_points <- subset(df_fam, variable != proxy_name) %>%
      group_by(variable) %>%
      filter(row_number() %% 30 == 0) 
    
    p <- p + geom_point(data = df_points, 
                        aes(color = variable, shape = variable), size = 2)
    
    # (4) スケール設定
    p <- p + facet_wrap(~ group)
    p <- p + scale_color_manual(values = model_cols)
    p <- p + scale_linetype_manual(values = setNames(rep(custom_ltys, length.out=length(models_in_fam)), models_in_fam))
    p <- p + scale_shape_manual(values = setNames(rep(custom_shapes, length.out=length(models_in_fam)), models_in_fam))
    
    p <- p + coord_cartesian(ylim = y_range)
    p <- p + xlab("") + ylab("") + scale_x_date(date_labels = "%Y-%m")
    
    # (5) 凡例
    if(i.legend){
      p <- p + theme(
        legend.position = leg_settings$pos,
        legend.justification = leg_settings$just,
        legend.title = element_blank(),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), color = "gray80"),
        legend.text = element_text(size = size * 0.5),
        legend.key.width = unit(2, "line") # 線種を見やすくするために幅を広げる
      )
      p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
    } else {
      p <- p + theme(legend.position = "none")
    }
    
    # 軸の制御
    if(i %% 2 == 0) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if(i <= 2) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    else p <- p + theme(axis.text.x = element_text(size = size * 0.7, angle = 45, hjust = 1))
    
    p <- p + theme(strip.text = element_text(size = size, face = "bold"))
    plot_list[[i]] <- p
  }
  
  g_combined <- wrap_plots(plot_list, ncol = 2)
  ggsave(file = out.name, plot = g_combined, dpi = 300, width = 16, height = 10)
}

# --- 凡例を外側（右側）にまとめる修正版ロジック ---

plot_vol_faceted_outer_legend <- function(f.vol, out.name, size, cols, i.legend = TRUE){
  
  f.vol.melt <- melt(f.vol, id.vars = "date")
  proxy_name <- names(f.vol)[2] 
  
  f.vol.melt <- f.vol.melt %>%
    mutate(group = case_when(
      variable == proxy_name ~ "Proxy",
      grepl("^SV", variable) ~ "SV",
      grepl("^RSV", variable) ~ "RSV",
      grepl("^EGARCH", variable) ~ "EGARCH",
      grepl("^REGARCH", variable) ~ "REGARCH",
      TRUE ~ "Other"
    ))
  
  y_range <- range(f.vol.melt$value, na.rm = TRUE)
  
  # --- 線種と形状の定義 ---
  # 各ファミリー内で最大7モデル程度と想定し、区別しやすいスタイルを用意
  custom_ltys <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "11")
  custom_shapes <- c(16, 17, 15, 18, 1, 2, 0) # 塗りつぶしと白抜き
  
  model_cols <- cols[-1]
  names(model_cols) <- model.names
  
  plot_list <- list()
  family_levels <- c("SV", "RSV", "EGARCH", "REGARCH")
  
  for(i in seq_along(family_levels)){
    fam <- family_levels[i]
    df_fam <- f.vol.melt %>% filter(group == fam | group == "Proxy")
    df_fam$group <- fam 
    
    # ファミリー内のモデル一覧を取得して線種等を割り当て
    models_in_fam <- setdiff(unique(df_fam$variable), proxy_name)
    
    p <- ggplot(df_fam, aes(x = date, y = value)) + theme_bw(base_size = size)
    
    # (1) 背景（Proxy）
    p <- p + geom_line(data = subset(df_fam, variable == proxy_name),
                       color = "gray70", linetype = "solid", size = 0.5)
    
    # (2) モデル（線種を分ける）
    # aes(linetype = variable) を追加
    p <- p + geom_line(data = subset(df_fam, variable != proxy_name),
                       aes(color = variable, linetype = variable), size = 0.8)
    
    # (3) モデル（形状を分ける：データが多いので間引いて表示）
    # 30日ごとに1点プロットする例
    df_points <- subset(df_fam, variable != proxy_name) %>%
      group_by(variable) %>%
      filter(row_number() %% 30 == 0) 
    
    p <- p + geom_point(data = df_points, 
                        aes(color = variable, shape = variable), size = 2)
    
    # (4) スケール設定
    p <- p + facet_wrap(~ group)
    p <- p + scale_color_manual(values = model_cols)
    p <- p + scale_linetype_manual(values = setNames(rep(custom_ltys, length.out=length(models_in_fam)), models_in_fam))
    p <- p + scale_shape_manual(values = setNames(rep(custom_shapes, length.out=length(models_in_fam)), models_in_fam))
    
    p <- p + coord_cartesian(ylim = y_range)
    p <- p + xlab("") + ylab("") + scale_x_date(date_labels = "%Y-%m")
    
    # (5) 凡例
    # ★ ここで各プロットの凡例は非表示にする
    p <- p + theme(legend.position = "none") 
    
    # ただし、凡例を抽出するために一つだけ情報を保持させる設定にする
    # (実際には後で guides = "collect" でまとめます)
    p <- p + theme(legend.title = element_blank(),
                   legend.text = element_text(size = size * 0.6))
    
    plot_list[[i]] <- p
  }
  
  # --- 3. 結合して凡例を右側に配置 ---
  # guides = "collect" で全パネルの凡例を一つにまとめます
  g_combined <- wrap_plots(plot_list, ncol = 2) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "right") # ここで右側に指定
  
  ggsave(file = out.name, plot = g_combined, dpi = 300, width = 16, height = 10)
}

plot_cld_faceted <- function(f.cld, out.name, size, cols, 
                             legend.pos.type = "topleft", # 4角指定用に追加
                             i.legend = TRUE, s.title = "", i.title = 0){
  
  # 1. データ準備
  f.cld.melt <- melt(f.cld, id.vars = "date")
  
  f.cld.melt <- f.cld.melt %>%
    mutate(group = case_when(
      grepl("^SV", variable) ~ "SV",
      grepl("^RSV", variable) ~ "RSV",
      grepl("^EGARCH", variable) ~ "EGARCH",
      grepl("^REGARCH", variable) ~ "REGARCH",
      TRUE ~ "Other"
    ))
  
  y_range <- range(f.cld.melt$value, na.rm = TRUE)
  
  # --- 線種と形状の定義（ボラティリティプロットと統一） ---
  custom_ltys <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "11")
  custom_shapes <- c(16, 17, 15, 18, 1, 2, 0)
  
  model_cols <- cols[-1]
  names(model_cols) <- model.names 
  
  # 凡例位置の座標変換
  leg_settings <- switch(legend.pos.type,
                         "topleft"     = list(pos = c(0.02, 0.98), just = c("left", "top")),
                         "topright"    = list(pos = c(0.98, 0.98), just = c("right", "top")),
                         "bottomleft"  = list(pos = c(0.02, 0.02), just = c("left", "bottom")),
                         "bottomright" = list(pos = c(0.98, 0.02), just = c("right", "bottom"))
  )
  
  plot_list <- list()
  family_levels <- c("SV", "RSV", "EGARCH", "REGARCH")
  
  # 2. 4つのグラフを作成
  for(i in seq_along(family_levels)){
    
    fam <- family_levels[i]
    df_fam <- f.cld.melt %>% filter(group == fam)
    df_fam$facet_label <- fam
    
    # ファミリー内のモデル一覧を取得
    models_in_fam <- unique(df_fam$variable)
    
    p <- ggplot(df_fam, aes(x = date, y = value)) + theme_bw(base_size = size)
    
    # (1) 基準線 (y=0)
    p <- p + geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5)
    
    # (2) モデルの累積差分（線種を分ける）
    if(nrow(df_fam) > 0){
      p <- p + geom_line(aes(color = variable, linetype = variable), size = 0.8)
      
      # (3) ポイントの追加（視認性向上のため間引く：40日ごと）
      df_points <- df_fam %>%
        group_by(variable) %>%
        filter(row_number() %% 40 == 0)
      
      p <- p + geom_point(data = df_points, 
                          aes(color = variable, shape = variable), size = 2)
    }
    
    # (4) デザイン設定
    p <- p + facet_wrap(~ facet_label)
    p <- p + scale_color_manual(values = model_cols)
    p <- p + scale_linetype_manual(values = setNames(rep(custom_ltys, length.out=length(models_in_fam)), models_in_fam))
    p <- p + scale_shape_manual(values = setNames(rep(custom_shapes, length.out=length(models_in_fam)), models_in_fam))
    
    p <- p + coord_cartesian(ylim = y_range) 
    p <- p + xlab("") + ylab("")
    p <- p + scale_x_date(date_labels = "%Y-%m")
    
    # 凡例の設定
    if(i.legend){
      p <- p + theme(
        legend.position = leg_settings$pos,
        legend.justification = leg_settings$just,
        legend.title = element_blank(),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), color = "gray80"),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = size * 0.55),
        legend.key.width = unit(1.8, "line")
      )
      p <- p + guides(color = guide_legend(ncol = 1))
    } else {
      p <- p + theme(legend.position = "none")
    }
    
    # 軸の表示制御
    if(i %% 2 == 0) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if(i <= 2) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    else p <- p + theme(axis.text.x = element_text(size = size * 0.7, angle = 45, hjust = 1))
    
    p <- p + theme(axis.text.y = element_text(size = size * 0.7))
    p <- p + theme(strip.text = element_text(size = size, face = "bold"))
    
    plot_list[[i]] <- p
  }
  
  # 3. 結合して出力
  g_combined <- wrap_plots(plot_list, ncol = 2)
  if(i.title) g_combined <- g_combined + plot_annotation(title = s.title, theme = theme(plot.title = element_text(size = size)))
  
  ggsave(file = out.name, plot = g_combined, dpi = 300, width = 16, height = 10)
}

plot_cld_faceted_outer_legend <- function(f.cld, out.name, size, cols, i.legend = TRUE, 
                             s.title = "", i.title = 0){
  
  # 1. データ準備
  f.cld.melt <- melt(f.cld, id.vars = "date")
  
  # グループ分け（モデルの接頭辞で判定）
  f.cld.melt <- f.cld.melt %>%
    mutate(group = case_when(
      grepl("^SV", variable) ~ "SV",
      grepl("^RSV", variable) ~ "RSV",
      grepl("^EGARCH", variable) ~ "EGARCH",
      grepl("^REGARCH", variable) ~ "REGARCH",
      TRUE ~ "Other"
    ))
  
  # y軸の範囲を固定
  y_range <- range(f.cld.melt$value, na.rm = TRUE)

  # 線種と形状の定義（モデル判別用）
  custom_ltys <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "11")
  custom_shapes <- c(16, 17, 15, 18, 1, 2, 0)
  
  # 色のマッピング（model.namesはグローバル変数を参照）
  model_cols <- cols[-1]
  names(model_cols) <- model.names 
  
  plot_list <- list()
  family_levels <- c("SV", "RSV", "EGARCH", "REGARCH")
  
  # 2. 4つのファミリーごとにプロット作成
  for(i in seq_along(family_levels)){
    
    fam <- family_levels[i]
    df_fam <- f.cld.melt %>% filter(group == fam)
    df_fam$facet_label <- fam # ヘッダー用
    
    models_in_fam <- unique(df_fam$variable)
    
    p <- ggplot(df_fam, aes(x = date, y = value)) + theme_bw(base_size = size)
    
    # (1) 基準線 (y=0)
    p <- p + geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.3)
    
    # (2) モデルの累積差分（線種・色を適用）
    if(nrow(df_fam) > 0){
      p <- p + geom_line(aes(color = variable, linetype = variable), size = 0.8)
      
      # (3) ポイントを追加（40日ごとに間引き）
      df_points <- df_fam %>%
        group_by(variable) %>%
        filter(row_number() %% 40 == 0)
      
      p <- p + geom_point(data = df_points, 
                          aes(color = variable, shape = variable), size = 2)
    }
    
    # (4) ファセットとスケールの設定
    p <- p + facet_wrap(~ facet_label)
    p <- p + scale_color_manual(values = model_cols)
    p <- p + scale_linetype_manual(values = setNames(rep(custom_ltys, length.out=length(models_in_fam)), models_in_fam))
    p <- p + scale_shape_manual(values = setNames(rep(custom_shapes, length.out=length(models_in_fam)), models_in_fam))
    
    # y軸を固定
    p <- p + coord_cartesian(ylim = y_range) 
    p <- p + xlab("") + ylab("")
    p <- p + scale_x_date(date_labels = "%Y-%m")
    
    # 個別プロットのテーマ設定
    p <- p + theme(
      legend.title = element_blank(),
      legend.text = element_text(size = size * 0.6),
      legend.key.width = unit(2, "line"),
      strip.text = element_text(size = size, face = "bold"),
      axis.text.y = element_text(size = size * 0.7)
    )
    
    # 軸表示の整理（内側の軸を消す）
    if(i %% 2 == 0) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if(i <= 2) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    else p <- p + theme(axis.text.x = element_text(size = size * 0.7, angle = 45, hjust = 1))
    
    plot_list[[i]] <- p
  }
  
  # 3. 結合して凡例を外側（右）に配置
  g_combined <- wrap_plots(plot_list, ncol = 2) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = if(i.legend) "right" else "none")
  
  if(i.title){
    g_combined <- g_combined + plot_annotation(title = s.title, theme = theme(plot.title = element_text(size = size)))
  }
  
  ggsave(file = out.name, plot = g_combined, dpi = 300, width = 18, height = 10, bg = "white")
}

plot_fz0_cld_faceted <- function(f.cld, out.name, size, cols, 
                                 legend.pos.type = "topleft", # 引数を追加
                                 i.legend = TRUE, s.title = "", i.title = 0){
  
  # 1. データ準備
  f.cld.melt <- melt(f.cld, id.vars = "date")
  
  f.cld.melt <- f.cld.melt %>%
    mutate(group = case_when(
      grepl("^SV", variable) ~ "SV",
      grepl("^RSV", variable) ~ "RSV",
      grepl("^EGARCH", variable) ~ "EGARCH",
      grepl("^REGARCH", variable) ~ "REGARCH",
      TRUE ~ "Other"
    ))
  
  y_range <- range(f.cld.melt$value, na.rm = TRUE)
  
  # 線種と形状の定義
  custom_ltys <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "11")
  custom_shapes <- c(16, 17, 15, 18, 1, 2, 0)
  
  model_cols <- cols[-1]
  names(model_cols) <- model.names 
  
  # 基本となる凡例位置の座標変換
  leg_settings <- switch(legend.pos.type,
                         "topleft"     = list(pos = c(0.02, 0.98), just = c("left", "top")),
                         "topright"    = list(pos = c(0.98, 0.98), just = c("right", "top")),
                         "bottomleft"  = list(pos = c(0.02, 0.02), just = c("left", "bottom")),
                         "bottomright" = list(pos = c(0.98, 0.02), just = c("right", "bottom"))
  )
  
  plot_list <- list()
  family_levels <- c("SV", "RSV", "EGARCH", "REGARCH")
  
  # 2. ループでグラフ作成
  for(i in seq_along(family_levels)){
    
    fam <- family_levels[i]
    df_fam <- f.cld.melt %>% filter(group == fam)
    df_fam$facet_label <- fam
    
    models_in_fam <- unique(df_fam$variable)
    
    # ★ 特定のファミリーだけ位置を逆にしたい場合の個別ロジック（任意）
    # グラフの推移を見て重なる場合は、ここで fam ごとに switch しても良いです。
    # 現在は引数 legend.pos.type の設定に従います。

    p <- ggplot(df_fam, aes(x = date, y = value)) + theme_bw(base_size = size)
    p <- p + geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5)
    
    if(nrow(df_fam) > 0){
      # 線種と色を適用
      p <- p + geom_line(aes(color = variable, linetype = variable), size = 0.8)
      
      # ポイントを追加（40日ごとに間引き）
      df_points <- df_fam %>%
        group_by(variable) %>%
        filter(row_number() %% 40 == 0)
      
      p <- p + geom_point(data = df_points, 
                          aes(color = variable, shape = variable), size = 2)
    }
    
    p <- p + facet_wrap(~ facet_label)
    p <- p + scale_color_manual(values = model_cols)
    p <- p + scale_linetype_manual(values = setNames(rep(custom_ltys, length.out=length(models_in_fam)), models_in_fam))
    p <- p + scale_shape_manual(values = setNames(rep(custom_shapes, length.out=length(models_in_fam)), models_in_fam))
    
    # p <- p + coord_cartesian(ylim = y_range) 
    p <- p + coord_cartesian(ylim = range(y_range, -max(abs(y_range)), max(abs(y_range))))    
    p <- p + xlab("") + ylab("")
    p <- p + scale_x_date(date_labels = "%Y-%m")
    
    if(i.legend){
      p <- p + theme(
        legend.position = leg_settings$pos,
        legend.justification = leg_settings$just,
        legend.title = element_blank(),
        legend.background = element_rect(fill = scales::alpha("white", 0.6), color = "gray80"),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = size * 0.55),
        legend.margin = margin(2, 2, 2, 2),
        legend.key.width = unit(1.8, "line")
      )
      p <- p + guides(color = guide_legend(ncol = 1))
    } else {
      p <- p + theme(legend.position = "none")
    }
    
    # 軸・ラベルの制御
    if(i %% 2 == 0) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if(i <= 2) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    else p <- p + theme(axis.text.x = element_text(size = size * 0.7, angle = 45, hjust = 1))
    
    p <- p + theme(axis.text.y = element_text(size = size * 0.7))
    p <- p + theme(strip.text = element_text(size = size, face = "bold"))
    
    plot_list[[i]] <- p
  }
  
  g_combined <- wrap_plots(plot_list, ncol = 2)
  if(i.title) g_combined <- g_combined + plot_annotation(title = s.title, theme = theme(plot.title = element_text(size = size)))
  
  ggsave(file = out.name, plot = g_combined, dpi = 300, width = 16, height = 10)
}

plot_fz0_cld_faceted_outer_legend <- function(f.cld, out.name, size, cols, i.legend = TRUE, 
                                 s.title = "", i.title = 0){
  
  # 1. データ準備
  f.cld.melt <- melt(f.cld, id.vars = "date")
  
  f.cld.melt <- f.cld.melt %>%
    mutate(group = case_when(
      grepl("^SV", variable) ~ "SV",
      grepl("^RSV", variable) ~ "RSV",
      grepl("^EGARCH", variable) ~ "EGARCH",
      grepl("^REGARCH", variable) ~ "REGARCH",
      TRUE ~ "Other"
    ))
  
  # y軸の範囲を固定
  y_range <- range(f.cld.melt$value, na.rm = TRUE)

  custom_ltys <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "11")
  custom_shapes <- c(16, 17, 15, 18, 1, 2, 0)
  
  model_cols <- cols[-1]
  names(model_cols) <- model.names 
  
  plot_list <- list()
  family_levels <- c("SV", "RSV", "EGARCH", "REGARCH")
  
  # 2. ループ
  for(i in seq_along(family_levels)){
    
    fam <- family_levels[i]
    df_fam <- f.cld.melt %>% filter(group == fam)
    df_fam$facet_label <- fam
    
    models_in_fam <- unique(df_fam$variable)
    
    p <- ggplot(df_fam, aes(x = date, y = value)) + theme_bw(base_size = size)
    p <- p + geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.3)
    
    if(nrow(df_fam) > 0){
      p <- p + geom_line(aes(color = variable, linetype = variable), size = 0.8)
      
      df_points <- df_fam %>%
        group_by(variable) %>%
        filter(row_number() %% 40 == 0)
      
      p <- p + geom_point(data = df_points, 
                          aes(color = variable, shape = variable), size = 2)
    }
    
    p <- p + facet_wrap(~ facet_label)
    p <- p + scale_color_manual(values = model_cols)
    p <- p + scale_linetype_manual(values = setNames(rep(custom_ltys, length.out=length(models_in_fam)), models_in_fam))
    p <- p + scale_shape_manual(values = setNames(rep(custom_shapes, length.out=length(models_in_fam)), models_in_fam))
    
    p <- p + coord_cartesian(ylim = y_range) 
    p <- p + xlab("") + ylab("")
    p <- p + scale_x_date(date_labels = "%Y-%m")
    
    p <- p + theme(
      legend.title = element_blank(),
      legend.text = element_text(size = size * 0.6),
      legend.key.width = unit(2, "line"),
      strip.text = element_text(size = size, face = "bold"),
      axis.text.y = element_text(size = size * 0.7)
    )
    
    if(i %% 2 == 0) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if(i <= 2) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    else p <- p + theme(axis.text.x = element_text(size = size * 0.7, angle = 45, hjust = 1))
    
    plot_list[[i]] <- p
  }
  
  # 3. 結合
  g_combined <- wrap_plots(plot_list, ncol = 2) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = if(i.legend) "right" else "none")
  
  if(i.title) g_combined <- g_combined + plot_annotation(title = s.title, theme = theme(plot.title = element_text(size = size)))
  
  ggsave(file = out.name, plot = g_combined, dpi = 300, width = 18, height = 10, bg = "white")
}
