# set colors
ggColorHue <- function(n, l=65) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=100)[1:n]
}

plot_vol <- function(f.vol, out.name, size, cols, ltys, i.legend, legend.xy = c(0,0), i.title = 0, s.title = ""){
  # melt data
  f.vol.melt <- melt(f.vol, id.vars = "date")
  
  # times-series plot
  g <- ggplot(f.vol.melt, aes(x = date, y = value, color = variable, linetype = variable)) + theme_bw(base_size = size) + geom_line()
  g <- g + theme(axis.text=element_text(size = size)) + theme(legend.text=element_text(size = size))
  # g <- g + theme(legend.position=c(0.85, 0.8)) + theme(legend.title = element_blank()) + xlab("") + ylab("")
  g <- g + xlab("") + ylab("")
  g <- g + scale_color_manual(values = cols) + scale_linetype_manual(values = ltys)
  if(i.legend){
    g <- g + theme(legend.title = element_blank())
    # if(length(legend.xy) == 1){
    #   
    #   g <- g + geom_text_repel(
    #     data = subset(f.vol.melt, conc == max(conc)),
    #     aes(label = names(f.vol.melt)[3:length]),
    #     nudge_x = 100,
    #     segment.alpha = 0.5,
    #     size = 3
    #   ) +
    #     lims(x = c(min(CO2$conc), max(CO2$conc)*1.1))
    # }
    if(length(legend.xy) == 1){
      g <- g + theme(legend.position = legend.xy) + theme(legend.title = element_blank())
      g <- g + guides(color = guide_legend(nrow = 3))
      # g <- g + theme(legend.background = element_rect(fill = "white", colour = "black"))
    }
    if(length(legend.xy) == 2){
      g <- g + theme(legend.position = legend.xy) + theme(legend.title = element_blank())
      g <- g + theme(legend.background = element_rect(fill = "white", colour = "black"))
    }
  }
  else{
    g <- g + theme(legend.position = "none")
  }
  if(i.title){
    g <- g + ggtitle(s.title)
  }
  g <- g + scale_x_date(date_labels = "%Y-%m")
#  g <- g + guides(fill=guide_legend(nrow = 1))
  # plot(g)
  ggsave(file = out.name, plot = g, dpi = 100, width = 16, height = 9)
}

plot_var_es <-   function(f.var.es, out.name, size, i.legend, legend.xy = c(0,0), range.y){
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
    g <- g + theme(legend.position = legend.xy) + theme(legend.title = element_blank())
    g <- g + theme(legend.background = element_rect(fill = "white", colour = "black"))
  }
  else
    g <- g + theme(legend.position = "none")
  g <- g + theme(axis.text=element_text(size = size))
  g <- g + scale_x_date(date_labels = "%Y-%m")
  g <- g + ylim(range.y)
  ggsave(file = out.name, plot = g, dpi = 100, width = 16, height = 9)
}


plot_var_es_cld <- function(f.vol, out.name, size, cols, ltys, i.legend, legend.xy = c(0,0), i.title = 0, s.title = ""){
  # melt data
  f.vol.melt <- melt(f.vol, id.vars = "date")
  
  # times-series plot
  g <- ggplot(f.vol.melt, aes(x = date, y = value, color = variable, linetype = variable)) + theme_bw(base_size = size) + geom_line()
  g <- g + theme(axis.text=element_text(size = size)) + theme(legend.text=element_text(size = size))
  # g <- g + theme(legend.position=c(0.85, 0.8)) + theme(legend.title = element_blank()) + xlab("") + ylab("")
  g <- g + xlab("") + ylab("")
  g <- g + scale_color_manual(values = cols) + scale_linetype_manual(values = ltys)
  if(i.legend){
    g <- g + theme(legend.title = element_blank())
    # if(length(legend.xy) == 1){
    #   
    #   g <- g + geom_text_repel(
    #     data = subset(f.vol.melt, conc == max(conc)),
    #     aes(label = names(f.vol.melt)[3:length]),
    #     nudge_x = 100,
    #     segment.alpha = 0.5,
    #     size = 3
    #   ) +
    #     lims(x = c(min(CO2$conc), max(CO2$conc)*1.1))
    # }
    if(length(legend.xy) == 1){
      g <- g + theme(legend.position = legend.xy) + theme(legend.title = element_blank())
      g <- g + guides(color = guide_legend(nrow = 3))
      # g <- g + theme(legend.background = element_rect(fill = "white", colour = "black"))
    }
    if(length(legend.xy) == 2){
      g <- g + theme(legend.position = legend.xy) + theme(legend.title = element_blank())
      g <- g + theme(legend.background = element_rect(fill = "white", colour = "black"))
      g <- g + guides(color = guide_legend(nrow = 3))
    }
  }
  else{
    g <- g + theme(legend.position = "none")
  }
  if(i.title){
    g <- g + ggtitle(s.title)
  }
  g <- g + scale_x_date(date_labels = "%Y-%m")
  #  g <- g + guides(fill=guide_legend(nrow = 1))
  # plot(g)
  ggsave(file = out.name, plot = g, dpi = 100, width = 16, height = 9)
}

# ---- 既存のライブラリ読み込みの後に追加 ----

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

# Facet対応版のプロット関数
plot_vol_faceted <- function(f.vol, out.name, size, cols, ltys, i.legend = TRUE, 
                             legend.xy = "bottom", i.title = 0, s.title = ""){
  
  # melt data
  f.vol.melt <- melt(f.vol, id.vars = "date")
  
  # ---- 追加: グルーピング処理 ----
  # 変数名に基づいてグループ列を作成
  f.vol.melt <- f.vol.melt %>%
    mutate(group = get_model_group(variable))
  
  # RV (Proxy) と モデル予測値を区別したい場合
  # 通常、Proxy（実測値）は黒線などで全パネルに表示したいことが多いですが、
  # ここでは単純にモデル群を分割表示する形式にします。
  
  # times-series plot
  g <- ggplot(f.vol.melt, aes(x = date, y = value, color = variable, linetype = variable)) + 
    theme_bw(base_size = size) + 
    geom_line(size = 0.8) # 線を少し太くして見やすく
  
  # ---- 追加: Facet (サブプロット) ----
  # グループごとに分割。scales="free_y"にするとY軸が個別に調整されますが、
  # 比較のためには "fixed" が望ましいことが多いです。状況に合わせて変更してください。
  g <- g + facet_wrap(~ group, ncol = 2, scales = "fixed")
  
  g <- g + theme(axis.text = element_text(size = size * 0.8), # 軸ラベル調整
                 strip.text = element_text(size = size))      # Facetタイトルのサイズ
  
  g <- g + xlab("") + ylab("")
  g <- g + scale_color_manual(values = cols) + scale_linetype_manual(values = ltys)
  
  # 凡例の調整
  if(i.legend){
    g <- g + theme(legend.title = element_blank())
    g <- g + theme(legend.position = legend.xy) 
    g <- g + guides(color = guide_legend(nrow = 4)) # 凡例の行数を調整
  } else {
    g <- g + theme(legend.position = "none")
  }
  
  if(i.title){
    g <- g + ggtitle(s.title)
  }
  
  g <- g + scale_x_date(date_labels = "%Y-%m")
  
  # 保存 (サイズは大きめに確保)
  ggsave(file = out.name, plot = g, dpi = 300, width = 18, height = 12)
}