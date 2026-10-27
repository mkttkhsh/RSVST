rm(list = ls())

# packages
library(readr)
library(readxl)

# define sets of forecasts
f.sets <- c("volf", "var", "es")
f.n <- length(f.sets)

# select data
data.name <- "N225"

# set sets of true returns and RVs
true.sets <- c("N225_rcc_rv5_rkth2_bv_medrv_hl_roll_2020")

# set model names
model.names <- c(paste0("SV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("RSV-", c("N", "T", "AZ-SN", "AZ-ST", "FS-SN", "FS-ST", "GH-ST")),
                 paste0("EGARCH-", c("N", "T")),
                 paste0("REGARCH-", c("N", "T"))
                 )
model.n <- length(model.names)

# set path
path.data <- "../Data/"
path.out <- "./prediction_results_n225_2020/"

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
path.regarch <- "../REGARCH_Bayes/"

# set input path
path.in <- paste0("../Ox/", data.name, "_prediction_2020/")
path.in.az <- if (nzchar(Sys.getenv("RSVST_AZ", ""))) paste0(Sys.getenv("RSVST_AZ_DIR", "../Ox/az_recomputation/aggregated/"), sub("_excl$", "", Sys.getenv("RSVST_AZ")), "/", data.name, "_prediction_2020/") else path.in  # 2026-10-01 (§86): AZ 系 4 モデルの読み込み先
path.in.fs <- paste0("../Ox/", data.name, "_prediction_2020/")

# set true returns and RVs
data.true <- read_csv(paste0(path.data, true.sets, ".csv"))
date <- data.true$yyyymmdd
rcc <- data.true$rcc * 100 # パーセント・リターンに
rv5.hl <- data.true$rv5_hl * 10000
rk.th2.hl <- data.true$rk_th2_hl * 10000
bv.hl <- data.true$bv_hl * 10000
medrv.hl <- data.true$medrv_hl * 10000

# read prediction results for EGARCH and REGARCH models
f.egarch.n <- read_excel(paste0(path.regarch, "EGARCHn_", data.name, "_2020.xlsx"))
f.egarch.t <- read_excel(paste0(path.regarch, "EGARCHt_", data.name, "_2020.xlsx"))
f.regarch.n <- read_excel(paste0(path.regarch, "REGARCHn_", data.name, "_2020.xlsx"))
f.regarch.t <- read_excel(paste0(path.regarch, "REGARCHt_", data.name, "_2020.xlsx"))

# read and merge prediction results
# f.i <- 3
for (f.i in 1:f.n) {
  # set forecast variable
  f.name <- f.sets[f.i]
  
  # read prediction results for SV and RSV models
  f.sv.n <- read_csv(paste0(path.in, f.name, "_sv_az_n.csv"))
  f.sv.t <- read_csv(paste0(path.in, f.name, "_sv_az_t.csv"))
  f.sv.azsn <- read_csv(paste0(path.in.az, f.name, "_sv_az_sn.csv"))
  f.sv.azst <- read_csv(paste0(path.in.az, f.name, "_sv_az_st.csv"))
  f.sv.fssn <- read_csv(paste0(path.in.fs, f.name, "_sv_fs_sn.csv"))
  f.sv.fsst <- read_csv(paste0(path.in.fs, f.name, "_sv_fs_st.csv"))
  f.sv.ghst <- read_csv(paste0(path.in, f.name, "_sv_gh_st.csv"))
  f.rsv.n <- read_csv(paste0(path.in, f.name, "_rsv_az_n.csv"))
  f.rsv.t <- read_csv(paste0(path.in, f.name, "_rsv_az_t.csv"))
  f.rsv.azsn <- read_csv(paste0(path.in.az, f.name, "_rsv_az_sn.csv"))
  f.rsv.azst <- read_csv(paste0(path.in.az, f.name, "_rsv_az_st.csv"))
  f.rsv.fssn <- read_csv(paste0(path.in.fs, f.name, "_rsv_fs_sn.csv"))
  f.rsv.fsst <- read_csv(paste0(path.in.fs, f.name, "_rsv_fs_st.csv"))
  f.rsv.ghst <- read_csv(paste0(path.in, f.name, "_rsv_gh_st.csv"))

  # merge volatility forecasts
  if(f.i == 1){
    df.vol <- data.frame(date, rv5.hl, rk.th2.hl, bv.hl, medrv.hl, 
                         f.sv.n$Var3, f.sv.t$Var3, f.sv.azsn$Var3, f.sv.azst$Var3,
                         f.sv.fssn$MeanExp, f.sv.fsst$MeanExp, f.sv.ghst$Var3,
                         f.rsv.n$Var3, f.rsv.t$Var3, f.rsv.azsn$Var3, f.rsv.azst$Var3,
                         f.rsv.fssn$MeanExp, f.rsv.fsst$MeanExp, f.rsv.ghst$Var3,
                         f.egarch.n$`volatility forecast (mean)`, 
                         f.egarch.t$`volatility forecast (mean)`,
                         f.regarch.n$`volatility forecast (mean)`,
                         f.regarch.t$`volatility forecast (mean)`)
    names(df.vol) <- c("date", "RV5", "RKth2", "BV", "MedRV", model.names)
    df.vol.med <- data.frame(date, rv5.hl, rk.th2.hl, bv.hl, medrv.hl, 
                             f.sv.n$Var4, f.sv.t$Var4, f.sv.azsn$Var4, f.sv.azst$Var4, 
                             f.sv.fssn$MedianExp, f.sv.fsst$MedianExp, f.sv.ghst$Var4, 
                             f.rsv.n$Var4, f.rsv.t$Var4, f.rsv.azsn$Var4, f.rsv.azst$Var4,
                             f.rsv.fssn$MedianExp, f.rsv.fsst$MedianExp, f.rsv.ghst$Var4,
                             f.egarch.n$`volatility forecast (median)`, 
                             f.egarch.t$`volatility forecast (median)`,
                             f.regarch.n$`volatility forecast (median)`,
                             f.regarch.t$`volatility forecast (median)`)
    names(df.vol.med) <- c("date", "RV5", "RKth2", "BV", "MedRV", model.names)
  }
  
  # merge VaR forecasts (1% and 5%)
  if(f.i == 2){
    df.var.1p <- data.frame(date, rcc,
                            f.sv.n$Var2, f.sv.t$Var2, f.sv.azsn$Var2, f.sv.azst$Var2,
                            f.sv.fssn$VaR_1p, f.sv.fsst$VaR_1p, f.sv.ghst$Var2,
                            f.rsv.n$Var2, f.rsv.t$Var2, f.rsv.azsn$Var2, f.rsv.azst$Var2, 
                            f.rsv.fssn$VaR_1p, f.rsv.fsst$VaR_1p, f.rsv.ghst$Var2,
                            f.egarch.n$`1% VaR`, f.egarch.t$`1% VaR`,
                            f.regarch.n$`1% VaR`, f.regarch.t$`1% VaR`)
    names(df.var.1p) <- c("date", "Return", model.names)
    df.var.5p <- data.frame(date, rcc,
                            f.sv.n$Var3, f.sv.t$Var3, f.sv.azsn$Var3, f.sv.azst$Var3, 
                            f.sv.fssn$VaR_5p, f.sv.fsst$VaR_5p, f.sv.ghst$Var3,
                            f.rsv.n$Var3, f.rsv.t$Var3, f.rsv.azsn$Var3, f.rsv.azst$Var3,
                            f.rsv.fssn$VaR_5p, f.rsv.fsst$VaR_5p, f.rsv.ghst$Var3,
                            f.egarch.n$`5% VaR`, f.egarch.t$`5% VaR`,
                            f.regarch.n$`5% VaR`, f.regarch.t$`5% VaR`)
    names(df.var.5p) <- c("date", "Return", model.names)
  }
  
  # merge ES forecasts (1% and 5%)
  if(f.i == 3){
    df.es.1p <- data.frame(date, rcc,
                           f.sv.n$Var2, f.sv.t$Var2, f.sv.azsn$Var2, f.sv.azst$Var2,
                           f.sv.fssn$ES_1p, f.sv.fsst$ES_1p, f.sv.ghst$Var2,
                           f.rsv.n$Var2, f.rsv.t$Var2, f.rsv.azsn$Var2, f.rsv.azst$Var2, 
                           f.rsv.fssn$ES_1p, f.rsv.fsst$ES_1p, f.rsv.ghst$Var2,
                           f.egarch.n$`1% ES`, f.egarch.t$`1% ES`,
                           f.regarch.n$`1% ES`, f.regarch.t$`1% ES`)
    names(df.es.1p) <- c("date", "Return", model.names)
    df.es.5p <- data.frame(date, rcc,
                           f.sv.n$Var3, f.sv.t$Var3, f.sv.azsn$Var3, f.sv.azst$Var3,
                           f.sv.fssn$ES_5p, f.sv.fsst$ES_5p, f.sv.ghst$Var3,
                           f.rsv.n$Var3, f.rsv.t$Var3, f.rsv.azsn$Var3, f.rsv.azst$Var3, 
                           f.rsv.fssn$ES_5p, f.rsv.fsst$ES_5p, f.rsv.ghst$Var3,
                           f.egarch.n$`5% ES`, f.egarch.t$`5% ES`,
                           f.regarch.n$`5% ES`, f.regarch.t$`5% ES`)
    names(df.es.5p) <- c("date", "Return", model.names)
  }
}

# save forecasts
# 2026-10-11 (IJF R3, PROGRESS.md §100-101): RSVST_MASK_DIR=<dir> drops the rows whose date is listed in <dir>/<NAME>_exclude_dates.csv
# (NAME = DJIA | N225 | DJIA_2020 | N225_2020) from every merged data frame, so that all models and proxies are evaluated on the same
# reduced sample (exclusion sensitivity for the windows that did not pass the MCMC confirmation pass; rule §7-6). Used with RSVST_AZ=<tag>_excl.
.mask.dir <- Sys.getenv("RSVST_MASK_DIR", "")
if (nzchar(.mask.dir)) {
  .mask.file <- file.path(.mask.dir, paste0(paste0(data.name, "_2020"), "_exclude_dates.csv")); .mask <- as.integer(read.csv(.mask.file)$date); .n0 <- nrow(df.vol)
  for (.d in c("df.vol", "df.vol.med", "df.var.1p", "df.var.5p", "df.es.1p", "df.es.5p")) if (exists(.d)) assign(.d, get(.d)[!(as.integer(get(.d)$date) %in% .mask), ])
  cat(sprintf("RSVST_MASK_DIR: %s, excluded %d of %d rows (%d dates listed)\n", .mask.file, .n0 - nrow(df.vol), .n0, length(.mask)))
}
write_csv(df.vol, paste0(path.out, data.name, "_forecast_vol.csv"))
write_csv(df.vol.med, paste0(path.out, data.name, "_forecast_vol_med.csv"))
write_csv(df.var.1p, paste0(path.out, data.name, "_forecast_var1p.csv"))
write_csv(df.var.5p, paste0(path.out, data.name, "_forecast_var5p.csv"))
write_csv(df.es.1p, paste0(path.out, data.name, "_forecast_es1p.csv"))
write_csv(df.es.5p, paste0(path.out, data.name, "_forecast_es5p.csv"))
