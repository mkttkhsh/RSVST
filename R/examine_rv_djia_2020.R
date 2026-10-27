rm(list = ls())

# packages
library(readr)
library(ggplot2)
library(reshape2)
library(psych)
library(tseries)

# functions
source("myplot.R")
source("LB_Diebold.R")

# path
path.figure <- "./figure_djia_2020/"
path.out <- "./stats/"

# set data
data.name <- "DJIA"

# select period
beg.y <- 2012
end.y <- 2020

# read data
data0 <- read_csv(paste0("../Data/", data.name, "_rcc_rv5_rkth2_bv_medrv_", beg.y, "_", end.y, "_tick.csv"), skip_empty_rows = TRUE)
data0 <- subset(data0, !is.na(data0$yyyymmdd))
data0$rcc <- data0$rcc * 100
data0$rv <- data0$rv5 * 10000

# edit dates
yyyymmdd <- data0$yyyymmdd
data0$date <- as.Date(as.character(yyyymmdd), format = '%Y%m%d')
data1 <- data.frame(yyyymmdd, date = data0$date, data0[, 7])
data1$lnrv = log(data1$rv)

# descriptive statistics
des <- describe(data1[, 3:4], type = 1)
des$kurtosis <- des$kurtosis + 3 #make kurtosis representing E(x^4)/E(x^2)^2
n <- nrow(data1)
se.skew.kurt <- c(sqrt(6 / n), sqrt(24 / n))
stats <- data.frame(des$mean, des$se, des$sd, des$skew, se.skew.kurt[1], des$kurtosis, se.skew.kurt[2], des$min, des$max)
names(stats) <- c("Mean", "SE", "SD", "Skew", "SE", "Kurt", "SE", "Min", "Max")

# Jarque-Bera test
jb.rcc <- jarque.bera.test(x = data1[, 3])
stats <- cbind(stats, rbind(jb.rcc$p.value))
names(stats)[ncol(stats)] <- "JB"

# Ljung-Box test
lb.rcc <- LB_Diebold(data1[, 3], 10)
stats <- cbind(stats, rbind(lb.rcc))
names(stats)[ncol(stats)] <- "LB"

# save
row.names(stats) <- c("RV", "logRV")
write_csv(stats, paste0(path.out, data.name, "_", beg.y, "_", end.y, "_stats_rv.csv"))
print(stats, digits = 3)

# output for latex
stats.mat <- as.matrix(stats)
nr <- nrow(stats.mat)
nc <- ncol(stats.mat)
for (i in 1:nr) {
  xxx <- rownames(stats.mat)[i]
  for (j in 1:nc) {
    if (j == 2 | j == 5 | j == 7) {
      xxx <- paste0(xxx, " & $(", format(round(stats.mat[i,j], 3), nsmall = 3), ")$")
    } else if (j >= 10) {
      xxx <- paste0(xxx, " & $", format(round(stats.mat[i,j], 2), nsmall = 2), "$")
    } else {
      xxx <- paste0(xxx, " & $", format(round(stats.mat[i,j], 3), nsmall = 3), "$")
    }
  }
  xxx <- paste0(xxx, " \\\\")
  cat(xxx, "\n")
}
# xtable(stats.mat, label="tab:stats", caption="Summary statistics", digits = 3)

# melt data
# data.melt <- melt(data1[, 2:7], id.vars = "date")

# time-series plot of 5-min RV
p <- ggplot(data1, aes(date, rv)) + geom_line() + xlab("") + ylab("") + theme_bw() +
  theme(axis.text=element_text(size=20)) + scale_x_date(date_breaks = "1 year", date_labels = "%Y")
print(p)
ggsave(file = paste0(path.figure, data.name, "_", beg.y, "_", end.y, "_rv_tsplot.png"), plot = p, dpi = 100, width=16, height=9)

# time-series plot of log of 5-min RV
p <- ggplot(data1, aes(date, lnrv)) + geom_line() + xlab("") + ylab("") + theme_bw() +
  theme(axis.text=element_text(size=20)) + scale_x_date(date_breaks = "1 year", date_labels = "%Y")
print(p)
ggsave(file = paste0(path.figure, data.name, "_", beg.y, "_", end.y, "_lnrv_tsplot.png"), plot = p, dpi = 100, width=16, height=9)

# histogram of log of 5-min RV
bw <- diff(range(data1$lnrv)) / 20
p <- ggplot(data1, aes(x = lnrv)) + geom_histogram(aes(y=..density..), binwidth = bw, fill = 'white', color = 'black') +
  geom_density(fill='black', alpha = 0.3) +
  xlab("") + ylab("") + theme_bw()
# + xlim(range(dens$x))
print(p)
ggsave(file = paste0(path.figure, data.name, "_", beg.y, "_", end.y, "_lnrv_histogram.png"), plot = p, dpi = 100, width=16, height=9)
