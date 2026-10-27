rm(list = ls())

# packages
library(readr)
library(ggplot2)

# download data
curl <- "https://realized.oxford-man.ox.ac.uk/images/oxfordmanrealizedvolatilityindices.zip"
cdestfile <- "realizedlibrary.zip"
download.file(curl,cdestfile)

# read data
data0 <- read_csv(cdestfile)
table(data0$Symbol)

# extract dates
date <- substr(data0$...1, 1, 10)
yyyymmdd <- as.numeric(sub("-", "", sub("-", "", date)))
date <- as.Date(date)
data1 <- data.frame(date, yyyymmdd, data0)

# extract 1 index
index <- ".SPX"
data2 <- subset(data1, data1$Symbol == index)

# times-series plot
p1 <- ggplot(data2, aes(date, nobs)) + geom_line() + xlab("") + ylab("Nobs")
p2 <- ggplot(data2, aes(date, open_to_close)) + geom_line() + xlab("") + ylab("Open-to-close Return")
p3 <- ggplot(data2, aes(date, rv5)) + geom_line() + xlab("") + ylab("5-min RV")
p4 <- ggplot(data2, aes(date, rv5_ss)) + geom_line() + xlab("") + ylab("5-min RV (sub-sampled)") # same as the 5-min RV
p5 <- ggplot(data2, aes(date, rk_parzen)) + geom_line() + xlab("") + ylab("RK (non-flat Parzen)")
p6 <- ggplot(data2, aes(date, rk_th2)) + geom_line() + xlab("") + ylab("RK (two-scale / Bartlet)")
p7 <- ggplot(data2, aes(date, rv10)) + geom_line() + xlab("") + ylab("10-min RV")
p8 <- ggplot(data2, aes(date, medrv)) + geom_line() + xlab("") + ylab("MedRV")
gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8)

+ 
  annotate("rect", xmin=as.Date("2010-01-01"), xmax=as.Date("2019-12-31"), ymin=0, ymax=max(data2$rv5), 
           alpha=0.1, fill="blue") +
  annotate("rect", xmin=as.Date("2017-01-01"), xmax=as.Date("2019-12-31"), ymin=0, ymax=max(data2$rv5), 
           alpha=0.2, fill="red")

# number of observations in 2010-2019
data3 <- subset.data.frame(data2, data2$yyyymmdd>=20100101 & data2$yyyymmdd<=20191231, c(date, yyyymmdd, rv5))
nall <- length(data3$rv5)
nblue <- length(data3$rv5[data3$yyyymmdd<20170101])
nred <- length(data3$rv5[data3$yyyymmdd>=20170101])
data.frame(All=nall, Blue=nblue, Red=nred)
