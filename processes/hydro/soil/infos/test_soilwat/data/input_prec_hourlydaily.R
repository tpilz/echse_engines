library(xts) 

file <- "inputs_ts.dat"

dat <- read.table(file, sep="\t", header=T)
dat_xts <- xts(dat$test_loc, as.POSIXct(dat$end_of_interval, timezone="UTC"))

# calculate daily sums
dat_daily <- apply.daily(dat_xts, mean)

write.table(dat_daily, "inputs_ts_daily.dat", col.names = T, row.names = format(index(dat_daily), "%Y-%m-%d 00:00:00"), quote=F, sep="\t")
