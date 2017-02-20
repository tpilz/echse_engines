library(xts)

setwd("/home/tobias/Dokumente/test_soilwat/run/analysis/")

# plot name names prefix
plot_file <- "Green-Ampt_test"

# model debug file
res_file <- "/home/tobias/Dokumente/test_soilwat/run/out/test1.dbg"

# horizon parameters
wc_par_file <- "/home/tobias/Dokumente/test_soilwat/data/hor_pars.dat"

# which dates to plot
#dates <- as.POSIXct(c("2001-01-01 01:00:00", "2001-01-01 02:00:00"), format="%Y-%m-%d %H:%M:%S", tz="UTC")



# read data
res_dat <- read.table(res_file, header=T, sep="\t")
wc_par_dat <- read.table(wc_par_file, header=T, sep="\t")

wc_hor_dat <- res_dat[res_dat$item_name == "wc",]
wc_hor_xts <- xts(as.numeric(wc_hor_dat$value), as.POSIXct(wc_hor_dat$end_of_interval, tz='UTC'))

dates <- unique(index(wc_hor_xts))

# plot
for (i in 1:length(dates)) {
  date <- dates[i]
  png(paste(plot_file, date, sep="_"))
  # initialise
  plot(0,0, type="n", xlim=c(0,1), ylim=c(-sum(wc_par_dat$hor_depth),0), ylab="Depth (m)", xlab="Water content (m3/m3)")
  
  # plot layers
  abline(h=-cumsum(wc_par_dat$hor_depth), lty=2)
  abline(h=0, lty=1)
  
  # plot parameters
  segments(x0=wc_par_dat$wc_res, y0=c(0,-head(cumsum(wc_par_dat$hor_depth),n=-1)), y1=-cumsum(wc_par_dat$hor_depth), col="red")
  segments(x0=wc_par_dat$wc_s, y0=c(0,-head(cumsum(wc_par_dat$hor_depth),n=-1)), y1=-cumsum(wc_par_dat$hor_depth), col="blue")
  segments(x0=wc_par_dat$wc_fc, y0=c(0,-head(cumsum(wc_par_dat$hor_depth),n=-1)), y1=-cumsum(wc_par_dat$hor_depth), col="green", lty=2)
  
  # plot actual water content
  segments(x0=wc_hor_xts[date,], y0=c(0,-head(cumsum(wc_par_dat$hor_depth),n=-1)), y1=-cumsum(wc_par_dat$hor_depth), col="black")
  dev.off()
}
