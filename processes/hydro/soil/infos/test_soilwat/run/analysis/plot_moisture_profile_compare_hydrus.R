library(xts)

setwd("/home/tobias/Dokumente/test_soilwat/run/analysis/")

# plot name names prefix
plot_file <- "compare_hydrus"

# model debug file
res_file <- "/home/tobias/Dokumente/test_soilwat/run/output_daily_res/test1.dbg"

# horizon parameters
wc_par_file <- "/home/tobias/Dokumente/test_soilwat/data/hor_pars.dat"

# initials
inits_file <- "/home/tobias/Dokumente/test_soilwat/data/init_vect.dat"

# hydrus 1d results, soil moisture
hydrus_states_file <- "/home/tobias/Dokumente/test_soilwat/hydrus_results/moisture_states_init_fc_long.dat"

# hydrus 1d results, fluxs
hydrus_flux_file <- "/home/tobias/Dokumente/test_soilwat/hydrus_results/flux_init_fc_long.dat"

# which dates to plot
#dates <- as.POSIXct(c("2001-01-01 01:00:00"), format="%Y-%m-%d %H:%M:%S", tz="UTC")#, "2001-01-01 02:00:00"), format="%Y-%m-%d %H:%M:%S", tz="UTC")



# read data
res_dat <- read.table(res_file, header=T, sep="\t")
wc_par_dat <- read.table(wc_par_file, header=T, sep="\t")
hydrus_states_dat <- read.table(hydrus_states_file, header = T, sep = "\t")
hydrus_flux_dat <- read.table(hydrus_flux_file, header = T, sep = "\t")
inits_dat <- read.table(inits_file, header=T, sep="\t")

# echse water content per horizon
start_date <- as.character(rep(format(as.POSIXct(head(res_dat$end_of_interval[grep("^wc$", res_dat$item_name)], n=1), tz='UTC')-86400, format="%Y-%m-%d %H:%M:%S"), nrow(wc_par_dat)))
echse_states_dat <- data.frame(date = c(start_date, as.character(res_dat$end_of_interval[grep("^wc$", res_dat$item_name)])),
                               wc = c(inits_dat$value[grep("^wc$", inits_dat$variable)], res_dat$value[grep("^wc$", res_dat$item_name)]),
                               mat_pot = c(inits_dat$value[grep("^mat_pot$", inits_dat$variable)], res_dat$value[grep("^mat_pot$", res_dat$item_name)]),
                               k_u = c(inits_dat$value[grep("^k_u$", inits_dat$variable)], res_dat$value[grep("^k_u$", res_dat$item_name)])*1000*3600)
echse_states_xts <- xts(echse_states_dat[,-1], as.POSIXct(echse_states_dat$date, tz='UTC', format="%Y-%m-%d %H:%M:%S"))

# echse fluxes
fac <- 1000*3600 # factor to convert units from ms-1 -> mmh-1
vol_wat_init <- xts(sum(inits_dat$value[grep("^wc$", inits_dat$variable)] * wc_par_dat$hor_depth)*1000, as.POSIXct(head(index(echse_states_xts),n=1), tz='UTC'))
echse_flux <- data.frame(run_surf = c(0,res_dat$value[grep("^run_surf$", res_dat$item_name)]*fac),
                         inflow = c(0,res_dat$value[grep("^inflow$", res_dat$item_name)]*fac),
                         vol_wat = c(vol_wat_init, tapply(res_dat$value[grep("^wc$", res_dat$item_name)], res_dat$end_of_interval[grep("^wc$", res_dat$item_name)], function(x) sum(x*wc_par_dat$hor_depth))*1000),
                         perc_out = c(0,res_dat$value[grep("^run_gw$", res_dat$item_name)]*fac))
echse_flux$inf_sum <- cumsum(echse_flux$inflow - echse_flux$run_surf)
echse_flux_xts <- xts(echse_flux, unique(index(echse_states_xts)))

# hydrus water content per layer
hydrus_states_xts <- xts(hydrus_states_dat[,-1], as.POSIXct(hydrus_states_dat$date, tz='UTC'))

# hydrus fluxes, volumetric states (whole soil profile)
hydrus_flux <- hydrus_flux_dat[,c("run_surf_mh.1", "inf_sum_m", "vol_wat_m", "perc_out_mh.1")]
hydrus_flux_xts <- xts(hydrus_flux, as.POSIXct(hydrus_flux_dat[,3], tz='UTC'))


dates <- unique(index(hydrus_states_xts))
#dates <- grep("00:00:00", dates, value = T) # only daily values in this case


for (i in 1:length(dates)) {
  date <- dates[i]
  
  # plot moisture states
  png(paste(plot_file, "wc", date, sep="_"), width=650, height=480)
  par(oma=c(0,0,0,10), mar=c(5,4,1,0))
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
  segments(x0=echse_states_xts[date,"wc"], y0=c(0,-head(cumsum(wc_par_dat$hor_depth),n=-1)), y1=-cumsum(wc_par_dat$hor_depth), col="black")
  lines(as.vector(hydrus_states_xts[date,"moisture"]),-seq(0,0.01*(length(hydrus_states_xts[date,"moisture"])-1), by=0.01), lty=2)
  
  # legend
  legend(1.06,0,legend = c("ECHSE", "HYDRUS-1D", "Residual WC", "Saturated WC", "Field capacity"), lty=c(1,2,1,1,1), 
         col=c("black", "black", "red", "blue", "green"), xpd=NA)
  dev.off()
  
  
  
  # plot matric potential
  png(paste(plot_file, "mat_pot", date, sep="_"), width=650, height=480)
  par(oma=c(0,0,0,10), mar=c(5,4,1,0))
  # initialise
  plot(0,0, type="n", xlim=c(0,5), ylim=c(-sum(wc_par_dat$hor_depth),0), ylab="Depth (m)", xlab="Matric potential (m)/(100 hPa)")
  
  # plot layers
  abline(h=-cumsum(wc_par_dat$hor_depth), lty=2)
  abline(h=0, lty=1)
  
  # plot actual water content
  segments(x0=echse_states_xts[date,"mat_pot"], y0=c(0,-head(cumsum(wc_par_dat$hor_depth),n=-1)), y1=-cumsum(wc_par_dat$hor_depth), col="black")
  lines(-1*as.vector(hydrus_states_xts[date,"head"]),-seq(0,0.01*(length(hydrus_states_xts[date,"head"])-1), by=0.01), lty=2)

  # legend
  legend(5.3,0,legend = c("ECHSE", "HYDRUS-1D"), lty=c(1,2), xpd=NA)
  dev.off()
  
  
  
  # plot conductivity
  png(paste(plot_file, "cond", date, sep="_"), width=650, height=480)
  par(oma=c(0,0,0,10), mar=c(5,4,1,0))
  # initialise
  plot(0,0, type="n", xaxt="n", xlim=log10(c(1e-6,100)), ylim=c(-sum(wc_par_dat$hor_depth),0), ylab="Depth (m)", xlab="")
  
  ats = log10(c(1e-6, 1e-4, 1e-2, 1, 100))
  labs = 10^ats
  axis(1, at=ats, labels=labs, las=1)
  title(xlab="Hydraulic conductivity (mm/h)", line=3.5)
  
  # plot layers
  abline(h=-cumsum(wc_par_dat$hor_depth), lty=2)
  abline(h=0, lty=1)
  
  # plot ksat
  segments(x0=log10(wc_par_dat$k_s*1000*3600), y0=c(0,-head(cumsum(wc_par_dat$hor_depth),n=-1)), y1=-cumsum(wc_par_dat$hor_depth), col="blue")
  
  # plot actual water content
  segments(x0=log10(echse_states_xts[date,"k_u"]), y0=c(0,-head(cumsum(wc_par_dat$hor_depth),n=-1)), y1=-cumsum(wc_par_dat$hor_depth), col="black")
  lines(log10(as.vector(hydrus_states_xts[date,"ku"])*1000),-seq(0,0.01*(length(hydrus_states_xts[date,"ku"])-1), by=0.01), lty=2)
  
  # legend
  legend(log10(250),0,legend = c("ECHSE", "HYDRUS-1D", "Satur. conduct."), lty=c(1,2,1), col=c("black", "black", "blue"), xpd=NA)
  dev.off()
}





# plot time series

# surface runoff
png(paste(plot_file, "ts_surf", sep="_"), width=650, height=480)
par(oma=c(0,0,0,10), mar=c(5,4,1,0))
plot(echse_flux_xts[,"run_surf"], ylim=c(0,max(echse_flux_xts[,"run_surf"],hydrus_flux_xts[,"run_surf_mh.1"]*1000)), type="l", main="", xlab="Time", ylab="Surface runoff [mm/h]")
lines(hydrus_flux_xts[,"run_surf_mh.1"]*1000, col="red")
legend(max(index(echse_flux_xts))+10*86400,80,legend = c("ECHSE", "HYDRUS-1D"), col=c("black", "red"), lty=1, xpd=NA)
dev.off()

# infiltration (cumulated)
png(paste(plot_file, "ts_inf", sep="_"), width=650, height=480)
par(oma=c(0,0,0,10), mar=c(5,4,1,0))
plot(echse_flux_xts[,"inf_sum"], ylim=c(0,max(hydrus_flux_xts[,"inf_sum_m"]*1000,echse_flux_xts[,"inf_sum"])), type="l", main="", xlab="Time", ylab="Cumulated infiltration [mm]")
lines(hydrus_flux_xts[,"inf_sum_m"]*1000, col="red")
legend(max(index(echse_flux_xts))+10*86400,60,legend = c("ECHSE", "HYDRUS-1D"), col=c("black", "red"), lty=1, xpd=NA)
dev.off()

# water in soil profile
png(paste(plot_file, "ts_volwat", sep="_"), width=650, height=480)
par(oma=c(0,0,0,10), mar=c(5,4,1,0))
plot(echse_flux_xts[,"vol_wat"], type="l", ylim=c(0,max(hydrus_flux_xts[,"vol_wat_m"]*1000,echse_flux_xts[,"vol_wat"])), main="", xlab="Time", ylab="Water in soil profile [mm]")
lines(hydrus_flux_xts[,"vol_wat_m"]*1000, col="red")
legend(max(index(echse_flux_xts))+10*86400,250,legend = c("ECHSE", "HYDRUS-1D"), col=c("black", "red"), lty=1, xpd=NA)
dev.off()

# percolation rate from bottom of profile
png(paste(plot_file, "ts_percout", sep="_"), width=650, height=480)
par(oma=c(0,0,0,10), mar=c(5,4,1,0))
plot(echse_flux_xts[,"perc_out"], ylim=c(0,max(echse_flux_xts[,"perc_out"],hydrus_flux_xts[,"perc_out_mh.1"]*-1000)), type="l", main="", xlab="Time", ylab="Percolation [mm/h]")
lines(hydrus_flux_xts[,"perc_out_mh.1"]*-1000, col="red")
legend(max(index(echse_flux_xts))+10*86400,0.06,legend = c("ECHSE", "HYDRUS-1D"), col=c("black", "red"), lty=1, xpd=NA)
dev.off()
