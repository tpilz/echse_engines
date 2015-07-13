rm(list=ls())

# BEGIN SETTINGS ###############################################################

station= "Fichtelberg"
#station= "KahlerAsten"

dir_base= "/home/dkneis/progress/echse_core/processes/snow/modeltest"

ifile_sim= paste(dir_base,"/out/",station,"/snow_sim.txt",sep="")
ifile_dbg= paste(dir_base,"/out/",station,"/snow_dbg.txt",sep="")

ifile_obs= paste(dir_base,"/data/DWD/raw/",station,"/snowwe.txt",sep="")

ofile= paste(dir_base,"/out/",station,"/show.pdf",sep="")

# END SETTINGS #################################################################

sim= read.table(file= ifile_sim, sep="\t", header=T)
dbg= read.table(file= ifile_dbg, sep="\t", header=T)
dat= merge(x=sim, y=dbg, by="time", all=T)
rm(sim)
rm(dbg)

obs= read.table(file= ifile_obs, sep="\t", header=T)
names(obs)= c("time","swe_obs")
obs$swe_obs= obs$swe_obs/1000
dat= merge(x=dat, y=obs, by="time", all.x=T, all.y=F)
rm(obs)

dat$time= as.POSIXct(strptime(dat$time,"%Y-%m-%d"))

pdf(file= ofile, width=6, height=4, onefile=T, pointsize=10)

scale= function(x, x_zero, x_one) {
  return((x-x_zero)/(x_one-x_zero))
}

ax= function(x, side, line, col) {
  z= pretty(range(x, na.rm=T))
  axis(side=side, line=line, col=col, at=scale(z, min(x, na.rm=T), max(x, na.rm=T)),
    labels=z)
}

# Plot states
op= par(no.readonly=T)
par(mar=c(4,4,1,6))
plot(range(dat$time), c(0,1.1), type="n", ylab="", yaxt="n", bty="n")
lines(dat$time, scale(dat$alb, min(dat$alb), max(dat$alb)), col="grey")
  ax(x=dat$alb, side=4, line=2.5, col="grey")
lines(dat$time, scale(dat$sec, min(dat$sec), max(dat$sec)), col="red")
  ax(x=dat$sec, side=4, line=0, col="red")
mi= min(c(dat$swe,dat$swe_obs),na.rm=T)
mx= max(c(dat$swe,dat$swe_obs),na.rm=T)
lines(dat$time, scale(dat$swe, mi, mx), col="blue")
points(dat$time, scale(dat$swe_obs, mi, mx), col="blue", pch=20)
  ax(x=dat$swe, side=2, line=0, col="blue")
legend("top",bty="n",lty=c(1,1,1,0), pch=c(NA,NA,NA,20), horiz=T,
  col=c("blue","red","grey","blue"), legend=c("swe","sec","alb","swe_obs"))
par(op)

# Plot derived values
op= par(no.readonly=T)
par(mar=c(4,4,1,6))
plot(range(dat$time), c(0,1.1), type="n", ylab="", yaxt="n", bty="n")
lines(dat$time, scale(dat$temp_mean, min(dat$temp_mean,na.rm=T), max(dat$temp_mean,na.rm=T)), col="red")
  ax(x=dat$temp_mean, side=4, line=2.5, col="red")
lines(dat$time, scale(dat$temp_surf, min(dat$temp_surf,na.rm=T), max(dat$temp_surf,na.rm=T)), col="orange")
  ax(x=dat$temp_surf, side=4, line=0, col="orange")
lines(dat$time, scale(dat$liqu_frac, min(dat$liqu_frac,na.rm=T), max(dat$liqu_frac,na.rm=T)), col="grey")
  ax(x=dat$liqu_frac, side=2, line=0, col="grey")
legend("top",bty="n",lty=c(1,1,1), horiz=T,
  col=c("red","orange","grey"), legend=c("T_mean","T_surf","LiquFrac"))
par(op)

# Plot swe for winters
op= par(no.readonly=T)
par(mar=c(4,2.5,1.2,0.5))
years= as.integer(unique(format(dat$time,"%Y")))
nc= ceiling(length(years)/2)
layout(matrix(1:(nc*2),ncol=nc,byrow=T),width=1,height=1)
mi= min(c(dat$swe,dat$swe_obs),na.rm=T)
mx= max(c(dat$swe,dat$swe_obs),na.rm=T)
for (i in 1:length(years)) {
  # Select data
  t0= ISOdatetime(years[i],11,1,0,0,0)
  t1= ISOdatetime(years[i]+1,5,15,0,0,0)
  # Empty plot
  plot(c(t0,t1), c(mi,mx), type="n", ylab="", main=paste(years[i],"/",years[i]+1,sep=""))
  # Plot background data
  inds= which(dat$windSpeed >= quantile(dat$windSpeed,probs=0.99,na.rm=T))
  abline(v=dat$time[inds], col="lightgrey")
  # Plot swe
  points(dat$time, dat$swe_obs, pch=20, col="darkblue")
  lines(dat$time, dat$swe, col="blue")
}
par(op)

# Plot meteo vars
layout(matrix(1:8,nrow=2,ncol=4,byrow=T),width=1,height=1)
plot(dat$time, dat$precipSumMM, type="l", col="blue", ylab="", main="PRECIP")
plot(dat$time, dat$shortRad, type="l", col="blue", ylab="", main="GLORAD")
plot(dat$time, dat$tempAir, type="l", col="blue", ylab="", main="TEMPER")
plot(dat$time, dat$pressAir, type="l", col="blue", ylab="", main="APRESS")
plot(dat$time, dat$relHumid, type="l", col="blue", ylab="", main="RHUMID")
plot(dat$time, dat$windSpeed, type="l", col="blue", ylab="", main="WINDSP")
plot(dat$time, dat$cloudCoverage, type="l", col="blue", ylab="", main="CLNESS")

# Plot rates and stoi factors
layout(matrix(1:12,nrow=3,ncol=4,byrow=T),width=1,height=1)
plot(dat$time, dat$flux_R_netS, type="l", col="blue", ylab="", main="flux_R_netS")
plot(dat$time, dat$flux_R_netL, type="l", col="blue", ylab="", main="flux_R_netL")
plot(dat$time, dat$flux_R_soil, type="l", col="blue", ylab="", main="flux_R_soil")
plot(dat$time, dat$flux_R_sens, type="l", col="blue", ylab="", main="flux_R_sens")
plot(dat$time, dat$flux_M_prec, type="l", col="blue", ylab="", main="flux_M_prec")
plot(dat$time, dat$flux_M_subl, type="l", col="blue", ylab="", main="flux_M_subl")
plot(dat$time, dat$flux_M_flow, type="l", col="blue", ylab="", main="flux_M_flow")
plot(dat$time, dat$rate_G_alb, type="l", col="blue", ylab="", main="rate_G_alb")
plot(dat$time, dat$stoi_f_prec, type="l", col="blue", ylab="", main="stoi_f_prec")
plot(dat$time, dat$stoi_f_subl, type="l", col="blue", ylab="", main="stoi_f_subl")
plot(dat$time, dat$stoi_f_flow, type="l", col="blue", ylab="", main="stoi_f_flow")

graphics.off()

