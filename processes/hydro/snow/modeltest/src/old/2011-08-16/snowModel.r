rm(list=ls())

# BEGIN SETTINGS ###############################################################

# Load library function
source("/home/dkneis/cpp/dklib/hydro_snow/r_wrapper/snow_wrapper.r")
source("/home/dkneis/progress/echse_tools/precipcorr/precipcorr.r")

dir_base= "/home/dkneis/progress/echse_core/processes/snow/modeltest"

station= "Fichtelberg"
#station= "KahlerAsten"

dir_meteo= paste(dir_base,"/data/DWD/use/",station,sep="")
ifile_obs= paste(dir_base,"/data/DWD/raw/",station,"/snowwe.txt",sep="")

ofile_sim= paste(dir_base,"/out/",station,"/snow_sim.txt",sep="")
ofile_dbg= paste(dir_base,"/out/",station,"/snow_dbg.txt",sep="")
replace= TRUE

# Parameters
params= c(
  precipSeconds= 86400, # for daily meteo data
  a0= 2./(92.6*3.6*1.29) /10,  # --> to be calibrated
  a1= 1.6/(92.6*3.6*1.29) /10, # --> to be calibrated
  kSatSnow= 0.02/3600,   # --> to be calibrated
  densDrySnow= 450.,
  specCapRet= 0.05,
  emissivitySnowMin= 0.84,
  emissivitySnowMax= 0.99,
  tempAir_crit= 0.,
  albedoMin= 0.4,
  albedoMax= 0.85,
  agingRate_tAirPos= -0.12/86400* log(0.85-0.4),
  agingRate_tAirNeg= -0.05/86400* log(0.85-0.4),
  soilDepth= 0.1,
  soilDens= 0.6*1500.+0.4*1000.,    # Saturated soil
  soilSpecHeat= 0.6*0.85+0.4*4.18,  # Saturated soil
  weightAirTemp= 0.3 # --> to be calibrated
)

# Parameters for precipitation correction
height_zeroWind= 0.4
if (station == "Fichtelberg") {
  sensorHeight_u= 29.
} else if (station == "KahlerAsten") {
  sensorHeight_u= 27.3
} else {
  sensorHeight_u= 10.
}
sensorHeight_p= 1.

# Initial values of the states
if (station == "Fichtelberg") {
  initials= c(sec= 0., swe= 0.45, alb= 0.8)
} else if (station == "KahlerAsten") {
  initials= c(sec= 0., swe= 0.21, alb= 0.8)
} else {
  stop("Station unknown.")
}

# Forcings
meteovars= c("precip", "glorad", "temper", "apress", "rhumid", "windsp", "clness")

# Number of sub-time steps per day
n_dtsub= 24

# END SETTINGS #################################################################

# Read forcings and merge into a single table
print("Reading meteo data...")
for (i in 1:length(meteovars)) {
  ifile= paste(dir_meteo,"/",meteovars[i],".txt",sep="")
  tmp= read.table(file=ifile, header=T, sep="\t", colClasses=c("POSIXct","numeric"))
  names(tmp)= c("time",meteovars[i])
  if (i == 1) {
    meteo= tmp
  } else {
    meteo= merge(x=meteo, y=tmp, by="time", all=T)
  }
}
if (any(is.na(meteo[,which(names(meteo) != "time")]))) stop("NAs in meteo data.")
meteo= meteo[sort.list(meteo$time),]

# Correct precipitation for wind error
meteo$precip= meteo$precip * factPrecipCorr(meteo$windsp, meteo$temper,
  sensorHeight_p, sensorHeight_u, height_zeroWind)
# Convert windspeed to values as standard height (10m)
meteo$windsp= meteo$windsp * log(10/height_zeroWind)/log(sensorHeight_p/height_zeroWind)

# Set simulation window
tmin= meteo$time[1]
tmax= meteo$time[nrow(meteo)]
delta_t= as.numeric(difftime(tmax,tmin,units="secs")/(nrow(meteo)-1))
print(paste("Simulation window: ",tmin," - ",tmax," (",
  as.numeric(difftime(tmax,tmin,units="days"))," days)",sep=""))
print(paste("Delta t (seconds): ",delta_t,sep=""))

# Init output files
if (file.exists(ofile_sim) && (!replace)) stop("Output already exists.")
ovect= c("time",names(initials))
write(x=ovect, file=ofile_sim, ncolumns=length(ovect), sep="\t", append=F)

# Time loop
states= initials
for (i in 1:nrow(meteo)) {
#for (i in 1:1280) {
  print(paste(meteo$time[i]," (",round(i/nrow(meteo)*100,1),"%)",sep=""))

  # Merge constant parameters with forcings for current time step
  met= c(
    precipSumMM= meteo$precip[i],
    shortRad= meteo$glorad[i],
    tempAir= meteo$temper[i],
    pressAir= meteo$apress[i],
    relHumid= meteo$rhumid[i],
    windSpeed= meteo$windsp[i],
    cloudCoverage= meteo$clness[i]
  )
  p= c(params, met)

  # Print debug info
  dbg= snowDebug(y=states, p=p)
  if (i == 1) {
    if (file.exists(ofile_dbg) && (!replace)) stop("Output already exists.")
    ovect= c("time",names(met),names(dbg))
    write(x=ovect, file=ofile_dbg, ncolumns=length(ovect), sep="\t", append=F)
  }
  ovect= c(format(meteo$time[i],"%Y-%m-%d"),sprintf("%0.2f",met),sprintf("%0.3e",dbg))
  write(x=ovect, file=ofile_dbg, ncolumns=length(ovect), sep="\t", append=T)

  # Integrate and update states
  for (k in 1:n_dtsub) {
    ddt_states= snowDerivs(t=NA, y=states, p=p)
    ddt_states= unlist(ddt_states)
    # Correct if SWE would become < 0
    if ((ddt_states["ddt_swe"] < 0.) &&
        (abs(ddt_states["ddt_swe"]*delta_t/n_dtsub) > states["swe"])) {
      states["swe"]= 0.0001 # this is 1/10 mm
      states["sec"]= 0.
      states["alb"]= as.numeric(params["albedoMax"])
      break
    } else {
      states= states + ddt_states * delta_t/n_dtsub
    }
  }
  
  # Output for this step
  if (i == 1) {
    ovect= c(format((meteo$time[i]-delta_t),"%Y-%m-%d"),
      sprintf("%0.2f",initials))
    write(x=ovect, file=ofile_sim, ncolumns=length(ovect), sep="\t", append=T)
  }
  ovect= c(format(meteo$time[i],"%Y-%m-%d"), sprintf("%0.2f",states))
  write(x=ovect, file=ofile_sim, ncolumns=length(ovect), sep="\t", append=T)

}

# Read output data
sim= read.table(file= ofile_sim, sep="\t", header=T)
dbg= read.table(file= ofile_dbg, sep="\t", header=T)
dat= merge(x=sim, y=dbg, by="time", all=T)
rm(sim)
rm(dbg)

obs= read.table(file= ifile_obs, sep="\t", header=T)
names(obs)= c("time","swe_obs")
obs$swe_obs= obs$swe_obs/1000
dat= merge(x=dat, y=obs, by="time", all.x=T, all.y=F)
rm(obs)

dat$time= as.POSIXct(strptime(dat$time,"%Y-%m-%d"))

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



