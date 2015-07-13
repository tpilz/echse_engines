rm(list=ls())

# This script creates ready-to-use meteo data from DWD raw data. These raw data
# are subject to missing records (for some dates) and also duplicate records
# (multiple records for the same date).

# BEGIN SETTINGS ###############################################################

#idir= paste(Sys.getenv("path_progress"),"/echse_core/processes/snow/modeltest/data/DWD/raw/Fichtelberg",sep="")
#odir= paste(Sys.getenv("path_progress"),"/echse_core/processes/snow/modeltest/data/DWD/use/Fichtelberg",sep="")
idir= paste(Sys.getenv("path_progress"),"/echse_core/processes/snow/modeltest/data/DWD/raw/KahlerAsten",sep="")
odir= paste(Sys.getenv("path_progress"),"/echse_core/processes/snow/modeltest/data/DWD/use/KahlerAsten",sep="")

vect_vars= c("apress", "clness", "precip", "rhumid", "sunhrs", "temper", "windsp")

time_first= ISOdatetime(2002,1,1,0,0,0)
time_last= ISOdatetime(2011,6,30,0,0,0)

# END SETTINGS #################################################################

# Make continous time axis
times= data.frame(
  time= seq(from= time_first, to=time_last, by=86400)
)
print(paste("Anzahl Zeiten:",nrow(times)))

# Loop through variables
for (var in vect_vars) {

  # Read
  print(paste("Reading '",var,"'",sep=""))
  dat= read.table(file= paste(idir,"/",var,".txt",sep=""), header=T, sep="\t",
    colClasses= c("POSIXct","numeric"))
  names(dat)= c("time","value")

  # Remove duplicate records (keep first record for a particular date)
  x= tapply(X=dat$time, INDEX=dat$time, FUN=length)
  x= which(x > 1)
  print(paste(" -- Duplicate dates: ",length(x),sep=""))
  if (length(x) > 0) {
    for (i in 1:length(x)) {
      rowinds_dup= which(as.character(dat$time) == names(x[i]))
      if (length(rowinds_dup) < 2) stop("Bug.")
      rowinds_del= rowinds_dup[2:length(rowinds_dup)]
      dat$value[rowinds_del]= NA
    }
    print(paste(" -- Records removed: ",sum(is.na(dat$value)),sep=""))  
    dat= dat[!is.na(dat$value),] 
  }

  # Make equidistant time axis
  nmiss= nrow(times)-nrow(dat)
  print(paste(" -- Missing dates:",nmiss))
  dat= merge(x=times, y=dat, by.x="time", by.y="time", all.x=T)
  if (sum(is.na(dat$value)) != nmiss) stop("Bug.")

  # Fill missing data
  if (nmiss > 0) {
    tmp= approx(x=dat$time, y=dat$value, xout=dat$time[which(is.na(dat$value))],
      method="linear", rule=1, ties=mean)
    for (i in 1:length(tmp$x)) {
      rowind= which(as.character(dat$time) == tmp$x[i])
      if (length(rowind) != 1) stop("Bug.")
      dat$value[rowind]= tmp$y[i]
    }
    print(paste(" -- All missing data filled?: ",
      as.logical(sum(is.na(dat$value)) == 0),sep=""))
  }

  # Convert cloudiness (from range 0...8 to range 0...1)
  if (var == "clness") dat$value= dat$value/8.

  # Create output
  print("Creating output...")
  dat$value= round(dat$value, 1)
  write.table(x=dat, file= paste(odir,"/",var,".txt",sep=""), sep="\t",
    col.names=T, row.names=F, quote=F)

}


