rm(list=ls())

# BEGIN SETTINGS ###############################################################

#ifile= paste(Sys.getenv("path_progress"),"/echse_core/processes/snow/modeltest/data/DWD/use/Fichtelberg/sunhrs.txt",sep="")
#ofile= paste(Sys.getenv("path_progress"),"/echse_core/processes/snow/modeltest/data/DWD/use/Fichtelberg/glorad.txt",sep="")
ifile= paste(Sys.getenv("path_progress"),"/echse_core/processes/snow/modeltest/data/DWD/use/KahlerAsten/sunhrs.txt",sep="")
ofile= paste(Sys.getenv("path_progress"),"/echse_core/processes/snow/modeltest/data/DWD/use/KahlerAsten/glorad.txt",sep="")

latitude_degree= 50.429

# END SETTINGS #################################################################

# Load conversion function
library("meteovars")

# Read data
dat= read.table(file=ifile, sep="\t", header=T, colClasses=c("POSIXct","numeric"))
names(dat)= c("time","sunhours")

# Convert
dat$daynum= as.integer(format(dat$time,"%j"))
dat$value= NA
for (i in 1:nrow(dat)) {
  dat$value[i]= glorad_from_sunhrs(dat$sunhours[i], dat$daynum[i], latitude_degree)
}

# Delete unnecessary columns
dat$daynum= NULL
dat$sunhours= NULL

# Show result
plot(dat$time, dat$value, pch=20)

# Create output
dat$value= round(dat$value, 1)
write.table(x=dat, file= ofile, sep="\t", col.names=T, row.names=F, quote=F)

