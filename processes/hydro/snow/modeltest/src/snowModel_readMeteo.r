
# Reads all meteo data into a data frame
readMeteo = function (
  dir_meteo,
  station
) {

  meteovars= c("precip", "glorad", "temper", "apress", "rhumid", "windsp", "clness")
  ext_meteo= ".txt"

  # Read forcings and merge into a single table
  for (i in 1:length(meteovars)) {
    ifile= paste(dir_meteo,"/",meteovars[i],ext_meteo,sep="")
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

  return(meteo)
}

