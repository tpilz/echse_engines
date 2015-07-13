rm(list=ls())

# Load functions
source("/home/dkneis/progress/echse_apps/hypsoRR/processes/snow/modeltest/src/snowModel_readMeteo.r")
source("/home/dkneis/progress/echse_apps/hypsoRR/processes/snow/modeltest/src/snowModel_simulate.r")
source("/home/dkneis/progress/echse_apps/hypsoRR/processes/snow/modeltest/src/snowModel_gof.r")
source("/home/dkneis/progress/echse_apps/hypsoRR/processes/snow/modeltest/src/snowModel_show.r")

source("/home/dkneis/progress/echse_tools/precipcorr/precipcorr.r")

# Get settings
source("/home/dkneis/progress/echse_apps/hypsoRR/processes/snow/modeltest/src/snowModel_settings.r")

# Read meteo data
meteo= readMeteo(
  dir_meteo,
  station
)
names(meteo)[which(names(meteo) == "precip")] = "precip_raw"

# Read sensor properties
sens= read.table(file=file_sensors, header=T, sep="\t")
sens= sens[as.character(sens$station)==station, which(names(sens) != "station")]
sens= unlist(sens)

# Convert windspeed to values as standard height (10m)
meteo$windsp= meteo$windsp * log(10/sens["height_zeroWind"])/
  log(sens["sensorHeight_windsp"]/sens["height_zeroWind"])

# Read initial data
initials= read.table(file=file_initials, header=T, sep="\t")
initials= initials[as.character(initials$station)==station, which(names(initials) != "station")]
initials= unlist(initials)

# Read observations
obs= read.table(file=file_obs, header=T, sep="\t")
obs[,col_time]= as.POSIXct(strptime(obs[,col_time],"%Y-%m-%d"))

# Read time intervals for GOF analysis and set observations filter
timesGOF= read.table(file=file_timesGOF, sep="\t", header=TRUE,
  colClasses= c("character","POSIXct","POSIXct"))
timesGOF= timesGOF[timesGOF$station == station,]
gof_selected= rep(FALSE, nrow(obs))
for (i in 1:nrow(timesGOF)) {
  gof_selected[(obs[,col_time] >= timesGOF$from[i]) & (obs[,col_time] <= timesGOF$till[i])]= TRUE
}

# Read parameters
tab_params= read.table(file=file_params, header=T, sep="\t")
if (any(tab_params$min > tab_params$max)) stop("Error in table of parameters.")

# MONT-CARLO MODE ##############################################################
if (mode == "mcs") {
  require("lhs")
  print("Drawing LHS sample...")
  lhSample= improvedLHS(n= mcs_nSamples, k=nrow(tab_params))
  print("Entering MC loop...")
  for (i in 1:mcs_nSamples) {
    print(paste("This is run",i,"of",mcs_nSamples))
    print("  Setting parameters...")    
    p= qunif(p= lhSample[i,], min= tab_params$min, max= tab_params$max)
    names(p)= tab_params[,"parameter"]
    print("  Correcting precipitation...")    
    meteo$precip= meteo$precip_raw * factPrecipCorr_richter(meteo$precip_raw,
      meteo$temper, tempCrit_lower= p["tempAir_crit"]-0.5,
      tempCrit_upper= p["tempAir_crit"]+0.5, protection_index=sens["gage_protectionIndex"])
    print("  Running simulation...")
    simulate(initials, meteo, params=p,
      n_dtsub, silent=TRUE, print_sim=TRUE, print_dbg=FALSE, ofile_sim, ofile_dbg, replace)
    print("  Reading result...")
    sim= read.table(file=ofile_sim, header=T, sep="\t")
    sim[,col_time]= as.POSIXct(strptime(sim[,col_time],"%Y-%m-%d"))
    print("  Computing GOF...")
    gof= computeGOF(obs=obs, sim=sim, col_sim="swe", col_obs=col_obs,
      col_time=col_time, fact_obs=fact_obs, selected_obs=gof_selected)
    print("  Storing GOF and corresponding parameters...")    
    if (i == 1) {
      if (file.exists(ofile_mcs_runs) & (!replace)) stop("Output exists.")
      ovect= c("gof",names(p))
      write(x=ovect, file=ofile_mcs_runs, ncolumns=length(ovect), sep="\t", append=F)
    }
    ovect= c(gof,p)    
    write(x=ovect, file=ofile_mcs_runs, ncolumns=length(ovect), sep="\t", append=T)
  }
  print("MC loop done.")
  print("Analyzing result...")
  res= read.table(file=ofile_mcs_runs, sep="\t", header=T)
  pdf(file=ofile_pdf, paper="a4r", width=25/2.54, height=17/2.54)
  nc= ceiling(length(p)/3)
  layout(matrix(1:(nc*3),ncol=nc,byrow=T),width=1,height=1)
  for (i in 1:length(p)) {
    plot(res[,names(p)[i]],res$gof,pch=20,xlab=names(p)[i],ylab="GOF",main="")
  }
  graphics.off()
  print("Print best set...")
  if (file.exists(ofile_mcs_best) & (!replace)) stop("Output exists.")
  write.table(
    x= data.frame(parameter=names(p), value=unlist(res[which.min(res$gof),names(p)])),
    file=ofile_mcs_best, sep="\t", col.names=T, row.names=F, quote=F)

# STANDARD MODE ################################################################
} else if (mode == "std") {
  print("Setting parameters...")
  p= as.numeric(tab_params[,"value"])
  names(p)= tab_params[,"parameter"]
  print("Correcting precipitation...")    
  meteo$precip= meteo$precip_raw * factPrecipCorr_richter(meteo$precip_raw,
    meteo$temper, tempCrit_lower= p["tempAir_crit"]-0.5,
    tempCrit_upper= p["tempAir_crit"]+0.5, protection_index=sens["gage_protectionIndex"])
  print("Running simulation...")
  simulate(initials, meteo, params=p,
    n_dtsub, silent=TRUE, print_sim=TRUE, print_dbg=FALSE, ofile_sim, ofile_dbg, replace)
  print("Reading result...")
  sim= read.table(file=ofile_sim, header=T, sep="\t")
  sim[,col_time]= as.POSIXct(strptime(sim[,col_time],"%Y-%m-%d"))
  print("Plotting...")
  pdf(file=ofile_pdf, paper="a4r", width=25/2.54, height=17/2.54)
#  postscript(file=ofile_pdf, paper="a4", width=25/2.54, height=17/2.54)
  show(obs=obs, sim=sim, col_sim="swe", col_obs=col_obs,
    col_time="date", fact_obs=fact_obs, selected_obs=gof_selected)
  graphics.off()
  print("Computing GOF...")
  gof= computeGOF(obs=obs, sim=sim, col_sim="swe", col_obs=col_obs,
    col_time=col_time, fact_obs=fact_obs, selected_obs=gof_selected)
  print(gof)
} else {
  stop("Unknown mode.")
}

