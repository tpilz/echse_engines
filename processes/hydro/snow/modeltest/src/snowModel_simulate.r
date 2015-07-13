
source("/home/dkneis/cpp/dklib/hydro_snow/r_wrapper/snow_wrapper.r")

# Simulation function
simulate = function (
  initials,
  meteo,
  params,
  n_dtsub,
  silent,
  print_sim,
  print_dbg,
  ofile_sim,
  ofile_dbg,
  replace
) {

  # Set simulation window
  tmin= meteo$time[1]
  tmax= meteo$time[nrow(meteo)]
  delta_t= as.numeric(difftime(tmax,tmin,units="secs")/(nrow(meteo)-1))
  if (!silent) print(paste("Simulation window: ",tmin," - ",tmax," (",
    as.numeric(difftime(tmax,tmin,units="days"))," days)",sep=""))
  if (!silent) print(paste("Delta t (seconds): ",delta_t,sep=""))

  # Set initials (--> order of states must match order of returned derivatives)
  states= c(
    sec= initials[["sec"]],
    swe= initials[["swe"]],
    alb= initials[["alb"]]
  )

  # Time loop
  for (i in 1:nrow(meteo)) {
    if (!silent) print(paste(meteo$time[i]," (",round(i/nrow(meteo)*100,1),"%)",sep=""))

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
    if (print_dbg) {
      dbg= snowDebug(y=states, p=p)
      if (i == 1) {
        if (file.exists(ofile_dbg) && (!replace)) stop("Output already exists.")
        ovect= c("date",names(met),names(dbg))
        write(x=ovect, file=ofile_dbg, ncolumns=length(ovect), sep="\t", append=F)
      }
      ovect= c(format(meteo$time[i],"%Y-%m-%d"),sprintf("%0.2f",met),sprintf("%0.3e",dbg))
      write(x=ovect, file=ofile_dbg, ncolumns=length(ovect), sep="\t", append=T)
    }

    # Integrate and update states
    for (k in 1:n_dtsub) {
      ddt_states= snowDerivs(t=NA, y=states, p=p)$derivs
      # Correct if SWE would become < 0
      if ((ddt_states["ddt_swe"] < 0.) &&
          (abs(ddt_states["ddt_swe"]*delta_t/n_dtsub) > states["swe"])) {
        states["swe"]= 0.00001
        states["sec"]= 0.
        states["alb"]= as.numeric(params["albedoMax"])
        break
      } else {
        states= states + ddt_states * delta_t/n_dtsub
      }
    }
    
    # Print states
    if (print_sim) {
      if (i == 1) {
        if (file.exists(ofile_sim) && (!replace)) stop("Output already exists.")
        ovect= c("date",names(initials))
        write(x=ovect, file=ofile_sim, ncolumns=length(ovect), sep="\t", append=F)
        ovect= c(format((meteo$time[i]-delta_t),"%Y-%m-%d"),
          sprintf("%0.2f",initials))
        write(x=ovect, file=ofile_sim, ncolumns=length(ovect), sep="\t", append=T)
      }
      ovect= c(format(meteo$time[i],"%Y-%m-%d"), sprintf("%0.2f",states))
      write(x=ovect, file=ofile_sim, ncolumns=length(ovect), sep="\t", append=T)
    }

  } # End time loop
} # End function


