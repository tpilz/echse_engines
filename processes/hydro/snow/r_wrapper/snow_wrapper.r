
shared_lib= "/home/dkneis/cpp/dklib/hydro_snow/r_wrapper/snow_wrapper.so"

# Load shared library (unload first as it may have changed)
if (is.loaded("snowModel_derivs") || is.loaded("snowModel_debug")) dyn.unload(shared_lib)
dyn.load(shared_lib)

# Wrapper for the derivatives function
if (!is.loaded("snowModel_derivs")) stop("Symbol 'snowModel_derivs' not available.")
snowDerivs= function(
  t,  # time (unused)
  y,  # states
  p   # parameters (includes forcings which are taken as const over a time step)
) {
  res= .C("snowModel_derivs",
    # Input part 1: Forcings (passed as parameters for a time step)
    as.double(p["precipSumMM"]),
    as.double(p["shortRad"]),
    as.double(p["tempAir"]),
    as.double(p["pressAir"]),
    as.double(p["relHumid"]),
    as.double(p["windSpeed"]),
    as.double(p["cloudCoverage"]),
    # Input part 2: True parameters
    as.double(p["precipSeconds"]),
    as.double(p["a0"]),
    as.double(p["a1"]),
    as.double(p["kSatSnow"]),
    as.double(p["densDrySnow"]),
    as.double(p["specCapRet"]),
    as.double(p["emissivitySnowMin"]),
    as.double(p["emissivitySnowMax"]),
    as.double(p["tempAir_crit"]),
    as.double(p["albedoMin"]),
    as.double(p["albedoMax"]), 
    as.double(p["agingRate_tAirPos"]),
    as.double(p["agingRate_tAirNeg"]),
    as.double(p["soilDepth"]),
    as.double(p["soilDens"]),
    as.double(p["soilSpecHeat"]),
    as.double(p["weightAirTemp"]),
    # Input part 3: States
    as.double(y["sec"]), # snowEnergyCont
    as.double(y["swe"]), # snowWaterEquiv
    as.double(y["alb"]), # albedo
    # Outputs
    change_sec= double(1),
    change_swe= double(1),
    change_alb= double(1),
    flux_melt= double(1),
    flux_subl= double(1)
  )
  return(list(
   derivs= c(ddt_sec=res$change_sec, ddt_swe=res$change_swe, ddt_alb=res$change_alb),
   others= c(flux_melt=res$flux_melt, flux_subl=res$flux_subl)
  ))
}

# Wrapper for the debug function
if (!is.loaded("snowModel_debug")) stop("Symbol 'snowModel_debug' not available.")
snowDebug= function(
  y,  # states
  p   # parameters (includes forcings which are taken as const over a time step)
) {
  res= .C("snowModel_debug",
    # Input part 1: Forcings (passed as parameters for a time step)
    as.double(p["precipSumMM"]),
    as.double(p["shortRad"]),
    as.double(p["tempAir"]),
    as.double(p["pressAir"]),
    as.double(p["relHumid"]),
    as.double(p["windSpeed"]),
    as.double(p["cloudCoverage"]),
    # Input part 2: True parameters
    as.double(p["precipSeconds"]),
    as.double(p["a0"]),
    as.double(p["a1"]),
    as.double(p["kSatSnow"]),
    as.double(p["densDrySnow"]),
    as.double(p["specCapRet"]),
    as.double(p["emissivitySnowMin"]),
    as.double(p["emissivitySnowMax"]),
    as.double(p["tempAir_crit"]),
    as.double(p["albedoMin"]),
    as.double(p["albedoMax"]), 
    as.double(p["agingRate_tAirPos"]),
    as.double(p["agingRate_tAirNeg"]),
    as.double(p["soilDepth"]),
    as.double(p["soilDens"]),
    as.double(p["soilSpecHeat"]),
    as.double(p["weightAirTemp"]),
    # Input part 3: States
    as.double(y["sec"]), # snowEnergyCont
    as.double(y["swe"]), # snowWaterEquiv
    as.double(y["alb"]), # albedo
    # Outputs
    TEMP_MEAN= double(1),
    TEMP_SURF= double(1),
    LIQU_FRAC= double(1),
    flux_R_netS= double(1),
    flux_R_netL= double(1),
    flux_R_soil= double(1),
    flux_R_sens= double(1),
    stoi_f_prec= double(1),
    stoi_f_subl= double(1),
    stoi_f_flow= double(1),
    flux_M_prec= double(1),
    flux_M_subl= double(1),
    flux_M_flow= double(1),
    rate_G_alb= double(1)
  )
  return( c(
    temp_mean= res$TEMP_MEAN,
    temp_surf= res$TEMP_SURF,
    liqu_frac= res$LIQU_FRAC,
    flux_R_netS= res$flux_R_netS,
    flux_R_netL= res$flux_R_netL,
    flux_R_soil= res$flux_R_soil,
    flux_R_sens= res$flux_R_sens,
    stoi_f_prec= res$stoi_f_prec,
    stoi_f_subl= res$stoi_f_subl,
    stoi_f_flow= res$stoi_f_flow,
    flux_M_prec= res$flux_M_prec,
    flux_M_subl= res$flux_M_subl,
    flux_M_flow= res$flux_M_flow,
    rate_G_alb= res$rate_G_alb
  ))
}


