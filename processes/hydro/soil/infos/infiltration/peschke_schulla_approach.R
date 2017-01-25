# test the approach to calculate infiltration and infiltration excess runoff
# as described in Schulla (1997) based on Green and Ampt (1911) and Peschke (1977, 1987)

# function to calculate infiltration (mm per delta_t)
infiltration <- function(
  prec,     # precipitation hitting the soil, mm (input)
  na,       # refillable porosity, m3/m3 (parameter)
  delta_t,  # time step length, s (control parameter)
  suc,      # suction at wetting front, mm (parameter)
  k_s      # saturated hydraulic conductivity, mm/s (parameter)
) {
  
  # precipitation flux, mm/s
  prec_f <- prec / delta_t
  
  # check if saturation will occur at all; write output if not
  if( k_s > prec_f ) {
    inf_vol <- prec
    runoff <- 0
    return(data.frame(infiltration = inf_vol, horton_runoff = runoff))
  }
  
  # depth of wetting front at saturation
  d_wet <- suc / (prec_f / k_s - 1)
  
  # infiltrated precipitation until saturation (mm)
  inf_sat <- d_wet * na
  
  # time until saturation (s)
  t_sat <- inf_sat / prec_f
  
  # check if saturation will occur within this time step; write output if not
  if ( t_sat > delta_t ) {
    inf_vol <- prec
    runoff <- 0
    return(data.frame(infiltration = inf_vol, horton_runoff = runoff))
  }
  
  # calculate half the amount of the remaining precipitation
  F_start <- (prec - inf_sat) / 2
  
  # calculate remaining infiltration iteratively (using F_start as start value)
  F <- F_start
  for (i in 1:12) { # not more than 12 iterations
    F_t <- k_s * (delta_t - t_sat) + na * suc * log( (F + na * suc) / (inf_sat + na * suc) )
    if(abs(F_t - F) < 0.1)
      break
    F <- F_t
  }
  
  # calculate output
  inf_vol <- inf_sat + F
  runoff <- prec - inf_vol
  return(data.frame(infiltration = inf_vol, horton_runoff = runoff))
}

# example
k <- 7.1e-4
psi <- 300
na <- 0.2
delta_t <- 3600
prec <- c(3,5,12,11,13,10,9,7,7,10,11,8,8,6,5)

infiltration(prec[1], na, delta_t, psi, k)


inf <- function(inf_i, t_sat, delta_t, k_s, na, suc, inf_sat) {
  return(inf_sat + k_s * (delta_t - t_sat) + na * suc * log( (inf_i + na * suc) / (inf_sat + na * suc) ))
}

F_t <- 22
F_res <- NULL
for(i in 1:100) {
  F_res[i] <- inf(F_t, 2000, 3600, 7e-4, 0.2, 300, 20)
#   if(abs(F_res[i]-F_t) < 0.1)
#     break
  F_t <- F_res[i]
}

