
# INPUT
ksat = 7.74777048045525E-006 # m/s
bubble = 0.0980958301 # m
pores_ind = 0.3617396378 # (-)
wc_sat = 0.413232826 # (-)
wc_res = 0.0744021222 # (-)
wc_pwp = 0.10477 # (-)
wc_fc = 0.202178076923077 # (-)
wc = seq(wc_res,wc_sat, by=0.01)




# MATRIC POTENTIAL
suc_brooks <- function(wc, wc_sat, wc_res, pores_ind, bubble) {
  sat_rel = (wc - wc_res) / (wc_sat - wc_res)
  sat_rel = max(sat_rel, 0.001)
  suction = bubble / sat_rel^{1./pores_ind}
  
  return(suction)
}

suc_gen<- function(wc, wc_sat, wc_res, pores_ind, bubble) {
  sat_rel = (wc - wc_res) / (wc_sat - wc_res)
  sat_rel = max(sat_rel, 0.001)
  m = pores_ind / (pores_ind + 1.)
  suction = (1./sat_rel^{1./m} - 1.)^{1./(pores_ind+1.)} * bubble
  
  return(suction)
}

suc_camp <- function(wc, wc_sat, pores_ind, bubble) {
  sat_rel = wc / wc_sat
  suction = bubble / sat_rel^{1./pores_ind}
  
  return(suction)
}


# CONDUCTIVITY
cond_brooks <- function(wc, wc_sat, wc_res, pores_ind, ksat) {
  # relative saturation (-)
  sat_rel = (wc - wc_res) / (wc_sat - wc_res)
  
  # Brooy and Corey parameter n from pore size index (-)
  n = 3. + 2. / pores_ind
  
  # unsaturated conductivity, m/s
  k_u = ksat * sat_rel^n
  
  return(k_u)
}

cond_gen <- function(wc, wc_sat, wc_res, pores_ind, ksat) {
  # relative saturation (-)
  sat_rel = (wc - wc_res) / (wc_sat - wc_res)
  
  # van-Genuchten parameter m from pore size index (-)
  m = pores_ind / (pores_ind + 1.);	
  
  # unsaturated conductivity (according to van Genuchten), m/s
  k_u = ksat * sat_rel^0.5 * (1. - (1. - sat_rel^{1./m})^m)^2;
  
  return(k_u)
}

cond_camp <- function(wc, wc_sat, pores_ind, ksat) {
  # relative saturation (-)
  sat_rel = wc / wc_sat
  
  # Campbell parameter n from pore size index (-)
  n = 3. + 2. / pores_ind
  
  # unsaturated conductivity, m/s
  k_u = ksat * sat_rel^n
  
  return(k_u)
}





# COMPARISON SUCTION
vals_corey <- sapply(wc, function(x) suc_brooks(x, wc_sat, wc_res, pores_ind, bubble))
vals_gen <- sapply(wc, function(x) suc_gen(x, wc_sat, wc_res, pores_ind, bubble))
vals_camp <- sapply(wc, function(x) suc_camp(x, wc_sat, pores_ind, bubble))

par(mar=c(4.5,4.5,0,0))

plot(wc, log10(vals_corey), yaxt="n", type="l", ylab="", xlab="Soil water content [m3/m3]",
     col="blue", las=1, ylim=c(log10(0.01), log10(1e4)), lwd=2)

ats = log10(c(0.1, 1, 10, 100, 1000, 10000))
labs = 10^ats
axis(2, at=ats, labels=labs, las=1)
title(ylab="Suction [m] / [100 hPa]", line=3.5)
abline(v=wc_sat, lty=2, col="blue")
abline(v=wc_fc, lty=2, col="brown")
abline(v=wc_pwp, lty=2, col="green")
abline(v=wc_res, lty=2, col="red")

lines(wc, log10(vals_gen), col="red", lwd=2)
lines(wc, log10(vals_camp), col="green", lwd=2)

legend("topright", legend = c("Van Genuchten", "Brooks and Corey", "Campbell"), col=c("red", "blue", "green"), lty=1, lwd=2)


# COMPARISON CONDUCTIVITY
vals_cond_corey <- sapply(wc, function(x) cond_brooks(x, wc_sat, wc_res, pores_ind, ksat)) * 3600 * 1000
vals_cond_gen <- sapply(wc, function(x) cond_gen(x, wc_sat, wc_res, pores_ind, ksat)) * 3600 * 1000
vals_cond_camp <- sapply(wc, function(x) cond_camp(x, wc_sat, pores_ind, ksat)) * 3600 * 1000

par(mar=c(4.5,4.5,0,0))

plot(wc, log10(vals_cond_corey), yaxt="n", type="l", ylab="", xlab="Soil water content [m3/m3]",
     col="blue", las=1, ylim=c(log10(0.0001), log10(1e2)), lwd=2)

ats = log10(c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100))
labs = 10^ats
axis(2, at=ats, labels=labs, las=1)
title(ylab="Conductivity [mm/h]", line=3.5)

lines(wc, log10(vals_cond_gen), col="red", lwd=2)
lines(wc, log10(vals_cond_camp), col="green", lwd=2)

legend("topleft", legend = c("Van Genuchten", "Brooks and Corey", "Campbell"), col=c("red", "blue", "green"), lty=1, lwd=2)

