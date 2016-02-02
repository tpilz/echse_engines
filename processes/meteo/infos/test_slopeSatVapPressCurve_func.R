# magnus formula to compute saturation water vapor pressure
magnus_dp_exp <- function(temp) {
  return(6.11 * exp(17.3 * temp / (237.3 + temp)))
}

# latent heat for condensation of water, https://en.wikipedia.org/wiki/Latent_heat
Lv_calc <- function(temp) {
  return( (2500.8 - 2.36*temp + 0.0016*temp^2 - 0.00006*temp^3) * 1000)
}

# solpe of saturation vapor pressure - temperature curve, i.e. Clausius–Clapeyron equation
# derived from thermodynamics for standard atmospheric conditions: https://en.wikipedia.org/wiki/Clausius–Clapeyron_relation
slopeSatVapCurve_cce <- function(es, Rv, Lv, temp) {
  return(Lv * es / (Rv*temp^2))
}

# derived from integration of magnus formula (Dyck & Peschke, 1995, eq. 11.14)
slopeSatVapCurve_dp <- function(es, temp) {
  return(4098 * es /(237.3+temp)^2)
}

# comparison
temp <- seq(-20,30, by=1)

Rv <- 461.5 # gas constant of water vapor, https://en.wikipedia.org/wiki/Water_vapor, [Jkg-1K-1]
Lv <- 2264.76e3 # specific latent heat of evaporation of water, [Jkg-1]
plot(temp, slopeSatVapCurve_cce(magnus_dp_exp(temp), Rv, Lv_calc(temp), temp+273.15), type="l", xlab="Temperature[°C]", ylab="Slope of saturation vapor pressure curve [hPaK-1]")
lines(temp, slopeSatVapCurve_dp(magnus_dp_exp(temp), temp), col="red")
