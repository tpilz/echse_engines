### INPUT ###
# current hour of day
hour <- seq(0,23)
# current day of year
doy <- 29
# latitude of location (deg)
lat <- 52.5 # Berlin
# Longitude of the location of interest (decimal degrees west of Greenwich)
L_m = 346.5 # Berlin
# deviation of local time zone from UTC (may vary over the year due to daylight saving time)
utc_add = 1


### FUNCTIONS ###
# function to compute sunset/sunrise hour angle
dayTime_fac <- function(doy,lat,delta) {
  # latitude in (rad)
  phi = abs(lat * pi / 180.)
  
  sunrise = -tan(delta) * tan(phi)
  
  if (sunrise > 1.)
    sunrise = 1.
  
  if (sunrise < -1.)
    sunrise = -1
  
  return(acos(sunrise))
}

# function to compute hourly raterrestrial radiation
rad_extraterr_hourly <- function(doy, lat, hour, utc_add, L_m) {

  # solar constant
  SOLAR_C = 1360.8
    
  # compute centre of current time zone
  L_z = (360 - utc_add * 15) %% 360
  
  # solar declination (rad)
  delta = asin(0.4 * sin(2. * pi / 365. * (doy - 82.)))
  
  # eccentricity correction factor (-)
  E_0 = 1. + 0.033 * cos(2. * pi * doy / 365.)
  
  # latitude in (rad)
  phi = lat * pi / 180.
  
  # seasonal correction factor, eqs. 32, 33
  b = 2. * pi * (doy - 81.) / 364.
  S_c = 0.1645 * sin(2.*b) - 0.1255 * cos(b) - 0.025 * sin(b)
  
  # solar time angle at the midpoint of the hour of day (rad), eq. 31
  hour_mid = hour + 0.5
  w_mid = pi / 12. * ( (hour_mid + 0.06667 * (L_z-L_m) + S_c) - 12. )
  
  # solar time angles and begin and end of hour, eqs. 29, 30
  w_ini = w_mid - pi / 24.
  w_end = w_mid + pi / 24.
  
  # sunset/sunrise hour angle (rad)
  w_s = dayTime_fac(doy,lat,delta)
  
  # hours without sunshine, i.e. incoming radiation (+/- half an hour accuracy)
  if ( (w_mid < -w_s) | (w_mid > w_s) )
    return(0)
  
  # calculate radiation
  return(1./pi * SOLAR_C * E_0 * ( (w_end-w_ini) * sin(delta) * sin(phi) + cos(delta) * cos(phi) * (sin(w_end)-sin(w_ini)) ))

}

### TEST: daily evolution of radex for Berlin at different doy ###
# at spring equinox
rad_vals_spring <- sapply(hour, function(x) rad_extraterr_hourly(79, lat, x, 1, L_m))
# at summer solstice
rad_vals_summer <- sapply(hour, function(x) rad_extraterr_hourly(172, lat, x, 2, L_m))
# at autumn equinox
rad_vals_autumn <- sapply(hour, function(x) rad_extraterr_hourly(266, lat, x, 2, L_m))
# at winter solstice
rad_vals_winter <- sapply(hour, function(x) rad_extraterr_hourly(356, lat, x, 1, L_m))

pdf("test_radex_hourly.pdf")

par(cex=1.25, mar=c(4.5,4.5,1,1), cex.lab=1.5, lwd=2)
plot(hour, rad_vals_summer, type="l", col="red", ylab="Extraterrestrial radiation [Wm-2]", xlab="Hour of day", xaxt="n", las=1)
axis(1, at=hour)
lines(hour, rad_vals_spring, col="green")
lines(hour, rad_vals_autumn, col="brown")
lines(hour, rad_vals_winter, col="blue")
legend("topleft", legend=c("Summer solstice", "Spring equinox", "Autumn equinox", "Winter solstice"),
       col=c("red", "green", "brown", "blue"), lty=1)

dev.off()
