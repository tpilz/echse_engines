# function implemented by David citing: http://www.fao.org/docrep/X0490E/x0490e07.htm#atmospheric pressure
apress_simple <- function(elev) {
  return( 1013.25 * ( 1 - (0.0065 * elev / 293.))^(5.255))
}

# by Tobias citing: SWAT manual (Neitsch, 2011) eq. 1:2.3.8 citing Doorenbos and Pruitt (1977) (cannot be found in here?!)
air_press <- function(elev) {
  return( (101.3 - 0.01152 * elev + 0.544e-06 * elev^2) * 10 )
}

# compare
elev <- seq(0,3000, by=1)
plot(elev, apress_simple(elev), type="l", ylab="Atmospheric pressure [hPa]", xlab="Elevation a.s.l. [m]")
lines(elev, air_press(elev), col="red")
legend("topright", legend=c("apress_simple", "air_press"), col=c("black", "red"), lty=1)

plot(elev, air_press(elev)-apress_simple(elev), type="l", ylab="Absolute deviation [hPa]", xlab="Elevation a.s.l. [m]")

# CONCLUSION:
# Absolute deviations relatively low (< 2.67 hPa until 1000m, < 10.8 hPa until 3000m)
