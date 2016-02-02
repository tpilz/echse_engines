# magnus formula in exponential form by Dyck & Peschke (1995), eq. 11.13
magnus_dp_exp <- function(temp) {
  return(6.11 * exp(17.3 * temp / (237.3 + temp)))
}

# magnus formula as power to 10 by Dyck & Peschke (1995), eq. 4.6
magnus_dp_pow10 <- function(temp) {
  return( 6.11 * 10.^(7.5*temp/(237.3+temp)) )
}

# magnus formula over ice by Dyck & Peschke (1995), eq. 4.6a
magnus_dp_ice <- function(temp) {
  return( 6.11 * 10.^(9.5*temp/(265.5+temp)) )
}

# magnus formula from https://en.wikipedia.org/wiki/Clausius–Clapeyron_relation (similar to magnus_dp_exp assuming different standard temperature?!)
magnus_wiki <- function(temp) {
  return(6.1094 * exp(17.625 * temp / (243.04 + temp)))
}

# magnus formula as in SWAT (SWAT manual (Neitsch, 2011) eq. 1:2.3.2)
magnus_swat <- function(temp) {
  return( exp( (16.87 * temp - 116.9) / (temp + 237.3) ) * 10)
}

# comparison
temp <- seq(-20,30, by=1)

pdf("compare_satVapPress_formulas.pdf")

par(cex=1.25, mar=c(4.5,4.5,1,1), cex.lab=1.5, lwd=2)
plot(temp, magnus_dp_exp(temp), type="l", xlab="Temperature [°C]", ylab="Saturation vapor pressure [hPa]", las=1)
lines(temp, magnus_dp_pow10(temp), col="black", lty=2)
lines(temp, magnus_dp_ice(temp), col="purple")
lines(temp, magnus_wiki(temp), col="red")
lines(temp, magnus_swat(temp), col="green")
legend("topleft", legend=c("Dyck & Peschke (1995) eq. 11.13", "Dyck & Peschke (1995) eq. 4.6", "Dyck & Peschke (1995) eq. 4.6a (ice)", "Wikipedia (en)", "SWAT 2011"),
       lty=c(1,2,1,1,1), col=c("black", "black", "purple", "red", "green"))

dev.off()
