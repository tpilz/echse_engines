# Input parameters
fi = 80 # mm/h
fe = 10 # mm/h
k = .5
t <- seq(1,24,by=1)# hours
wc_res <- 0.1
wc_sat <- 0.6
wc <- 0.2

inf <- function(fi,fe,k,t,wc,wc_sat,wc_res) {
  return(fe + (fi-fe) * exp(-k*t*(wc-wc_res)/(wc_sat-wc_res)))
}

plot(t, inf(fi,fe,k,t,wc,wc_sat,wc_res), type="l")
