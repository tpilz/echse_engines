rm(list=ls())

# Package for ODE solving
require("deSolve")

# Parameters
nRes= 3                       # Number of reservoirs per cascade
nCas= 5                       # Number of parallel cascades
vIni= 0.                      # Initial storage
dt= 3600                      # Time step

# Outflow-storage relation
pow= function(x,e) {x^e}    # Wrapper for C++'s pow function
h_reservoir= function(v, shapeFac, shapeExp, s, v_crit) {
  if (v > v_crit) {
    return( pow(v*(shapeExp+1.)*s/shapeFac, 1./(shapeExp+1.)) );
  } else {
    return( pow(v_crit*(shapeExp+1.)*s/shapeFac, 1./(shapeExp+1.)) * v/v_crit);
  }
}
q_hole= function(h, dimOut) {
  if (h > dimOut) {
    return( 0.58 * pow(dimOut, 2.) * sqrt(2*9.81*(h-dimOut/2.)) );
  } else if (h > 0.) {
    return( 0.58 * pow(dimOut, 2.) * sqrt(2*9.81*(dimOut/2.)) * pow(h/dimOut, 1.5));
  } else {
    return( 0. );
  }
}
q_over= function(h, hCrest, wCrest) {
  return(0.65 * 2./3. * wCrest * sqrt(2.*9.81) * pow(max(0., h-hCrest), 1.5));
}
qx= function(v) {
  # Fixed parameters
  shapeFac= 1.78
  shapeExp= 2.09
  s= 0.1
  hCrest= 4
  wCrest= 8
  dimOut= 0.1
  v_crit= 0.001
  h= h_reservoir(v, shapeFac, shapeExp, s, v_crit)
  return(
    q_over(h, hCrest, wCrest)
    +
    q_hole(h, dimOut)
  )
}

# Time series of inflow
hyd= data.frame(
#  qi= c(0,0,1,2,3,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  qi= c(0,0,0.00023,0.00023,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
)

# Split inflow --> Inflow of a single cascade
hyd$qi= hyd$qi / nCas

# Function returning dv/dt for all reservoirs in a cascade (for use with 'lsoda')
derivs= function(t,y,p) {
  nres= p[["nres"]]
  dydt= rep(NA, nres)

  # Derivatives of storage volume
  dydt[1]= p[["qi"]] - qx(y[1])                       # First reservoir
  if (nres > 1) {
    for (i in 2:nres) dydt[i]= qx(y[i-1]) - qx(y[i])  # Other reservoirs
  }
  return(list(dydt))
}

# Names for output
qxend.names= paste("qxend",1:nRes,sep="")
qxavg.name= "qxavg"


# Initialize vector of state variables
v.names= paste("v",1:nRes,sep="")
y= rep(vIni, nRes)
names(y)= c(v.names)

# Instatiate vector for individual flow rates
qx_end= rep(NA, length(v.names))

# Compute outflow hydrograph of a single cascade; Results go into tempfile
tmpfl= tempfile() 
for (i in 1:nrow(hyd)) {

  # Integrate over a time step --> compute volumes
##  new= lsoda(y=y, times=c(0, dt), func=derivs, parms=c(nres=nRes, qi=hyd$qi[i]))
##  if (attr(new,"istate")[1] != 2)
##    stop("LSODA failed.")
  new= rk(y=y, times=c(0, dt), func=derivs, parms=c(nres=nRes, qi=hyd$qi[i]), method="rk45ck")

  for (k in 1:length(v.names)) qx_end[k]= qx(new[2,v.names[k]])

  # Calculate time-step-averaged outflow of last reservoir (from water balance)
  qx_avg= hyd$qi[i]  - sum(new[2,v.names] - new[1,v.names]) / dt

  # Save flow rates
  if (i == 1) {
    ovect= c("timeStep", qxend.names, qxavg.name)
    write(x=ovect, file=tmpfl, ncolumns=length(ovect), sep="\t")
  }
  ovect= c(i, signif(qx_end, 4), signif(qx_avg, 4))
  write(x=ovect, file=tmpfl, ncolumns=length(ovect), sep="\t", append=T)

  # Set initial volumes for next time step
  y[v.names]= new[2,v.names]

}

# Read simulated flow rates from tempfile
res= read.table(file=tmpfl, header=TRUE, sep="\t")

# Upscale for the number of parallel cascades
res[,qxend.names]= res[,qxend.names] * nCas
res[,qxavg.name]= res[,qxavg.name] * nCas
hyd$qi= hyd$qi * nCas

# Merge tables of in- and outputs
if (nrow(res) != nrow(hyd))
  stop("Mismatch in length of tables.")
hyd= cbind(res,hyd)
rm(res)

# Plot hydrographs
# Inflow
plot(hyd$timeStep,hyd$qi, type="S", ylim=c(0,max(hyd$qi)*1.25), col="black", xlab="end of time step", ylab="q")
# Instantaneous outflow rates of individual reservoirs
clr= rainbow(length(qxend.names))
for (i in 1:length(qxend.names)) {
  lines(hyd$timeStep,hyd[,qxend.names[i]], type="S", lty=3, col=clr[i])
}
# Time-step averaged outflow of last reservoir
lines(hyd$timeStep,hyd[,qxavg.name], type="S", lty=2, col=clr[length(clr)])
# Legend
legend("topright",bty="n", cex=0.8,
  lty=c(1, rep(3,length(qxend.names)),2), col=c(rgb(0,0,0),clr,clr[length(clr)]),
  legend=c("qi",qxend.names,qxavg.name))

# Check water balance
vin= sum(hyd$qi)*dt
vex= sum(hyd[qxavg.name])*dt
vst= sum(y[v.names]) * nCas
print(paste("Cumulated inflow: ",vin))
print(paste("Cumulated outflow:",vex))
print(paste("Storage:          ",vst))
print(paste("Balance (in-out): ",vin - vex - vst))




