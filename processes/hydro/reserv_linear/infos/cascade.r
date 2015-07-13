
# dt-Impulsantwort der linearen Speicherkaskade
# Maniak S. 264, Gl. 6.35
# auch Dyck 1989, Bd. 1, S. 105, Gl. 4.37
dt_response= function (n,k,dt,m) {
  return(dt/prod(1:(n-1))/k * (m*dt/k)^(n-1) * exp(-m*dt/k))
}

dt=60         # dt in Sekunden
n= 40         # Anzahl Speicher              ---> typ. Wert aus Profilen
k= 20         # k in Sekunden                ---> typ. Wert aus Profilen
m= 0:round(3600/dt) # Index des Zeitintervalls
dat= data.frame(m= m*dt, r= dt_response(n,k,dt,m))
#dat$r= dat$r * 1/sum(dat$r)

# Man sieht an diesem Bsp., dass dt in einem sinnvollen Verhältnis zu n und k
# stehen muss. Bei den gegebenen Werten von n und k erhält man für 

plot(dat,type="S",xlab="t (s)", ylab= paste("h(dt=",dt," ,t)",sep=""),
  main=paste("Summe:",sum(dat$r)))


# Outflow of a single linear reservoir at time 0+dt with inflow varying
# linearily over the time step of length dt.
# This solution is, for example, used by
#    (1) Maniak (2005), p. 261
#    (2) Bremicker (2000), p. 39
#    (3) Pfützner (hydroskript - KalMil)
# Qi1: Inflow at t=0
# Qi2: Inflow at t=dt
# Qx1: Outflow at t=0
# k:   Retention constant
# dt:  Length of time step
Qx2_single= function (Qi1,Qi2,Qx1,k,dt) {
  c1= 1 - exp(-k/dt)
  c2= 1 - c1 * k/dt
  return(
    Qx1 + (Qi1-Qx1) * c1 + (Qi2-Qi1) * c2
  )
}

nres=2
dt= 60
ndtsub= 1
k= 20

q_in= c(0,0,1,0,0)

dtsub= dt/ndtsub

V= rep(0,nres)
V_corr= rep(NA,nres)

q_out= rep(NA,length(q_in))

# Primary time loop
for (idt in 1:length(q_in)) {

  v_out= 0
  # Secondary time loop
  for (idtsub in 1:ndtsub) {
  
    # Reservoir loop
    Qx2= rep(NA,nres)
    for (ires in 1:nres) {
      if (ires == 1) {
        Qi1= q_in[idt]
        Qi2= q_in[idt]
      } else {
        Qi1= 1/k*V[ires-1]
        Qi2= Qx2[ires-1]
      }
      Qx1= 1/k*V[ires]
      Qx2[ires]= Qx2_single(Qi1,Qi2,Qx1,k,dtsub)

print(paste("idt",idt,"ires",ires,"v",V[ires],"Qi1",Qi1,"Qi2",Qi2,"Qx1",Qx1,"Qx2",Qx2[ires]))

      V_corr[ires]= V[ires] + (Qi1+Qi2)/2*dtsub - (Qx1+Qx2[ires])/2*dtsub
    }
    v_out= v_out + (Qx1+Qx2[nres])/2*dtsub # dont move because of Qx1
    V= V_corr

print(paste("after idt",idt,"v=",V))

  }
  q_out[idt]= v_out / dt
}

plot(1:length(q_in),q_in,type="S")
lines(1:length(q_out),q_out,type="S",lwd=2)


