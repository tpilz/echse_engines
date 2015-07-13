
# Response function of a single linear reservoir to constant inflow
qout = function (k,dt,qi,v0) {
  qi - qi * exp(-dt/k) + v0/k * exp(-dt/k)
}

qi= c(0,0,1,1,0,3,4,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,)

dt= 1
k=5
qx_ini= 0

qx= rep(NA,length(qi))
qxm= rep(NA,length(qi))
for (i in 1:length(qi)) {
  if (i == 1) {
    v0= qx_ini*k
  } else {
    v0= qx[i-1]*k
  }
  qx[i]= qout(k,dt,qi[i],v0)
  qxm[i]= (v0 + qi[i] * dt - qx[i]*k) / dt
}

plot(1:length(qi),qi,type="S",lwd=1)
lines(1:length(qi),qx,type="S",lwd=2)
lines(1:length(qi),qxm,type="S",col="blue")

print(sum(qi))
print(sum(qx))
print(sum(qxm))

