

rel_sat= seq(0, 1, 0.01)
expo= c(0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)

f= function (rel_sat, expo, type) {
  if (type==1) return ( rel_sat^expo )
  else if (type==2) return ( 1 - rel_sat^expo )
  else if (type==3) return ( (1-rel_sat)^expo )
  else if (type==4) return ( 1 - (1-rel_sat)^expo )
  else stop("Bad type.")
}

show= function(rel_sat, type, expo) {
  plot(range(rel_sat), c(0,1), type="n", xlab="rel. sat.", ylab="term", main=paste("Type",type))
  for (i in 1:length(expo)) {
    lines(rel_sat, f(rel_sat, expo[i], type), col=i)
  }
  legend("right", lty=rep(1,length(expo)), col=1:length(expo), legend=expo, cex=0.5)
}

split.screen(c(2,2))
for (i in 1:4) {
  screen(i)
  show(rel_sat, i, expo)
}
close.screen(all=T)


