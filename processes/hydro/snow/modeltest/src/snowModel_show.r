
show = function (obs, sim, col_sim, col_obs, col_time, fact_obs, selected_obs) {
  # Plot swe for winters
  op= par(no.readonly=T)
  par(mar=c(4,2.5,1.2,0.5))
  years= as.integer(unique(format(sim[,col_time],"%Y")))
  nc= ceiling(length(years)/2)
  layout(matrix(1:(nc*2),ncol=nc,byrow=T),width=1,height=1)
  mi= min(c(sim[,col_sim],obs[,col_obs]*fact_obs),na.rm=T)
  mx= max(c(sim[,col_sim],obs[,col_obs]*fact_obs),na.rm=T)
  for (i in 1:length(years)) {
    # Select data
    t0= ISOdatetime(years[i],11,1,0,0,0)
    t1= ISOdatetime(years[i]+1,5,15,0,0,0)
    # Empty plot
    plot(c(t0,t1), c(mi,mx), type="n", xlab="", ylab="", main=paste(years[i],"/",years[i]+1,sep=""))
    # Plot data
    points(obs[,col_time][selected_obs], obs[,col_obs][selected_obs]*fact_obs, pch=20, col="darkblue")
    points(obs[,col_time][!selected_obs], obs[,col_obs][!selected_obs]*fact_obs, pch=20, col="lightblue")
    lines(sim[,col_time], sim[,col_sim], col="blue")
    # Legend
    if (i == 1) {
      legend("top",bty="n",lty=c(0,0,1),pch=c(20,20,NA),col=c("darkblue","lightblue","blue"),
        legend=c("Calib. data", "Other data", "Simul."))
    }
  }
  par(op)
}

