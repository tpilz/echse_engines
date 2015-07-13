
computeGOF = function (obs, sim, col_sim, col_obs, col_time, fact_obs, selected_obs) {
  # Filter
  obs= obs[selected_obs,]
  # Make column names unique
  names(sim)[which(names(sim) == col_sim)]= paste("sim",col_sim,sep="_")
  names(obs)[which(names(obs) == col_obs)]= paste("obs",col_obs,sep="_")
  col_sim= paste("sim",col_sim,sep="_")
  col_obs= paste("obs",col_obs,sep="_")
  # Find corresponding data by merging
  all= merge(x=sim, y=obs, by=col_time, all=FALSE)
  if (nrow(all) == 0) stop("No corresponding data.")
  # Compute GOF
  return( sqrt(mean((all[,col_sim] - all[,col_obs]*fact_obs)^2)) )
}

