 
saturation <- function(wc_s, wc, var1, var2, frac1, frac2, frac3, frac4, frac5) {
  
  node <- vector(length=5)
  node[1] = wc_s - wc_s * var2
  node[2] = wc_s - wc_s * var1
  node[3] = wc_s
  node[4] = wc_s + wc_s * var1
  node[5] = wc_s + wc_s * var2
  frac <- c(frac1,frac2,frac3,frac4,frac5)
  
  if(wc <= node[1])
    return(0)
  
  tempalt <- 0
  testi = F
  for (i in 1:4) {
    if(testi == F){
      temp <- tempalt + 
              (node[i+1] - node[i]) * (frac[i+1] - frac[i]) / 2 + 
              (node[i+1] - node[i]) * (frac[5] - (frac[i+1] - frac[i]))
    
      if( (wc > tempalt) & (wc < temp) ) {
        frac_sat <- ( (frac[i] - frac[1]) + (wc - tempalt) / (temp - tempalt) * (frac[i+1] - frac[i]) ) * 
        testi = T
      }
      tempalt = temp
    }
  }
  
  if(exists("frac_sat"))
    return(frac_sat)
  else
    return(0)
}


saturation_new <- function(wc_s, wc, var1, var2, frac1, frac2, frac3, frac4, frac5) {
  
  node <- vector(length=5)
  node[1] = wc_s - wc_s * var2
  node[2] = wc_s - wc_s * var1
  node[3] = wc_s
  node[4] = wc_s + wc_s * var1
  node[5] = wc_s + wc_s * var2
  frac <- c(frac1,frac2,frac3,frac4,frac5)
  
  if(wc <= node[1])
    return(0)
  
  if(wc > node[5])
    return(1)
  
  for (i in 1:4) {
    if( wc > node[i] && wc <= node[i+1] )
      return( frac[i] + (frac[i+1] - frac[i]) * (1 - (node[i+1] - wc) / (node[i+1] - node[i])) )
  }
}


wc_sat = 0.6
v1 = 0.05
v2 = 0.1
f1 = 0
f2 = 0.1
f3 = 0.5
f4 = 0.9
f5 = 1

wc_a = seq(0.5,0.7, by=0.01)
h=0
res=vector(length=length(wc_a))
res_new = vector(length=length(wc_a))
for (i in wc_a) {
  h = h+1
  res[h] <- saturation(wc_sat, i, v1, v2, f1, f2, f3, f4, f5)
  res_new[h] <- saturation_new(wc_sat, i, v1, v2, f1, f2, f3, f4, f5)
}
