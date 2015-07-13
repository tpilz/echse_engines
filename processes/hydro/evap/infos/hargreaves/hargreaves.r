Hargreaves= function(Tmax,Tmin,Ra) {
    CH=0.0023 # Emperical Hargreaves coefficient
    CT= 17.8 # Emperical temperature coefficient
    return(CH*(((Tmax+Tmin)/2) + CT)*(sqrt(Tmax-Tmin))*Ra)
}

Tmin=29
Tmax=33
Rad= 1.28

Hargreaves(Tmax,Tmin,Rad)



