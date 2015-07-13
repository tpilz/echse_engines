
// Rate of precipitation input in m3/s
double p= inputExt(precip) * paramNum(fac_precip) / 1000. / delta_t *
  paramNum(area_max);

// Rate of evaporation loss in m3/s
double e= lakeEvap_makkink(inputExt(tavg), inputExt(glorad)) *
  paramFun(h2a, paramFun(v2h,u[INDEX_v]));

// Estimation of inflow rates at begin and end of the time step with priority
// on a proper mass balance
double qi_ini= max(0., 2 * inputSim(qi_avg) - inputSim(qi_end));
double qi_end= 2 * inputSim(qi_avg) - qi_ini;

// Derivatives with respect to time

// Storage volume of the lake
dudt[INDEX_v]= (qi_ini + (qi_end-qi_ini) * t/delta_t) // Inflow (linear f(time))
  + p                                                 // Precip
  - paramFun(h2q, paramFun(v2h,u[INDEX_v]))           // Outflow
  - e;                                                // Evaporation

// Volumes of precipitation and evaporation
dudt[INDEX_vp]= p;
dudt[INDEX_ve]= e;
