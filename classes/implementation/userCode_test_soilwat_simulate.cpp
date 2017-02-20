
// save initial values of runoff components
double surf_sat_init = stateScal(runst_surf_sat);
double surf_inf_init = stateScal(runst_surf_inf);
double sub_init = stateScal(runst_sub);
double gw_init = stateScal(runst_gw);

// scaling of ksat
vector<double> ksat_scale(stateVect(wc).size());
for (unsigned int i=0; i<stateVect(wc).size(); i++)
	ksat_scale[i] = paramFun(ksat,i+1) / sharedParamNum(scale_ks);

// compute initial states of hydraulic properties at very first time step (-9999. given as initial value; no need to compute within pre-processing)
// vector<double> mat_pot_t(stateVect(wc).size());
// vector<double> ku_t(stateVect(wc).size());
// vector<double> wc_a(stateVect(wc).size());
// if( abs(stateVect(mat_pot)[0] + 9999.) < 0.01 ){
// 	// get actual soil moisture states
// 	wc_a = stateVect(wc);
// 	
// 	// calculate matric potential
// 	for (unsigned int i=0; i<stateVect(wc).size(); i++)
// 		mat_pot_t[i] = matric_pot(sharedParamNum(ch_soilmod), wc_a[i], paramFun(wc_sat,i+1), paramFun(wc_res,i+1), paramFun(pores_ind,i+1), paramFun(bubble,i+1), sharedParamNum(na_val));
// 
// 	set_stateVect(mat_pot) = mat_pot_t;
// 	
// 	// calculate conductivity
// 	for (unsigned int i=0; i<stateVect(wc).size(); i++)
// 		ku_t[i] = k_unsat(sharedParamNum(ch_soilmod), wc_a[i], paramFun(wc_sat,i+1), paramFun(wc_res,i+1), paramFun(pores_ind,i+1), paramFun(ksat,i+1), sharedParamNum(na_val));
// 
// 	set_stateVect(k_u) = ku_t;
// }

// vector with all state variables
vector<double> states_all;
states_all.insert(states_all.end(), stateScal_all().begin(), stateScal_all().end());
// Put wc and hydraulic properties at the end to access the other variables via their pre-defined index in derivsScal()
states_all.insert(states_all.end(), stateVect(wc).begin(), stateVect(wc).end());
// states_all.insert(states_all.end(), stateVect(mat_pot).begin(), stateVect(mat_pot).end());
// states_all.insert(states_all.end(), stateVect(k_u).begin(), stateVect(k_u).end());

// update soil water storage by Numerical ODE integration (see *serivsScal.cpp)
odesolve_nonstiff(
  states_all,    // Initial value(s) of state variable(s)
  delta_t,            // Length of time step
  sharedParamNum(ode_accuracy),// Accuracy (in case of soil moisture 1e-3 is equivalent to 1 mm for a horizon of 1 m thickness)
  sharedParamNum(ode_max_iter),// Maximum number of sub-steps
  this,               // Pointer to this object
  states_all,// New value(s) of state variable(s)
	sharedParamNum(ch_odesolve) // choice flag of method for numerical integration
);


// UPDATE STORAGES
vector<double> wc_a(stateVect(wc).size());
vector<double> state_t(stateScal_all().size());
for (unsigned int i=0; i<stateScal_all().size(); i++)
	state_t[i] = states_all[i];

set_stateScal_all() = state_t;

for (unsigned int i=0; i<stateVect(wc).size(); i++)
	wc_a[i] = states_all[stateScal_all().size()+i];

set_stateVect(wc) = wc_a;


// soil hydraulic head / matric potential / capillary suction (m of water)/(100 hPa)
vector<double> mat_pot_t(stateVect(wc).size());
for (unsigned int i=0; i<stateVect(wc).size(); i++)
	mat_pot_t[i] = matric_pot(sharedParamNum(ch_soilmod), wc_a[i], paramFun(wc_sat,i+1), paramFun(wc_res,i+1), paramFun(pores_ind,i+1), paramFun(bubble,i+1), sharedParamNum(na_val));

set_stateVect(mat_pot) = mat_pot_t;

// hydraulic conductivity (m/s)
vector<double> ku_t(stateVect(wc).size());
for (unsigned int i=0; i<stateVect(wc).size(); i++)
	ku_t[i] = k_unsat(sharedParamNum(ch_soilmod), wc_a[i], paramFun(wc_sat,i+1), paramFun(wc_res,i+1), paramFun(pores_ind,i+1), ksat_scale[i], sharedParamNum(na_val));

set_stateVect(k_u) = ku_t;



// OUTPUT (m/s)
// calculate averages of wc and wc_sat
double w = -9999.;
double wc_sat_av = 0.;
double wc_a_av = 0.;
for (unsigned int i=0; i<wc_a.size(); i++) {
	w = paramFun(hor_depth,i+1)/paramNum(soil_depth);
	wc_sat_av += paramFun(wc_sat,i+1) * w;
	wc_a_av += wc_a[i] * w;
}
// calculate saturation of soil (-)
double f_sat = f_saturation(
	wc_a_av,
	wc_sat_av,
	sharedParamNum(var1),
	sharedParamNum(var2),
	sharedParamNum(var3),
	sharedParamNum(var4),
	sharedParamNum(var5),
	sharedParamNum(frac1),
	sharedParamNum(frac2),
	sharedParamNum(frac3),
	sharedParamNum(frac4),
	sharedParamNum(frac5)
);

set_output(saturation) = f_sat;

set_output(run_surf_sat) = (stateScal(runst_surf_sat) - surf_sat_init) / delta_t;
set_output(run_surf_inf) = (stateScal(runst_surf_inf) - surf_inf_init) / delta_t;
set_output(run_surf) = ( (stateScal(runst_surf_sat)+stateScal(runst_surf_inf)) - (surf_sat_init+surf_inf_init) ) / delta_t;

set_output(run_sub) = (stateScal(runst_sub) - sub_init) / delta_t ;

set_output(run_gw) = (stateScal(runst_gw) - gw_init) / delta_t ;
