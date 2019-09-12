
// INITIALS AND INTERNAL VARIABLES

// set evapotranspiration states to zero
set_stateScal(et_p) = 0.;
set_stateScal(et_a) = 0.;
set_stateScal(et_i) = 0.;
set_stateScal(r_interc) = 0.;

// set runoff states to zero
set_stateScal(runst_surf_sat) = 0.;
set_stateScal(runst_surf_inf) = 0.;
set_stateScal(runst_surf) = 0.;
set_stateScal(runst_sub) = 0.;
set_stateScal(runst_gw) = 0.;

// apply calibration parameters
vector<double> pores_ind_c(stateVect(wc).size());
for (unsigned int i=0; i<stateVect(wc).size(); i++) {
	pores_ind_c[i] = min(1., paramFun(pores_ind,i+1) * sharedParamNum(cal_pores) );
}

// scaling of ksat (as in WASA)
vector<double> ksat_scale(stateVect(wc).size());
double kfcorr_a = -9999.;
double kfcorr_b = -9999.;
double kfcorr = sharedParamNum(scale_ks);
// implicit choice of scaling method (depending on precipitation or not)
if( (sharedParamNum(scale_ks_a) > 1e-6) && (delta_t > 3600) ) {
	kfcorr_a = sharedParamNum(scale_ks_a);
	kfcorr_b = sharedParamNum(scale_ks_b);
} else {
	kfcorr_a = 0.;
	kfcorr_b = 0.;
}
// scaling depending on precipitation in case of daily resolution 
if(inputExt(precip) > 1e-6) {
	kfcorr = kfcorr * ( kfcorr_a/inputExt(precip) + kfcorr_b + 1 );
} else {
	kfcorr = 1.;
}
// update ksat values (and apply calibration factor)
for (unsigned int i=0; i<stateVect(wc).size(); i++)
	ksat_scale[i] = paramFun(ksat,i+1) * sharedParamNum(cal_ks) / kfcorr;

// in case of SUB-DAILY temporal resolution: assess energy budget over nighttimes in case energy balance components have to be calculated
double i_radex = inputExt(radex);
double i_glorad_max = inputExt(glorad_max);
double i_longrad = inputExt(rad_long);
if(delta_t < 86400) {
	// account for sub-daily resolution in calculation of energy budget (i.e. cloudiness correction factor); only needed for long-wave radiation calculation
	if ( (abs(i_longrad - sharedParamNum(na_val)) < 0.01) && (abs(inputExt(rad_net) - sharedParamNum(na_val)) < 0.01) ) {
		// calculate short-wave radiation under clear sky if not given
		if ( abs(i_glorad_max - sharedParamNum(na_val)) < 0.01 ) {
			// calculate radex if not given
			if ( abs(i_radex - sharedParamNum(na_val)) < 0.01 )
				i_radex = rad_extraterr_hourly(inputExt(doy), paramNum(lat), inputExt(hour), inputExt(utc_add), paramNum(lon));
			
			i_radex = max(i_radex, inputExt(glorad)); // minimum equal to measured glorad; might be smaller during sunrise/sunset hours due to uncertainties in calculation
			i_glorad_max = calc_glorad_max(sharedParamNum(choice_gloradmax), i_radex, sharedParamNum(radex_a), sharedParamNum(radex_b), paramNum(elev));
			i_glorad_max = max(i_glorad_max, inputExt(glorad)); 
		}
		// if i_glorad_max is very small assume nighttime and do not update state of longrad (i.e. assume persistance of cloudiness correction factor over nighttime)
		// take 5 W/m2 as "small" (although this is somewhat arbitrary) to avoid impact of atmospheric disturbances during low sun angle
		if ( i_glorad_max > 5. ) {
			i_longrad = net_longrad(inputExt(temper),inputExt(rhum),inputExt(glorad),i_glorad_max,
																						sharedParamNum(emis_a),sharedParamNum(emis_b),
																						sharedParamNum(fcorr_a),sharedParamNum(fcorr_b));
		} else
			i_longrad = stateScal(s_longrad);
	}
}
set_stateScal(s_longrad) = i_longrad;

// vector with state variables
vector<double> states_all;
states_all.insert(states_all.end(), stateScal_all().begin(), stateScal_all().end());
// Put wc at the end to access the other variables via their pre-defined index in derivsScal()
states_all.insert(states_all.end(), stateVect(wc).begin(), stateVect(wc).end());




// SOIL WATER PROCESSES AND NUMERICAL INTEGRATION

// update soil water storage by Numerical ODE integration (see *serivsScal.cpp)
odesolve_nonstiff(
  states_all,    // Initial value(s) of state variable(s)
  delta_t,            // Length of time step
  sharedParamNum(ode_accuracy),// Accuracy (in case of soil moisture 1e-3 is equivalent to 1 mm for a horizon of 1 m thickness)
  sharedParamNum(ode_max_iter),// Maximum number of sub-steps
  this,               // Pointer to this object
  states_all,// New value(s) of state variable(s)
	sharedParamNum(choice_odesolve) // choice flag of method for numerical integration
);




// UPDATE STORAGES
vector<double> wc_a(stateVect(wc).size());
vector<double> state_t(stateScal_all().size());
for (unsigned int i=0; i<stateScal_all().size(); i++) {
	state_t[i] = states_all[i];
	//cout << "states_all[" << i << "] = " << states_all[i] << endl;
}

set_stateScal_all() = state_t;

// runst_sub sometimes is negative although values of subsurf in derivsScal() are positive; don't know the reason, ignore this as long as value is negligible
// if(stateScal(runst_sub) < -1e-6) {
// 	stringstream errmsg;
// 	errmsg << "runst_sub is significantly negative (runst_sub = " << stateScal(runst_sub) << ")!";
// 	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 	throw(e); 
// }

// actual water content (m3/m3)
for (unsigned int i=0; i<stateVect(wc).size(); i++) {
	wc_a[i] = states_all[stateScal_all().size()+i];
	// limit water content within physical boundaries
	if(sharedParamNum(choice_constraint) > 0.01 ) {
		wc_a[i] = max(wc_a[i], paramFun(wc_res,i+1));
		wc_a[i] = min(wc_a[i], paramFun(wc_sat,i+1));
	}
}

set_stateVect(wc) = wc_a;


// soil hydraulic head / matric potential / capillary suction (m of water)/(100 hPa)
vector<double> mat_pot_t(stateVect(wc).size());
for (unsigned int i=0; i<stateVect(wc).size(); i++)
	mat_pot_t[i] = matric_pot(sharedParamNum(choice_soilmod), wc_a[i], paramFun(wc_sat,i+1), paramFun(wc_res,i+1), pores_ind_c[i], paramFun(bubble,i+1), sharedParamNum(na_val));

set_stateVect(mat_pot) = mat_pot_t;

// hydraulic conductivity (m/s)
vector<double> ku_t(stateVect(wc).size());
for (unsigned int i=0; i<stateVect(wc).size(); i++)
	ku_t[i] = k_unsat(sharedParamNum(choice_soilmod), wc_a[i], paramFun(wc_sat,i+1), paramFun(wc_res,i+1), pores_ind_c[i], ksat_scale[i], sharedParamNum(na_val));

set_stateVect(k_u) = ku_t;

// prevent interception storage from becoming negative
if(sharedParamNum(choice_constraint) > 0.01 ) {
	set_stateScal(v_interc) = max(0., stateScal(v_interc));
}



// OUTPUT; area-weighted, for investigation divide value by areal fraction of SVC in TC (m/s)

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

// plant available water within rooted zone (m)
double const rootdepth = min(inputExt(rootd)*sharedParamNum(cal_rootd), paramNum(soil_depth));
double wc_root = 0.;
double wcp_root = 0.;
double cum_depth = 0.;
double wv_plant = -9999.;
double wv_soil = 0.;
double wc_plant = 0.;
vector<double> wc_actplant(wc_a.size());
vector<double> w_root(wc_a.size());
vector<double> w_eta_t(wc_a.size());
for (unsigned int i=0; i<wc_a.size(); i++) {
	// weight according to root depth and depth of current horizon; accounting fraction of layer containing roots (in last rooted horizon)
	cum_depth += paramFun(hor_depth, i+1);
	if (rootdepth > cum_depth)
		w_root[i] = 1.;
	else if( (cum_depth - rootdepth) < paramFun(hor_depth, i+1) )
		w_root[i] = (1. - (cum_depth - rootdepth) / paramFun(hor_depth, i+1));
	else
		w_root[i] = 0.;
	// calc plant available water content and parameters for root zone
	wc_actplant[i] = min(paramFun(wc_fc, i+1), wc_a[i]); // plant available water content limited to field capacity
	wc_plant += max(0., wc_actplant[i] - paramFun(wc_pwp, i+1)) * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
	wc_root += wc_a[i] * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
	wcp_root += paramFun(wc_pwp, i+1) * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
	// total water in soil profile
	wv_soil += wc_a[i] * paramFun(hor_depth, i+1);
}
// actual plant available water (m)
wv_plant = wc_plant * rootdepth;
// calc fraction of total plant available water for each horizon (-)
double plant_wat = -9999.;
for (unsigned int i=0; i<wc_a.size(); i++) {
	if (wv_plant < 1e-6)
		if (i == 0)
			w_eta_t[i] = 1.;
		else
			w_eta_t[i] = 0.;
	else {
		plant_wat =  max(0., (wc_actplant[i] - paramFun(wc_pwp, i+1)) * w_root[i] * paramFun(hor_depth, i+1) );
		w_eta_t[i] = plant_wat / wv_plant;
	}
}

double vect_sum = std::accumulate(w_eta_t.begin(), w_eta_t.end(), 0.);
if( abs(vect_sum - 1.) > 1e-12 ) {
	stringstream errmsg;
	errmsg << "Sum of w_eta = " << vect_sum << " but should be equal to one!";
	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
	throw(e); 
}
set_output(v_plantwat) = wv_plant * paramNum(frac_area);
set_output(v_soilwat) = wv_soil * paramNum(frac_area);
set_stateVect(w_eta) = w_eta_t;

// evapotranspiration (m/s)
set_output(eta) = stateScal(et_a) / delta_t * paramNum(frac_area);
set_output(etp) = stateScal(et_p) / delta_t * paramNum(frac_area);
set_output(eti) = stateScal(et_i) / delta_t * paramNum(frac_area);
set_output(interc) = stateScal(r_interc) / delta_t * paramNum(frac_area);

// surface runoff (m/s)
// lateral TC outflow
set_output(run_surf_sat_tc) = stateScal(runst_surf_sat) / delta_t * paramNum(frac_area) * paramNum(frac_area);
set_output(run_surf_inf_tc) = stateScal(runst_surf_inf) / delta_t * paramNum(frac_area) * paramNum(frac_area);
set_output(run_surf_tc) = (stateScal(runst_surf_sat)+stateScal(runst_surf_inf)) / delta_t * paramNum(frac_area) * paramNum(frac_area);
// lateral SVC outflow (redistribution within TC to neighbour SVCs)
set_output(run_surf_sat_svc) = stateScal(runst_surf_sat) / delta_t * (1. - paramNum(frac_area)) * paramNum(frac_area);
set_output(run_surf_inf_svc) = stateScal(runst_surf_inf) / delta_t * (1. - paramNum(frac_area)) * paramNum(frac_area);
set_output(run_surf_svc) = (stateScal(runst_surf_sat)+stateScal(runst_surf_inf)) / delta_t * (1. - paramNum(frac_area)) * paramNum(frac_area);

// subsurface runoff (m/s)
// lateral TC outflow
set_output(run_sub_tc) = stateScal(runst_sub) / delta_t * paramNum(frac_area) * paramNum(frac_area);
// lateral SVC outflow (redistribution within TC to neighbour SVCs)
set_output(run_sub_svc) = stateScal(runst_sub) / delta_t * (1. - paramNum(frac_area)) * paramNum(frac_area);

// percolation out of soil profile; area-weighted for summation in TC (m/s)
set_output(run_gw) = stateScal(runst_gw) / delta_t * paramNum(frac_area);
