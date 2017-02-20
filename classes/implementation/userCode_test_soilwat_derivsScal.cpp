
// INTERNAL VARIABLES AND CONSTANTS

// NA values
const double naval = sharedParamNum(na_val);

// surface water input (m/s)
double in = inputExt(inflow);

// no. of horizons (u contains scalar state variables (5) and the water content vector)
const unsigned int ns = stateScal_all().size();
const unsigned int nh = u.size()-ns;

// sum of flows (i.e. in- and outflows for each horizon to update soil moisture) (m/s)
vector<double> flows(nh, 0.);

// groundwater recharge (m/s)
double gw_rch = 0.;

// scaling of ksat
vector<double> ksat_scale(nh);
for (unsigned int i=0; i<nh; i++)
	ksat_scale[i] = paramFun(ksat,i+1) / sharedParamNum(scale_ks);

// hydraulic properties based on current soil moisture state
vector<double> ku(nh+1, -9999.);
vector<double> mat_pot(nh+1, -9999.);
for (unsigned int i=0; i<nh; i++) {
	//mat_pot[i] = u[ns+nh+i];
	mat_pot[i] = matric_pot(sharedParamNum(ch_soilmod), u[ns+i], paramFun(wc_sat,i+1), paramFun(wc_res,i+1), paramFun(pores_ind,i+1), paramFun(bubble,i+1), naval);
	//ku[i] = u[ns+(2*nh)+i];
	ku[i] = k_unsat(sharedParamNum(ch_soilmod), u[ns+i], paramFun(wc_sat,i+1), paramFun(wc_res,i+1), paramFun(pores_ind,i+1), ksat_scale[i], naval);
}
// lower boundary conditions (at the moment free drainage, i.e. unit gradient (water movement due to gravitation) and persistency of conductivity)
ku[nh] = ku[nh-1];
mat_pot[nh] = mat_pot[nh-1];

//////////////////////////////////////////////////////////
// SURFACE RUNOFF
//////////////////////////////////////////////////////////
// calculate averages of wc and wc_sat
double w = -9999.;
double wc_sat_av = 0.;
double wc_a_av = 0.;
for (unsigned int i=0; i<nh; i++) {
	w = paramFun(hor_depth,i+1)/paramNum(soil_depth);
	wc_sat_av += paramFun(wc_sat,i+1) * w;
	wc_a_av += u[ns+i] * w;
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

// calculate saturation excess runoff (m/s) and input remaining for infiltration (m/s)
double surf_sat = f_sat * in;
in *= (1. - f_sat);

// infiltration (m/s)
// double inf = infiltration(
// // General input and parameters
// 	sharedParamNum(ch_inf),			// Choice flag for selection of method
// 	in,													// Water hitting top of soil to be infiltrated (ms-1)
// 	u[ns+0],										// Actual volumetric water content at start of infiltration (m3/m3)
// 	paramFun(wc_sat,1.),					// Volumetric water content at saturation (m3/m3)
// 	paramFun(ksat,1.),						// Saturated hydraulic conductivity (ms-1)
// 	delta_t,										// time step length (s)
// 	naval,											// NA value
// // Horton-specific parameters
// 	paramNum(Hort_ini),					// Horton parameter: initial infiltration rate (m/s)
// 	paramNum(Hort_end),					// Horton parameter: final infiltration rate (m/s)
// 	paramNum(Hort_k),						// Horton parameter: decay constant (1/s)
// // Philip specific parameters
// 	paramNum(Phil_s),						// Philip parameter: Sorptivity (ms^{-1/2}); calculated internally if set to NA
// 	paramNum(Phil_a),						// Philip parameter: second term parameter (m/s); calculated internally if set to NA
// 	sharedParamNum(Phil_cal),		// Philip parameter: calibration paremter for Phil_a if this is NA (fraction of ksat), should be within [0.2..1.0]
// // Green-Ampt after Peschke specific parameters
// 	paramFun(suc,1.)							// Suction at wetting front (m)
// );



// account for case of wetting front depth being larger than top soil layer (distribute excess infiltration into deeper horizons)
// double inf_t = -9999.;
// double inf_cum = 0.;
// double vol_ref = -9999.;
// for (unsigned int i=0; i<nh; i++) {
// 	
// 	//cout << "Horizon " << i << ", wc = " << u[ns+i] << endl;
// 	
// 	// calculate infiltration (m/s)
// 	inf_t = infiltration(
// 	// General input and parameters
// 		sharedParamNum(ch_inf),			// Choice flag for selection of method
// 		in,													// Water hitting top of soil to be infiltrated (ms-1)
// 		u[ns+i],										// Actual volumetric water content at start of infiltration (m3/m3)
// 		paramFun(wc_sat,i+1),					// Volumetric water content at saturation (m3/m3)
// 		paramFun(ksat,i+1),						// Saturated hydraulic conductivity (ms-1)
// 		delta_t,										// time step length (s)
// 		naval,											// NA value
// 	// Horton-specific parameters
// 		paramNum(Hort_ini),					// Horton parameter: initial infiltration rate (m/s)
// 		paramNum(Hort_end),					// Horton parameter: final infiltration rate (m/s)
// 		paramNum(Hort_k),						// Horton parameter: decay constant (1/s)
// 	// Philip specific parameters
// 		paramNum(Phil_s),						// Philip parameter: Sorptivity (ms^{-1/2}); calculated internally if set to NA
// 		paramNum(Phil_a),						// Philip parameter: second term parameter (m/s); calculated internally if set to NA
// 		sharedParamNum(Phil_cal),		// Philip parameter: calibration paremter for Phil_a if this is NA (fraction of ksat), should be within [0.2..1.0]
// 	// Green-Ampt after Peschke specific parameters
// 		paramFun(suc,i+1)							// Suction at wetting front (m)
// 	);
// 	
// 	// calculate refillable volume of soil layer (m)
// 	vol_ref = max( 0., (paramFun(wc_sat,i+1) - u[ns+i]) * paramFun(hor_depth,i+1) );
// 	
// 	// infiltration of current iteration must not be larger than wat is left to be infiltrated (m/s)
// 	inf_t = min(in - inf_cum, inf_t);
// 	
// 	// compare with calculated amount of infiltration and calculate excess for next horizon if necessary (m/s)
// 	if ( (inf_t * delta_t) > vol_ref ) {
// 		inf_cum += vol_ref/delta_t;
// 		// update sum of flows for ith horizon
// 		flows[i] += vol_ref/delta_t;
// 		
// 		// if this is the last horizon all remaining water goes into groundwater recharge storage
// 		if (i == (nh-1)) {
// 			gw_rch += inf_t - vol_ref/delta_t;
// 			inf_cum += inf_t - vol_ref/delta_t;
// 			break;
// 		}
// 		
// 	} else { // no excess for next horizon -> update sum of flows and leave loop
// 		inf_cum += inf_t;
// 		flows[i] += inf_t;
// 		break;
// 	}
// 	
// // 	// calculate infiltration of excess flow for next horizon as if it was a top horizon (to account for different soil properties) (m/s)
// // 	 inf_t = infiltration(
// // 	// General input and parameters
// // 		sharedParamNum(ch_inf),			// Choice flag for selection of method
// // 		in,													// Water hitting top of soil to be infiltrated (ms-1)
// // 		u[ns+i+1],									// Actual volumetric water content at start of infiltration (m3/m3)
// // 		paramFun(wc_sat,i+2),				// Volumetric water content at saturation (m3/m3)
// // 		paramFun(ksat,i+2),					// Saturated hydraulic conductivity (ms-1)
// // 		delta_t,										// time step length (s)
// // 		naval,											// NA value
// // 	// Horton-specific parameters
// // 		paramNum(Hort_ini),					// Horton parameter: initial infiltration rate (m/s)
// // 		paramNum(Hort_end),					// Horton parameter: final infiltration rate (m/s)
// // 		paramNum(Hort_k),						// Horton parameter: decay constant (1/s)
// // 	// Philip specific parameters
// // 		naval,											// Philip parameter: Sorptivity (ms^{-1/2}); calculated internally if set to NA
// // 		naval,											// Philip parameter: second term parameter (m/s); calculated internally if set to NA
// // 		sharedParamNum(Phil_cal),		// Philip parameter: calibration paremter for Phil_a if this is NA (fraction of ksat), should be within [0.2..1.0]
// // 	// Green-Ampt after Peschke specific parameters
// // 		paramFun(suc,i+2)						// Suction at wetting front (m)
// // 	);
// }
/*
// calculate infiltration excess runoff (m/s)
double surf_inf = in - inf_cum;*/

// cout << "ns = " << ns << endl;
// cout << "nh = " << nh << endl;
// cout << "mat_pot = " << mat_pot[0] << " m" << endl;
// cout << "ku = " << ku[0] << " m/s" << endl;
// cout << "infiltration input:" << endl;
// cout << "water soil surface flux: " << in << " m/s" << endl;
// cout << "soil moisture: " << u[ns] << endl;

// calculate infiltration (m/s)
double inf = infiltration(
// General input and parameters
	sharedParamNum(ch_inf),			// Choice flag for selection of method
	in,													// Water hitting top of soil to be infiltrated (ms-1)
	u[ns],											// Actual volumetric water content at start of infiltration (m3/m3)
	paramFun(wc_sat,1),					// Volumetric water content at saturation (m3/m3)
	ksat_scale[0],						// Saturated hydraulic conductivity (ms-1)
	delta_t,										// time step length (s)
	naval,											// NA value
// Horton-specific parameters
	paramNum(Hort_ini),					// Horton parameter: initial infiltration rate (m/s)
	paramNum(Hort_end),					// Horton parameter: final infiltration rate (m/s)
	paramNum(Hort_k),						// Horton parameter: decay constant (1/s)
// Philip specific parameters
	paramNum(Phil_s),						// Philip parameter: Sorptivity (ms^{-1/2}); calculated internally if set to NA
	paramNum(Phil_a),						// Philip parameter: second term parameter (m/s); calculated internally if set to NA
	sharedParamNum(Phil_cal),		// Philip parameter: calibration paremter for Phil_a if this is NA (fraction of ksat), should be within [0.2..1.0]
// Green-Ampt after Peschke specific parameters
	paramFun(suc,1)							// Suction at wetting front (m)
);

flows[0] += inf;

// calculate infiltration excess runoff (m/s)
double surf_inf = in - inf;



//////////////////////////////////////////////////////////
// SUB-SURFACE RUNOFF
//////////////////////////////////////////////////////////

// lateral flows for every horizon
double subsurf = 0.;
double perc = -9999.;
double latfl = -9999.;
double r_sum_flows = -9999.;
double r_water_free = -9999.;
double hor_depth_next = -9999.;
for (unsigned int i=0; i<nh; i++) {
	
	// percolation from ith horizon (m/s)
	if(i == (nh-1)) // always allow percolation into vadose zone under soil profile
		hor_depth_next = 0.;
	else 
		hor_depth_next = paramFun(hor_depth,i+2);
	perc = percolation(sharedParamNum(ch_perc), paramFun(hor_depth,i+1), ku[i], u[ns+i],  paramFun(wc_fc,i+1), delta_t,
										 ku[i+1], mat_pot[i], mat_pot[i+1], hor_depth_next);
	
	
	// lateral outflow from ith horizon (m/s)
	latfl = latflow(paramFun(hor_depth,i+1), u[ns+i], paramFun(wc_fc,i+1), paramFun(wc_sat,i+1),
										ksat_scale[i], paramNum(slopelength), paramNum(slope));
	
	
	// check outflows (latflow and percolation); only saturation exccess (wc-wc_fc) can flow out due to gravitation
	// based on Guentner (2002) eqs. 4.47 and 4.48; ATTENTION: the latter is erroneous (but correct in WASA code)!
	r_sum_flows = (latfl + perc) * delta_t; // estimated total outflow
	if( r_sum_flows > 0. ) {
		// maximum possible outflow
		r_water_free = paramFun(hor_depth,i+1) * (u[ns+i] - paramFun(wc_fc,i+1));
		// if estimated outflow greater than possible reduce outflows accordingly
		if(r_sum_flows > r_water_free) {
			perc *= r_water_free / r_sum_flows;
			latfl *= r_water_free / r_sum_flows;
		}
	}
	
	
	// add to flow vector (percolation to groundwater recharge if last horizon)
	subsurf += latfl;
	if (i == (nh-1)) {
		flows[i] -= (perc+latfl);
		gw_rch += perc;
	} else {
		flows[i] -= (perc+latfl);
		flows[i+1] += perc;
	}
}

// for (unsigned int i=0; i<nh; i++) 
// 	cout << "Sum flows hor " << i << ": " << flows[i] << " m/s." << endl;

//////////////////////////////////////////////////////////
// SAVE DERIVATIVES
//////////////////////////////////////////////////////////
dudt[INDEX_runst_surf_sat] = surf_sat;
dudt[INDEX_runst_surf_inf] = surf_inf;
dudt[INDEX_runst_surf] = surf_sat+surf_inf;
dudt[INDEX_runst_sub] = subsurf;
dudt[INDEX_runst_gw] = gw_rch;

for (unsigned int i=0; i<nh; i++)
	dudt[ns+i] = flows[i] / paramFun(hor_depth,i+1);
