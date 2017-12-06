//////////////////////////////////////////////////////////
// INTERNAL VARIABLES AND CONSTANTS
//////////////////////////////////////////////////////////
// convert precipitation sum into a flux (averaged over timestep) (mm -> m/s)
const double prec_flux = inputExt(precip) / 1000. / delta_t;

// NA values
const double naval = sharedParamNum(na_val);

// no. of horizons (u contains scalar state variables and vectors states, i.e. water content, conductivity, and matric potential)
const unsigned int ns = stateScal_all().size();
const unsigned int nh = u.size()-ns;

// sum of flows (i.e. in- and outflows for each horizon to update soil moisture) (m/s)
vector<double> flows(nh, 0.);

// groundwater recharge (vertical percolation out of soil profile) (m/s)
double gw_rch = 0.;

// sub-surface outflow (m/s)
double subsurf = 0.;

// scaling of ksat (as in WASA)
vector<double> ksat_scale(nh);
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
	kfcorr = kfcorr * ( sharedParamNum(scale_ks_a)/inputExt(precip) + sharedParamNum(scale_ks_b) + 1 );
} else {
	kfcorr = 1.;
}
// update ksat values (and apply calibration factor)
for (unsigned int i=0; i<nh; i++)
	ksat_scale[i] = paramFun(ksat,i+1) * sharedParamNum(cal_ks) / kfcorr;

// hydraulic properties based on current soil moisture state
vector<double> ku(nh+1, -9999.);
vector<double> mat_pot(nh+1, -9999.);
for (unsigned int i=0; i<nh; i++) {
	mat_pot[i] = matric_pot(sharedParamNum(choice_soilmod), u[ns+i], paramFun(wc_sat,i+1), paramFun(wc_res,i+1), paramFun(pores_ind,i+1), paramFun(bubble,i+1), naval);
	ku[i] = k_unsat(sharedParamNum(choice_soilmod), u[ns+i], paramFun(wc_sat,i+1), paramFun(wc_res,i+1), paramFun(pores_ind,i+1), ksat_scale[i], naval);
}
// lower boundary conditions (at the moment free drainage, i.e. unit gradient (water movement due to gravitation) and persistency of conductivity)
ku[nh] = ku[nh-1];
mat_pot[nh] = mat_pot[nh-1];

// set effective root depth to total soil depth
double const rootdepth = min(inputExt(rootd), paramNum(soil_depth));

// water contents (-)
// top soil water content
double wc_top = u[ns];
// actual water content and averaged parameters of root zone (according to actual root depth)
double wc_root = 0.;
double wcs_root = 0.;
double wcf_root = 0.;
double wcp_root = 0.;
double wcr_root = 0.;
double wv_plant = 0.;
double wc_plant = 0.;
double bubble_root = 0.;
double porei_root = 0.;
double cum_depth = 0.;
vector<double> w_root(nh);
vector<double> w_eta(nh); // needed later to relate total evapotranspiration to horizons
for (unsigned int i=0; i<nh; i++) {
	// weight according to root depth and depth of current horizon; accounting fraction of layer containing roots (in last rooted horizon)
	cum_depth += paramFun(hor_depth, i+1);
	if (rootdepth > cum_depth)
		w_root[i] = 1.;
	else if( (cum_depth - rootdepth) < paramFun(hor_depth, i+1) )
		w_root[i] = (1. - (cum_depth - rootdepth) / paramFun(hor_depth, i+1));
	else
		w_root[i] = 0.;
	// calc water content and parameters for root zone
	wc_plant += max(0., u[ns+i] - paramFun(wc_pwp, i+1)) * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
	wc_root += u[ns+i] * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
	wcs_root += paramFun(wc_sat, i+1) * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
	wcf_root += paramFun(wc_fc, i+1) * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
	wcp_root += paramFun(wc_pwp, i+1) * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
	wcr_root += paramFun(wc_res, i+1) * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
	bubble_root += paramFun(bubble, i+1) * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
	porei_root += paramFun(pores_ind, i+1) * w_root[i] * paramFun(hor_depth, i+1) / rootdepth;
}

// actual plant available water (m)
wv_plant = wc_plant * rootdepth;
// calc fraction of total plant available water for each horizon (-)
double plant_wat = -9999.;
for (unsigned int i=0; i<nh; i++) {
	if (wv_plant < 1e-6)
		if (i == 0)
			w_eta[i] = 1.;
		else
			w_eta[i] = 0.;
	else {
		plant_wat =  max(0., (u[ns+i] - paramFun(wc_pwp, i+1)) * w_root[i] * paramFun(hor_depth, i+1) );
		w_eta[i] = plant_wat / wv_plant;
	}
}




//////////////////////////////////////////////////////////
// SURFACE INFLOW (m/s)
//////////////////////////////////////////////////////////
// calculate interception of precipitation (m/s)
double r_inter = intercept(prec_flux, inputExt(lai), paramNum(intfc), max(0., u[INDEX_v_interc]), delta_t);

// calculate soil surface inflow (m/s)
double in = prec_flux - r_inter + inputSim(r_surf_in);



//////////////////////////////////////////////////////////
// SUB-SURFACE INFLOW (m/s)
// TODO: 
//  - Lateral subsurface inflows area-weighted added to
//    the soil layers. Better approach preserving information
//    of source layer?!
//  - More physically-based approach than just adding to flows,
//    e.g. taking k_u and mat_pot into account?
//////////////////////////////////////////////////////////
vector<double> w_depth(nh, -9999.);
double inflow = -9999.;
double vol_ref = -9999.;
// distribute inflow to soil layers
for (unsigned int i=0; i<nh; i++) {
	w_depth[i] = paramFun(hor_depth,i+1)/paramNum(soil_depth);
	inflow = inputSim(r_sub_in) * w_depth[i];
	// excess goes to subsurface outflow
	vol_ref = max( 0., (paramFun(wc_sat,i+1) - u[ns+i]) * paramFun(hor_depth,i+1) );
	if(inflow > vol_ref){
		flows[i] += vol_ref;
		subsurf += inflow - vol_ref;
	} else {
		flows[i] += inflow;
	}
}




//////////////////////////////////////////////////////////
// EVAPOTRANSPIRATION (m/s)
//////////////////////////////////////////////////////////
// potential evapotranspiration (m/s)
double r_etp = et_pot (
// Common meteorological variables
	inputExt(temper),									// Air temperature (°C)
	inputExt(temper_min),								// Minimum temperature within time step (°C)
	inputExt(temper_max),								// Maximum temperature within time step (°C)
  inputExt(glorad),   						// Downward short-wave radiation (W/m2)
	inputExt(rhum),										// relative humidity (%)
	inputExt(wind),										// wind speed (m/s)
	inputExt(apress),       						// Air pressure (hPa)
	inputExt(sundur),									// Sunshine duration of current day (h)
	inputExt(cloud),										// Cloudiness (%)
// Common site-secific parameters
	paramNum(lat),											// Latitude (decimal degree)
	paramNum(lon),											// Longitude of the location of interest (decimal degrees west of Greenwich, e.g. Greenwich: 0°, Berlin: 345°, New York: 75°)
	paramNum(elev),										// Elevation above sea level (m)
// Specific meteorological quantities (commonly calculated internally)
	inputExt(radex),									// Incoming extraterrestrial radiation (i.e. at top of atmosphere) (W/m2)
	inputExt(glorad_max), 							// Downward short-wave radiation under clear (cloudless) sky (Wm-2)
	inputExt(rad_net),										// net incoming ( (1-alb) * short-wave + long-wave) radiation (Wm-2)
	inputExt(rad_net_soil),									// net incoming (short-wave + long-wave) radiation hitting the soil surface (Wm-2)
	u[INDEX_s_longrad],									// Net incoming long-wave radiation (Wm-2)
	inputExt(soilheat),								// Soil heat flux (Wm-2)
	inputExt(totalheat),								// heat conduction into soil AND plants due to physical and biochemical energy storage (Wm-2)
// Meteorological parameters
	sharedParamNum(h_tempMeas),							// height of temperature measurement above ground (m)
	sharedParamNum(h_humMeas),								// height of humidity measurement above ground (psychrometer) (m)
	sharedParamNum(h_windMeas),							// height of windspeed measurement above ground (m)
	sharedParamNum(emis_a),									// Coefficient a for calculating net emissivity (-)
	sharedParamNum(emis_b),									// Coefficient b for calculating net emissivity (-)
	sharedParamNum(fcorr_a),									// Coefficient a for calculating cloudiness correction factor (-)
	sharedParamNum(fcorr_b),									// Coefficient b for calculating cloudiness correction factor (-)
	sharedParamNum(radex_a),									// Angstrom coefficient (share of radex on glorad under clouds) (-)
	sharedParamNum(radex_b),									// Angstrom coefficient (radex_a+radex_b = share of radex on glorad under clear sky) (-)
	sharedParamNum(f_day),										// soil heat flux calculation: Fraction of net_rad over daytime (in case of sub-daily application) (-)
	sharedParamNum(f_night),									// soil heat flux calculation: Fraction of net_rad over nighttime (in case of sub-daily application) (-)
// Vegetation and land-cover parameters and variables
	paramNum(crop_makk),     					// Crop-factor after Makkink (-)
	paramNum(crop_faoref),     				// Crop-factor for FAO reference approach (-)
	inputExt(cano_height),							// canopy height (m)
	inputExt(lai),											// leaf area index (m2/m2)
	inputExt(alb),											// albedo (-)
	sharedParamNum(ext),											// Canopy extinction coefficient, Beer's law (-) -  in original WASA model code set to 0.5
	paramNum(res_leaf_min),						// Plant-specific minimum (i.e. no stress occurs) stomatal resistance of a single leaf (sm-1)
	paramNum(soil_dens),								// bulk density of soil of topmost soil horizon (kg/m3)
	paramNum(glo_half),								// Solar radiation at which stomatal conductance is half of its maximum (W/m2) - in WASA model a value of 100
	sharedParamNum(res_b),										// Mean boundary layer resistance (sm-1) - in SW and WASA a value of 25
	sharedParamNum(drag_coef),								// Effective value of mean drag coef. of vegetative elements (-) - in SG and WASA a value of 0.07
	sharedParamNum(rough_bare),							// Roughness length of bare substrate (m) - in SW, SG and WASA a value of 0.01
	sharedParamNum(eddy_decay),							// Eddy diffusivity decay constant (-) - by SW and SG set to 2.5 for agricultural crops
	sharedParamNum(rss_a),										// Empirical constant for calculation of soil surface resistance (-) - in original code set to 26
	sharedParamNum(rss_b),										// Empirical constant for calculation of soil surface resistance (-) - in original code set to -1
// Computational parameters
	inputExt(doy),											// Current day of year
	inputExt(hour),												// hour of day in local time (including daylight daving time), range of [0..23]
	inputExt(utc_add),										// Deviation of local time zone from UTC; may vary over the year due to daylight saving time; range of [-12..14] (hours)
	naval,															// NA value
	delta_t,						// time step length (s)
// Choice flags
	sharedParamNum(choice_et),											// flag which method to use
	sharedParamNum(choice_rcs),											// Flag: canopy stomatal resistance by up-scaling of single leaf resistance, 1: Shuttlewort & Wallace (1985) eq. 19, 2: Saugier and Katerji (1991) eq. 4 used in WASA
	sharedParamNum(choice_roughLen), 								// Flag: Roughness length for momentum transfer, 1: SWAT manual (2011) eqs. 2:2.2.4, 2:2.2.5, 2: Shuttleworth & Gurney (1990) eq. 43
	sharedParamNum(choice_plantDispl),							// Flag: Displacement height for a plant, 1: SWAT manual (2011) eq. 2:2.2.7, 2: Shuttleworth & Gurney (1990) eq. 42
	sharedParamNum(choice_gloradmax)								// Flag: calculation of maximum incoming short-wave radiation (clear sky), 1: Angstroem, 2: Allen (2005), ASCE standard etp, eq. 19 (based on elevation)
);

// interception evaporation rate (m/s)
double r_eti = min(r_etp, u[INDEX_v_interc]/delta_t);

// actual evapotranspiration (m/s)
double r_eta = -9999.;
double r_etas = -9999.;
double r_etac = -9999.;
et_act (
// Parameters and variables needed specifically for et_act()
	max(min(wc_top, paramFun(wc_sat, 1)), paramFun(wc_res, 1)) - paramFun(wc_res, 1),													// Actual volumetric soil water content at topmost horizon (m3/m3)
	max(min(wc_root, wcs_root), wcr_root),					// Actual volumetric soil water content of the root zone (m3/m3)
	wcs_root,									// Volumetric water content at saturation of the root zone (m3/m3)
	wcp_root,									// Volumetric water content of the root zone at permanent wilting point (m3/m3)
	wcr_root,									// Residual volumetric water content of the root zone (m3/m3)
	paramNum(f_etmax)*wcf_root,								// Parameter giving the volumetric water content where et_act equals et_pot, typically wc_etmax / wc_fk = [0.5..0.8] (m3/m3)
	bubble_root*100.,									// Bubbling pressure of the root zone (cm) or (hPa)
	porei_root,								// Pore-size-index of the root zone (-)
	paramNum(wstressmin)*100.,							// Threshold for water stress effect on resistance (begin of stomata closure) (m) OR (100 hPa)
	paramNum(wstressmax)*100.,							// Threshold for water stress effect on resistance (total stomata closure, wilting point) (m) OR (100 hPa)
	paramNum(par_stressHum),						// Parameter to calculate water vapour deficit stomatal conductance stress factor (hPa-1) - in WASA a value of 0.03
// Common meteorological variables
	inputExt(temper),									// Air temperature (°C)
	inputExt(temper_min),								// Minimum temperature within time step (°C)
	inputExt(temper_max),								// Maximum temperature within time step (°C)
  inputExt(glorad),   						// Downward short-wave radiation (W/m2)
	inputExt(rhum),										// relative humidity (%)
	inputExt(wind),										// wind speed (m/s)
	inputExt(apress),       						// Air pressure (hPa)
	inputExt(sundur),									// Sunshine duration of current day (h)
	inputExt(cloud),										// Cloudiness (%)
// Common site-secific parameters
	paramNum(lat),											// Latitude (decimal degree)
	paramNum(lon),											// Longitude of the location of interest (decimal degrees west of Greenwich, e.g. Greenwich: 0°, Berlin: 345°, New York: 75°)
	paramNum(elev),										// Elevation above sea level (m)
// Specific meteorological quantities (commonly calculated internally)
	inputExt(radex),														// Incoming extraterrestrial radiation (i.e. at top of atmosphere) (W/m2)
	inputExt(glorad_max), 							// Downward short-wave radiation under clear (cloudless) sky (Wm-2)
	inputExt(rad_net),										// net incoming (short-wave + long-wave) radiation (Wm-2)
	inputExt(rad_net_soil),									// net incoming (short-wave + long-wave) radiation hitting the soil surface (Wm-2)
	u[INDEX_s_longrad],									// Net incoming long-wave radiation (Wm-2)
	inputExt(soilheat),								// Soil heat flux (Wm-2)
	inputExt(totalheat),								// heat conduction into soil AND plants due to physical and biochemical energy storage (Wm-2)
// Meteorological parameters
	sharedParamNum(h_tempMeas),							// height of temperature measurement above ground (m)
	sharedParamNum(h_humMeas),								// height of humidity measurement above ground (psychrometer) (m)
	sharedParamNum(h_windMeas),							// height of windspeed measurement above ground (m)
	sharedParamNum(emis_a),									// Coefficient a for calculating net emissivity (-)
	sharedParamNum(emis_b),									// Coefficient b for calculating net emissivity (-)
	sharedParamNum(fcorr_a),									// Coefficient a for calculating cloudiness correction factor (-)
	sharedParamNum(fcorr_b),									// Coefficient b for calculating cloudiness correction factor (-)
	sharedParamNum(radex_a),									// Angstrom coefficient (share of radex on glorad under clouds) (-)
	sharedParamNum(radex_b),									// Angstrom coefficient (radex_a+radex_b = share of radex on glorad under clear sky) (-)
	sharedParamNum(f_day),										// soil heat flux calculation: Fraction of net_rad over daytime (in case of sub-daily application) (-)
	sharedParamNum(f_night),									// soil heat flux calculation: Fraction of net_rad over nighttime (in case of sub-daily application) (-)
// Vegetation and land-cover parameters and variables
	paramNum(crop_makk),     					// Crop-factor after Makkink (-)
	paramNum(crop_faoref),     				// Crop-factor for FAO reference approach (-)
	inputExt(cano_height),							// canopy height (m)
	inputExt(lai),											// leaf area index (m2/m2)
	inputExt(alb),											// albedo (-)
	sharedParamNum(ext),											// Canopy extinction coefficient, Beer's law (-) -  in original WASA model code set to 0.5
	paramNum(res_leaf_min),						// Plant-specific minimum (i.e. no stress occurs) stomatal resistance of a single leaf (sm-1)
	paramNum(soil_dens),								// bulk density of soil of topmost soil horizon (kg/m3)
	paramNum(glo_half),								// Solar radiation at which stomatal conductance is half of its maximum (W/m2) - in WASA model a value of 100
	sharedParamNum(res_b),										// Mean boundary layer resistance (sm-1) - in SW and WASA a value of 25
	sharedParamNum(drag_coef),								// Effective value of mean drag coef. of vegetative elements (-) - in SG and WASA a value of 0.07
	sharedParamNum(rough_bare),							// Roughness length of bare substrate (m) - in SW, SG and WASA a value of 0.01
	sharedParamNum(eddy_decay),							// Eddy diffusivity decay constant (-) - by SW and SG set to 2.5 for agricultural crops
	sharedParamNum(rss_a),										// Empirical constant for calculation of soil surface resistance (-) - in original code set to 26
	sharedParamNum(rss_b),										// Empirical constant for calculation of soil surface resistance (-) - in original code set to -1
// Computational parameters
	inputExt(doy),											// Current day of year
	inputExt(hour),												// hour of day in local time (including daylight daving time), range of [0..23]
	inputExt(utc_add),										// Deviation of local time zone from UTC; may vary over the year due to daylight saving time; range of [-12..14] (hours)
	naval,														// NA value
	delta_t,						// time step length (s)
// Choice flags
	sharedParamNum(choice_et),											// flag which method to use
	sharedParamNum(choice_rcs),											// Flag: canopy stomatal resistance by up-scaling of single leaf resistance, 1: Shuttlewort & Wallace (1985) eq. 19, 2: Saugier and Katerji (1991) eq. 4 used in WASA
	sharedParamNum(choice_roughLen), 								// Flag: Roughness length for momentum transfer, 1: SWAT manual (2011) eqs. 2:2.2.4, 2:2.2.5, 2: Shuttleworth & Gurney (1990) eq. 43
	sharedParamNum(choice_plantDispl),							// Flag: Displacement height for a plant, 1: SWAT manual (2011) eq. 2:2.2.7, 2: Shuttleworth & Gurney (1990) eq. 42
	sharedParamNum(choice_gloradmax),								// Flag: calculation of maximum incoming short-wave radiation (clear sky), 1: Angstroem, 2: Allen (2005), ASCE standard etp, eq. 19 (based on elevation)
// output
	r_eta,												// total actual evapotranspiration (m/s)
	r_etas,												// actual evapotranspiration from soil surface; SW approach (m/s)
	r_etac												// actual evapotranspiration from canopy; SW approach (m/s)
);

// SW approach was selected, divide eta into flux from bare soil and vegetation
if ( abs(sharedParamNum(choice_et) - 13.) < 0.01 ) {
	// limit to water content available for evaporation
	// TODO: shouldn't that implicitly be done by scaling of resistance parameters? But still more water is substracted from the soil than is available. Maybe a numerical (or resolution) problem?!
	double wc_et = max(min(wc_top, paramFun(wc_sat, 1)), paramFun(wc_res, 1)) - paramFun(wc_res, 1);
	r_etas = min(r_etas, wc_et * paramFun(hor_depth,1)/delta_t);
	//r_etas = 0.;
	// update soil moisture of uppermost horizon
	flows[0] -= r_etas;
	// limit to actual usable field capacity of root zone
	// TODO: shouldn't that implicitly be done by scaling of resistance parameters? But still more water is substracted from the soil than is available. Maybe a numerical (or resolution) problem?!
	double wc_nfc = min(wc_root, wcf_root);
	wc_nfc = max(wc_nfc - wcp_root, 0.);
	r_etac = min(r_etac, wc_nfc*rootdepth/delta_t);
	// treat r_etac like r_eta for the other approaches
	r_eta = r_etac;
}

// first evaporates interception storage; may limit actual et
r_eta = min(r_eta, r_etp - r_eti);

// distribute flux of actual et to rooted horizons, according to actual plant available water in each horizon
for (unsigned int i=0; i<nh; i++)
	flows[i] -= r_eta * w_eta[i];


//////////////////////////////////////////////////////////
// SURFACE RUNOFF (m/s)
//////////////////////////////////////////////////////////
// calculate averages of wc and wc_sat
double wc_sat_av = 0.;
double wc_a_av = 0.;
for (unsigned int i=0; i<nh; i++) {
	wc_sat_av += paramFun(wc_sat,i+1) * w_depth[i];
	wc_a_av += u[ns+i] * w_depth[i];
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

// calculate infiltration (m/s)
double inf = infiltration(
// General input and parameters
	sharedParamNum(choice_inf),			// Choice flag for selection of method
	in,													// Water hitting top of soil to be infiltrated (ms-1)
	max(min(u[ns], paramFun(wc_sat,1)), paramFun(wc_res,1)),											// Actual volumetric water content at start of infiltration (m3/m3)
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
// SUB-SURFACE RUNOFF (m/s)
//////////////////////////////////////////////////////////

// lateral flows for every horizon
double perc = -9999.;
double latfl = -9999.;
//double r_sum_flows = -9999.;
//double r_water_free = -9999.;
double hor_depth_next = -9999.;
for (unsigned int i=0; i<nh; i++) {
	
	// percolation from ith horizon (m/s)
	if(i == (nh-1)) // always allow percolation into vadose zone under soil profile
		hor_depth_next = 0.;
	else 
		hor_depth_next = paramFun(hor_depth,i+2);
	perc = percolation(sharedParamNum(choice_perc), paramFun(hor_depth,i+1), ku[i], u[ns+i],  paramFun(wc_fc,i+1), delta_t,
										 ku[i+1], mat_pot[i], mat_pot[i+1], hor_depth_next);
	
	if( (abs(perc - naval) < 0.01) || !isfinite(perc) || std::isnan(perc)) {
		stringstream errmsg;
		errmsg << "Calculated percolation is NA or not finite!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}

	
	
	// lateral outflow from ith horizon (m/s)
	latfl = latflow(paramFun(hor_depth,i+1), u[ns+i], paramFun(wc_fc,i+1), paramFun(wc_sat,i+1),
										ksat_scale[i], paramNum(slopelength), paramNum(slope));
	
	if( (abs(latfl - naval) < 0.01) || !isfinite(latfl) || std::isnan(latfl)) {
		stringstream errmsg;
		errmsg << "Calculated lateral flow is NA or not finite!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	
	// check outflows (latflow and percolation); only saturation exccess (wc-wc_fc) can flow out due to gravitation
	// based on Guentner (2002) eqs. 4.47 and 4.48; ATTENTION: the latter is erroneous (but correct in WASA code)!
// 	r_sum_flows = (latfl + perc) * delta_t; // estimated total outflow
// 	if( r_sum_flows > 0. ) {
// 		// maximum possible outflow
// 		r_water_free = paramFun(hor_depth,i+1) * (u[ns+i] - paramFun(wc_fc,i+1));
// 		// if estimated outflow greater than possible reduce outflows accordingly
// 		if(r_sum_flows > r_water_free) {
// 			perc *= r_water_free / r_sum_flows;
// 			latfl *= r_water_free / r_sum_flows;
// 		}
// 	}
	
	
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




//////////////////////////////////////////////////////////
// CHECK FLOWS
//////////////////////////////////////////////////////////
// if( (abs(surf_sat - naval) < 0.01) || !isfinite(surf_sat) || (surf_sat < -1e-10) ) {
// 	stringstream errmsg;
// 	errmsg << "Problem with calculated saturation excess surface runoff: surf_sat = " << surf_sat;
// 	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 	throw(e); 
// }
// if( (abs(surf_inf - naval) < 0.01) || !isfinite(surf_inf) || (surf_inf < -1e-10) ) {
// 	stringstream errmsg;
// 	errmsg << "Problem with calculated infiltration excess surface runoff: surf_inf = " << surf_inf;
// 	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 	throw(e); 
// }
// if( (abs(subsurf - naval) < 0.01) || !isfinite(subsurf) || (subsurf < -1e-10) ) {
// 	stringstream errmsg;
// 	errmsg << "Problem with calculated subsurface runoff: subsurf = " << subsurf;
// 	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 	throw(e); 
// }
// if( (abs(gw_rch - naval) < 0.01) || !isfinite(gw_rch) || (gw_rch < -1e-10) ) {
// 	stringstream errmsg;
// 	errmsg << "Problem with calculated groundwater recharge: gw_rch = " << gw_rch;
// 	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 	throw(e); 
// }
// if( (abs(r_eta - naval) < 0.01) || !isfinite(r_eta) || (r_eta < -1e-10) ) {
// 	stringstream errmsg;
// 	errmsg << "Problem with calculated actual evapotranspiration: r_eta = " << r_eta;
// 	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 	throw(e); 
// }
// if( (abs(r_etp - naval) < 0.01) || !isfinite(r_etp) || (r_etp < -1e-10) ) {
// 	stringstream errmsg;
// 	errmsg << "Problem with calculated potential evapotranspiration: r_etp = " << r_etp;
// 	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 	throw(e); 
// }
// if( (abs(r_eti - naval) < 0.01) || !isfinite(r_eti) || (r_eti < -1e-10) ) {
// 	stringstream errmsg;
// 	errmsg << "Problem with calculated interception evaporation: r_eti = " << r_eti << ", v_interc = " << u[INDEX_v_interc] << ", r_inter = " << r_inter << ", r_etp = " << r_etp << ", r_eta = " << r_eta;
// 	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 	throw(e); 
// }
// 
// for (unsigned int i=0; i<nh; i++) {
// 	if( (abs(flows[i] - naval) < 0.01) || !isfinite(flows[i])) {
// 		stringstream errmsg;
// 		errmsg << "Calculated flow is NA or not finite!";
// 		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 		throw(e); 
// 	}
// }


//////////////////////////////////////////////////////////
// SAVE DERIVATIVES
//////////////////////////////////////////////////////////
// runoff storages
dudt[INDEX_runst_surf_sat] = surf_sat;
dudt[INDEX_runst_surf_inf] = surf_inf;
dudt[INDEX_runst_surf] = surf_sat+surf_inf;
dudt[INDEX_runst_sub] = subsurf;
dudt[INDEX_runst_gw] = gw_rch;
// et storages
dudt[INDEX_v_interc] = r_inter - r_eti;
if(r_etas < -99.)
	dudt[INDEX_et_a] = r_eta;
else
	dudt[INDEX_et_a] = r_eta + r_etas;
dudt[INDEX_et_p] = r_etp;
dudt[INDEX_et_i] = r_eti;
dudt[INDEX_r_interc] = r_inter;
// soil moisture
for (unsigned int i=0; i<nh; i++)
	dudt[ns+i] = flows[i] / paramFun(hor_depth,i+1);
// dummy
dudt[INDEX_s_longrad] = 0.;
