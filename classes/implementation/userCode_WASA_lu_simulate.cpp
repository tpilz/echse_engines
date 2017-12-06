
///////////////////////////////////////////////////////////////////////////////
// Runoff concentration (i.e. delay of generated runoff)
// Calculate runoff storages at the end of the time step (m).
///////////////////////////////////////////////////////////////////////////////
double v_surf = -9999.;
double v_inter = -9999.;
double v_base = -9999.;
double r_base = -9999.;

// conceptual linear reservoir storage approach; makes sense only if lateral re-distribution approach by TCs and SVCs is not used (i.e. one SVC per TC per LU)
if ( abs(sharedParamNum(choice_runconc) - 1.) < 0.01) {
	// Compute storage constants (s) from calib. parameters
	double k_surf =  sharedParamNum(str_surf)  * paramNum(ct_index);
	double k_inter = sharedParamNum(str_inter) * paramNum(ct_index);

	// New volumes in the reservoirs (m)
	v_surf =  v_new(stateScal(vol_surf), k_surf, delta_t, inputSim(r_river_surf));
	v_inter = v_new(stateScal(vol_inter), k_inter, delta_t, inputSim(r_river_sub));
	
} else if ( abs(sharedParamNum(choice_runconc) - 2.) < 0.01) {
// more physically-based approach by Guentner (2004) using a hierarchical flow redistribution system (actual runoff concentration already calculated TC level!)
	v_surf = 0.;
	v_inter = 0.;
} else {
	stringstream errmsg;
	errmsg << "Invalid choice to calculate runoff concentration! Currently supported is one of {1,2}.";
	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
	throw(e); 
}


// groundwater:
if ( sharedParamNum(choice_gw) < 0.01) {
	// groundwater disabled (groundwater recharge collected for output but in fact it is "lost" as it is not furtherly treated)
	v_base = stateScal(vol_base) + inputSim(r_gw_rch) * delta_t;
	r_base = 0.;
} else if ( abs(sharedParamNum(choice_gw) - 1.) < 0.01) {
	// simple linear storage approach
	double k_base =  sharedParamNum(str_base)  * paramNum(ct_index);
	v_base =  v_new(stateScal(vol_base), k_base, delta_t, inputSim(r_gw_rch)); // calculate groundwater reservoir volume at the end of time step
	r_base = (inputSim(r_gw_rch) - (v_base - stateScal(vol_base)) / delta_t) * paramNum(frac_area); // groundwater runoff, area-weighted for summation to Subbasin level
} else {
	stringstream errmsg;
	errmsg << "Invalid choice to calculate groundwater flows! Currently supported is one of {<=0,1}.";
	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
	throw(e); 
}




///////////////////////////////////////////////////////////////////////////////
// Output: Runoff (m/s)
// Area-weighted for summation to Subbasin level; for investigation divide by
// areal fraction of LU in Subbasin.
///////////////////////////////////////////////////////////////////////////////
// average outflows
double r_surf = (inputSim(r_river_surf) - (v_surf - stateScal(vol_surf)) / delta_t) * paramNum(frac_area);
double r_inter = (inputSim(r_river_sub) - (v_inter - stateScal(vol_inter)) / delta_t) * paramNum(frac_area);

set_output(r_out_surf) = r_surf;
set_output(r_out_inter) = r_inter;
set_output(r_out_base) = r_base;



///////////////////////////////////////////////////////////////////////////////
// Update runoff storages (m)
///////////////////////////////////////////////////////////////////////////////
set_stateScal(vol_surf) = v_surf;
set_stateScal(vol_inter) = v_inter;
set_stateScal(vol_base) = v_base;



///////////////////////////////////////////////////////////////////////////////
// general output, area-weighted for summation at Subbasin level
///////////////////////////////////////////////////////////////////////////////
set_output(v_plantwat) = inputSim(v_plantwat_tc) * paramNum(frac_area);
set_output(v_soilwat) = inputSim(v_soilwat_tc) * paramNum(frac_area);
set_output(v_runstor) = (v_surf+v_inter+v_base) * paramNum(frac_area);
set_output(run_surf_inf) = inputSim(run_surf_inf_tc) * paramNum(frac_area);
set_output(run_surf_sat) = inputSim(run_surf_sat_tc) * paramNum(frac_area);
set_output(run_sub) = inputSim(run_sub_tc) * paramNum(frac_area);
set_output(run_gw) = inputSim(r_gw_rch) * paramNum(frac_area);
set_output(etp) = inputSim(etp_tc) * paramNum(frac_area);
set_output(eta) = inputSim(eta_tc) * paramNum(frac_area);
set_output(eti) = inputSim(eti_tc) * paramNum(frac_area);
