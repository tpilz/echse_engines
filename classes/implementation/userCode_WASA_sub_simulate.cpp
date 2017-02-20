
// area in m2
double area_m = paramNum(area)*1.e6;

///////////////////////////////////////////////////////////////////////////////
// Runoff (m3/s)
///////////////////////////////////////////////////////////////////////////////
set_output(r_out_surf) = inputSim(r_river_surf)*area_m;
set_output(r_out_inter) = inputSim(r_river_sub)*area_m;
set_output(r_out_base) = inputSim(r_river_gw)*area_m;
set_output(r_out_total) = (inputSim(r_river_surf) + inputSim(r_river_sub) + inputSim(r_river_gw))*area_m;


///////////////////////////////////////////////////////////////////////////////
// general output
///////////////////////////////////////////////////////////////////////////////
set_output(v_plantwat) = inputSim(v_plantwat_lu);
set_output(v_soilwat) = inputSim(v_soilwat_lu);
set_output(v_runstor) = inputSim(v_runstor_lu);
set_output(run_surf_inf) = inputSim(run_surf_inf_lu);
set_output(run_surf_sat) = inputSim(run_surf_sat_lu);
set_output(run_sub) = inputSim(run_sub_lu);
set_output(run_gw) = inputSim(run_gw_lu);
set_output(etp) = inputSim(etp_lu);
set_output(eta) = inputSim(eta_lu);
set_output(eti) = inputSim(eti_lu);
