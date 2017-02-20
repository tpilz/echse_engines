
// TODO: special treatment for impervious surfaces (tc area = sum(svc areas) + impervious surfaces area in pre-processing of LUMP R-Package)

///////////////////////////////////////////////////////////////////////////////
// Initial calcvulations and preparation
///////////////////////////////////////////////////////////////////////////////
// get areal fraction of this TC and the other downslope TCs
double frac_tc = paramFun(pos2area, paramNum(position));
int no_tc_ds = (int)((paramNum(no_tc) - paramNum(position)) + 0.5); // no of downslope TCs as integer
int pos = (int)(paramNum(position)+0.5); // position of this TC as integer
double frac_ds = 0.;
for (int i = 1; i <= no_tc_ds; i = i+1) {
	frac_ds += paramFun(pos2area, i+pos);
}

// total surface flow generated on this TC
double r_latsur_tc = inputSim(r_latsur_inf_tc) + inputSim(r_latsur_sat_tc);
double r_latsur_svc = inputSim(r_latsur_inf_svc) + inputSim(r_latsur_sat_svc);

// re-scale input from upstream areas (as flows are in m/s, i.e area-weighted to upstream TC, they now have to be related to the area of this TC)
double latsur_up = inputSim(r_latsur_up) / frac_tc;
double latsub_up = inputSim(r_latsub_up) / frac_tc;

// general output, area-weighted for summation in TC !!!
set_output(v_plantwat) = inputSim(plantwat_svc) * frac_tc;
set_output(v_soilwat) = inputSim(soilwat_svc) * frac_tc;
set_output(run_surf_inf) = (inputSim(r_latsur_inf_tc) + inputSim(r_latsur_inf_svc)) * frac_tc;
set_output(run_surf_sat) = (inputSim(r_latsur_sat_tc) + inputSim(r_latsur_sat_svc)) * frac_tc;
set_output(run_sub) = (inputSim(r_latsub_tc) + inputSim(r_latsub_svc)) * frac_tc;
set_output(run_gw) = inputSim(run_gw_svc) * frac_tc;
set_output(etp) = inputSim(etp_svc) * frac_tc;
set_output(eta) = inputSim(eta_svc) * frac_tc;
set_output(eti) = inputSim(eti_svc) * frac_tc;



///////////////////////////////////////////////////////////////////////////////
// Compute flow distribution to downslope TCs or to the river (i.e. subbasin 
// outlet) based on relative position along hillslope and areal fractions.
// Cf. Fig. 4.7 Guentner, 2002.
///////////////////////////////////////////////////////////////////////////////

// subsurface flow: completely distribute to next downslope TC or river
// surface flow: distribute according to areal fractions and TC position
if ( (paramNum(no_tc) - paramNum(position)) < 0.01) { // downslope TC
	
	// river outflow
	// TODO: respect riverbed depth for r_latsub_tc
	set_output(r_river_surf) = r_latsur_tc * frac_tc;
	set_output(r_river_sub) = inputSim(r_latsub_tc) * frac_tc;
	
	// no further downslope TC -> set output to zero
	set_output(r_latsur_tc_out) = 0.;
	set_output(r_latsub_tc_out) = 0.;
	
	// set output for re-distribution among SVCs
	set_output(r_latsur_svc_out) = latsur_up + r_latsur_svc;
	set_output(r_latsub_svc_out) = latsub_up + inputSim(r_latsub_svc);
	
} else if (paramNum(position) < 1.5) { // upslope TC
	
	// surface flow
	// river outflow equal to areal fraction of this TC within LU (only surface flow as preferential flow)
	set_output(r_river_surf) = r_latsur_tc * frac_tc * frac_tc;
	set_output(r_river_sub) = 0.;
	
	// rest of surface flow to downslope TCs
	set_output(r_latsur_tc_out) = r_latsur_tc * (1. - frac_tc) * frac_tc;
	
	// subsurface flow goes completely into next TC
	set_output(r_latsub_tc_out) = inputSim(r_latsub_tc) * frac_tc;
	
	// set output for re-distribution among SVCs
	set_output(r_latsur_svc_out) = r_latsur_svc;
	set_output(r_latsub_svc_out) = inputSim(r_latsub_svc);
	
} else { // middle TC
	
	// aggregate subsurface flows to be distributed among SVCs:
	// flow from upstream TC and flow generated from SVCs for re-distribution
	set_output(r_latsub_svc_out) = inputSim(r_latsub_svc) + latsub_up;
	
	// aggregate surface flows to be distributed among SVCs:
	// portion of surface flow from upslope TCs draining into this TC
	double r_latsur_in = latsur_up * frac_tc / (frac_tc+frac_ds);
	// and flow generated from SVCs for re-distribution
	set_output(r_latsur_svc_out) = r_latsur_svc + r_latsur_in;
	
	// surface flows to be directly routed to river (preferential flow):
	// portion of surface flow generated from SVCs for distribution to downstream TCs
	set_output(r_river_surf) = r_latsur_tc * frac_tc / (frac_tc+frac_ds) * frac_tc;
	set_output(r_river_sub) = 0.;
	
	// surface flows to be routed to downstream TCs:
	// portion of surface flow from upslope TCs
	double r_latsur_out = latsur_up * frac_ds / (frac_tc+frac_ds);
	// and portion of surface flow generated from SVCs for distribution to downstream TCs
	double r_latsur_out2 = r_latsur_tc * frac_ds / (frac_tc+frac_ds);
	set_output(r_latsur_tc_out) = (r_latsur_out + r_latsur_out2) * frac_tc;
	
	// subsurface flows to be routed to downstream TCs:
	// subsurface flow generated from SVCs for distribution to downstream TCs
	set_output(r_latsub_tc_out) = inputSim(r_latsub_tc) * frac_tc;
	
}
