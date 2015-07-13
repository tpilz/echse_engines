// Switch to turn balance checks on (CHECK_BALANCE != 0) and off (CHECK_BALANCE == 0)
// NOTE: Effect on computation time appears to be negligible (tested on 2013-01-14)
#define CHECK_BALANCE 1
#define MAX_BALANCE_ERROR 1.e-06  // max. tolerated error in meters / time step

// Other constants
const double SWE_SMALL= 0.001;
const double ZERO= 0.0;
const double ZERO_FLOW= 1.0e-09;

////////////////////////////////////////////////////////////////////////////
// Compute areas (m2)
////////////////////////////////////////////////////////////////////////////

double area_water= paramNum(area) * paramNum(frac_water);
double area_noinf= paramNum(area) * paramNum(frac_noinf);
double area_soil= max(0., paramNum(area) - area_water - area_noinf);

#if CHECK_BALANCE
  // Compute initial storage (for balance check)
  double storage_ini_resv= stateScal(vol_surf) + stateScal(vol_pref) +
    stateScal(vol_inter) + stateScal(vol_base);
  double storage_ini_soil= stateScal(wc) * area_soil * paramNum(soildepth);
  double storage_ini_snow= stateScal(snow_swe) * (area_soil + area_noinf);
#endif

////////////////////////////////////////////////////////////////////////////
// Compute meteo data from residuals
////////////////////////////////////////////////////////////////////////////
double airtemp= inputExt(temper_resid) + inputExt(temper_slope) *
  paramNum(elev) + inputExt(temper_inter);

double apress= inputExt(apress_resid) + inputExt(apress_slope) *
  paramNum(elev) + inputExt(apress_inter);

double precip= max(0., inputExt(precip_resid) + inputExt(precip_slope) *
  paramNum(elev) + inputExt(precip_inter));

////////////////////////////////////////////////////////////////////////////
// Convert and correct precipitation
////////////////////////////////////////////////////////////////////////////

// Conversion of precip: mm/delta_t --> m/s
// Also correction of wind error for sensor heights 1m (rain) and 10m (wind)
double precip_raw= precip * paramNum(fac_precip) / 1000. / delta_t *
  factPrecipCorr(inputExt(windsp), airtemp, 1., 10., sharedParamNum(heightZeroWind));
// Init precipitation falling on soil (not trapped in the snow cover)
double precip_soil= precip_raw;

////////////////////////////////////////////////////////////////////////////
// Compute melt rate, adjust precipitation and update snow variables
////////////////////////////////////////////////////////////////////////////

double snow_melt= 0.;
double snow_subl= 0.;
// Update snow state only if either a snow cover exists or a new cover is
// building up due to snowfall. This also prevents summer rainfall to be
// trapped in a zero snow cover and circumvents the re-estimation of
// rainfall from the (non reasonable) melt rate.
if ((stateScal(snow_swe) > 0.) | (airtemp <=
     sharedParamNum(snow_tempAir_crit))) {
  // Let precipitation be trapped in snow cover
  precip_soil= 0.;
  // Compute derivatives of snow states
  double ddt_sec;
  double ddt_swe;
  double ddt_alb;
  snowModel_derivs(
    // Forcings
    precip_raw * delta_t * 1000.,     // precipSumMM
    inputExt(glorad) * (1. - min(1.,inputExt(lai)/sharedParamNum(snow_fullShadowLAI))),
    airtemp, apress,
    inputExt(rhumid), inputExt(windsp), inputExt(clness),
    delta_t,                 // precipSeconds
    // Parameters
    sharedParamNum(snow_a0), sharedParamNum(snow_a1),
    sharedParamNum(snow_kSat), sharedParamNum(snow_densDry),
    sharedParamNum(snow_specCapRet), sharedParamNum(snow_emissivityMin),
    sharedParamNum(snow_emissivityMax), sharedParamNum(snow_tempAir_crit),
    sharedParamNum(snow_albedoMin), sharedParamNum(snow_albedoMax),
    sharedParamNum(snow_agingRate_tAirPos), sharedParamNum(snow_agingRate_tAirNeg),
    sharedParamNum(snow_soilDepth), sharedParamNum(snow_soilDens),
    sharedParamNum(snow_soilSpecHeat), sharedParamNum(snow_weightAirTemp),
    // States
    stateScal(snow_sec), stateScal(snow_swe), stateScal(snow_alb),
    // Derivatives
    ddt_sec, ddt_swe, ddt_alb,
    // Rate
    snow_melt, snow_subl
  );
  // Pragmatic correction of negative Euler solution for swe
  if ((ddt_swe < 0.) & (abs(ddt_swe * delta_t) > stateScal(snow_swe))) {
    snow_melt= stateScal(snow_swe) / delta_t + precip_raw;
    snow_subl= 0.;
    set_stateScal(snow_swe)= 0.;
    set_stateScal(snow_sec)= 0.;
    set_stateScal(snow_alb)= sharedParamNum(snow_albedoMax);
  } else {
    set_stateScal(snow_swe)= stateScal(snow_swe) + ddt_swe * delta_t;
    set_stateScal(snow_sec)= stateScal(snow_sec) + ddt_sec * delta_t;
    set_stateScal(snow_alb)= stateScal(snow_alb) + ddt_alb * delta_t;
  }
  #if CHECK_BALANCE
    // Check snow balance
    double balance = storage_ini_snow / (area_soil + area_noinf)
      + precip_raw * delta_t - stateScal(snow_swe)
      - snow_melt * delta_t - snow_subl * delta_t;
    if (abs(balance) > MAX_BALANCE_ERROR) {
      stringstream errmsg;
      errmsg << "Error in snow balance of " << scientific << setprecision(5) <<
        balance*1000. << " mm.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);
    }
  #endif
}

////////////////////////////////////////////////////////////////////////////
// Evapotranspiration (m/s)
////////////////////////////////////////////////////////////////////////////

double et_pot;
double et_real;
if (stateScal(snow_swe) < SWE_SMALL) {
  // Potential rate
  et_pot= et_pot_makkink(inputExt(glorad), airtemp, apress,
    sharedParamNum(cropFac_slope) * inputExt(lai) + sharedParamNum(cropFac_inter));
  // Actual rate depending on the soil's relative saturation
  et_real= et_pot * min(1., max(0.,
     (stateScal(wc)/paramNum(wc_max) - paramNum(relsat_etmin)) / 
     (paramNum(relsat_etmax)-paramNum(relsat_etmin)) ));
  // Avoid potential soil underfilling due to 1st order solution
  et_real= min(et_real,
    precip_soil + snow_melt + stateScal(wc) * paramNum(soildepth) / delta_t);
  if (et_real < 1.e-12) // Avoids precision problems in 'runoff_4comp'; threshold is less than 1/1000 mm/day
    et_real= 0.;
} else {
  et_pot= 0.;  // Evaporation already considered in snow balance
  et_real= 0.;
}

////////////////////////////////////////////////////////////////////////////
// Runoff on areas with permeable surface, i.e. soil storage (m/s)
////////////////////////////////////////////////////////////////////////////

double r_surf, r_pref, r_inter, r_base;
runoff_4comp(
  // In
  precip_soil+snow_melt, et_real, stateScal(wc), paramNum(wc_max), paramNum(soildepth),
  paramNum(exp_satfrac), paramNum(thr_surf), paramNum(relsat_inter),
  paramNum(rate_inter), paramNum(rate_base), delta_t,
  // Out
  r_surf, r_pref, r_inter, r_base
);

////////////////////////////////////////////////////////////////////////////
// Update soil water content
////////////////////////////////////////////////////////////////////////////

set_stateScal(wc)= max(ZERO, stateScal(wc) + 
  (precip_soil + snow_melt - r_surf - r_pref - r_inter - r_base - et_real) * delta_t / paramNum(soildepth));

////////////////////////////////////////////////////////////////////////////
// New storage volumes in the reservoirs describing runoff concentration (m3)
////////////////////////////////////////////////////////////////////////////

// Compute storage constants (s) from calib. parameters
double k_surf=  paramNum(str_surf)  * paramNum(ct_index);
double k_pref=  paramNum(str_pref)  * paramNum(ct_index);
double k_inter= paramNum(str_inter) * paramNum(ct_index);
double k_base=  paramNum(str_base)  * paramNum(ct_index);

// Inflow rates (m3/s)
double q_surf=
  r_surf * area_soil +                   // Surface runoff from soil
  precip_raw * area_water +              // Precipitation on water
  (precip_soil + snow_melt) * area_noinf; // Rain + melt on impervious areas
double q_pref= r_pref * area_soil;
double q_inter= r_inter * area_soil;
double q_base= r_base * area_soil;

// New volumes in the reservoirs (m3)
double v_surf=  v_new(stateScal(vol_surf), k_surf, delta_t, q_surf);
double v_pref=  v_new(stateScal(vol_pref), k_pref, delta_t, q_pref);
double v_inter= v_new(stateScal(vol_inter), k_inter, delta_t, q_inter);
double v_base=  v_new(stateScal(vol_base), k_base, delta_t, q_base);

////////////////////////////////////////////////////////////////////////////
// Set output variables
////////////////////////////////////////////////////////////////////////////

// Average catchment outflow (m3/s)
set_output(qx_avg)= (q_surf -  (v_surf -  stateScal(vol_surf))  / delta_t) * sharedParamNum(mult_surf) +
                    (q_pref - (v_pref - stateScal(vol_pref)) / delta_t) * sharedParamNum(mult_pref) +
                    (q_inter - (v_inter - stateScal(vol_inter)) / delta_t) * sharedParamNum(mult_inter) +
                    (q_base -  (v_base -  stateScal(vol_base))  / delta_t) * sharedParamNum(mult_base);
double qx_avg_4balance= (q_surf -  (v_surf -  stateScal(vol_surf))  / delta_t) +
                    (q_pref - (v_pref - stateScal(vol_pref)) / delta_t) +
                    (q_inter - (v_inter - stateScal(vol_inter)) / delta_t) +
                    (q_base -  (v_base -  stateScal(vol_base))  / delta_t);

// Correct presicions problems (negligibly small negative average runoff)
if (output(qx_avg) < ZERO) {
  if (abs(output(qx_avg)) <= ZERO_FLOW) {
    set_output(qx_avg)= 0.;
  } else {
    stringstream errmsg;
    errmsg << "Non-negligible negative catchment runoff computed.";
    except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
    throw(e);
  }
}


// Catchment outflow at end of time step (m3/s)
set_output(qx_end)= v_surf/k_surf * sharedParamNum(mult_surf) +
                    v_pref/k_pref * sharedParamNum(mult_pref) +
                    v_inter/k_inter * sharedParamNum(mult_inter) +
                    v_base/k_base * sharedParamNum(mult_base);

// Snow water equivalent (m)
set_output(swe)= stateScal(snow_swe);

// Evapotranspiration (from pervious surfaces or snow cover)
if (stateScal(snow_swe) < SWE_SMALL) {
  set_output(etr)= et_real;
  set_output(etp)= et_pot;
} else {
  set_output(etr)= snow_subl;
  set_output(etp)= snow_subl;
}

////////////////////////////////////////////////////////////////////////////
// Save state of linear reservoirs
////////////////////////////////////////////////////////////////////////////

set_stateScal(vol_surf)= v_surf;
set_stateScal(vol_pref)= v_pref;
set_stateScal(vol_inter)= v_inter;
set_stateScal(vol_base)= v_base;

#if CHECK_BALANCE
  // Check water balance
  double storage_end_resv=
    stateScal(vol_surf) + stateScal(vol_pref) + stateScal(vol_inter) + stateScal(vol_base);
  double storage_end_soil= stateScal(wc) * area_soil * paramNum(soildepth);
  double storage_end_snow= stateScal(snow_swe) * (area_soil + area_noinf);
  double inputs= precip_raw * paramNum(area) * delta_t;
  double outputs= (et_real * area_soil + snow_subl * (area_soil + area_noinf) + qx_avg_4balance) * delta_t;
  double balance= storage_ini_resv + storage_ini_soil + storage_ini_snow +
       inputs - outputs - storage_end_resv - storage_end_soil - storage_end_snow;
  if (abs(balance / paramNum(area)) > MAX_BALANCE_ERROR) {
    stringstream errmsg;
    errmsg << "Error in catchment water balance of " << scientific <<
      setprecision(5) << balance*1000. << " mm.";
    except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
    throw(e);
  }
#endif
