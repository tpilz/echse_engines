
#ifndef SOIL_H
#define SOIL_H

// Surface runoff from saturated areas
//
// This is Eqn. 3.9 and 3.10 in Bremicker (2000).
//
// Originally, the equations are from the ARNO model but in the paper
// of Todini (1996), some signs are wrong.
//
// Returns: Amount of surface runoff generated in time step (m)
inline double directRunoffHeight_arno(
  const double input, // Amount of water supply in time step (m)
  const double w,     // Current filling of soil reservoir (m)
  const double wmax,  // Max. capacity of soil reservoir (m)
  const double b      // Shape parameter of Xinanjiang approach
) {
  double x= pow((1.-w/wmax),(1./(b+1.))) - input/(1.+b)/wmax;
  if (x > 0.) {
   return(input - (wmax - w) + wmax*pow(x,(b+1.)));
  } else {
   return(input - (wmax - w));
  }
}

////////////////////////////////////////////////////////////////////////////////

#define DEBUGMODE 0

// Calculates the runoff components (as in LARSIM)
void runoff_4comp(
  // In
  const double input,
  const double et_real,
  const double wc,
  const double wc_max,
  const double soildepth,
  const double exp_satfrac,
  const double thr_surf,
  const double relsat_inter,
  const double rate_inter,
  const double rate_base,
  const double delta_t,
  // Out
  double &r_surf,
  double &r_pref,
  double &r_inter,
  double &r_base
) {

  // Fixed parameters (as in LARSIM)
  const double RELSAT_BASE= 0.05;
  const double EXP_INTER= 1.5;
  const double EXP_BASE= 1.;

  // Relative soil saturation
  double relSat= wc / wc_max;

  // Direct runoff generated on saturated areas
  double r_direct= directRunoffHeight_arno(input * delta_t,
    wc*soildepth, wc_max*soildepth, exp_satfrac) / delta_t;
  // Surface runoff
  r_surf= max(0., r_direct - thr_surf);
  // Quick subsurface runoff
  r_pref= r_direct - r_surf;
  // Interflow
  if (relSat > relsat_inter) {
    r_inter= rate_inter * pow( (relSat-relsat_inter) /
      (1.-relsat_inter), EXP_INTER);
  } else {
    r_inter= 0.;
  }
  // Baseflow  --> Same expression as for interflow but with a fixed constant
  if (relSat > RELSAT_BASE) {
    r_base= rate_base * pow( (relSat-RELSAT_BASE) /
      (1.-RELSAT_BASE), EXP_BASE);
  } else {
    r_base= 0.;
  }

  #if DEBUGMODE
    if (!isfinite(r_surf)) cout << "raw r_surf" << r_surf << endl;
    if (!isfinite(r_pref)) cout << "raw r_pref" << r_pref << endl;
    if (!isfinite(r_inter)) cout << "raw r_inter" << r_inter << endl;
    if (!isfinite(r_base)) cout << "raw r_base" << r_base << endl;
  #endif

  // Correction of runoff rates to avoid soil underfilling
  // --> This is necessary since we are using a 1st order solution

  // Change in soil water content using the estimates
  double change= (input - et_real - r_surf - r_pref - r_inter - r_base) * delta_t / soildepth;

  // Case 1 (overfilling to due precision problems)
  if (change > (wc_max - wc)) {
      stringstream errmsg;
      errmsg << "Overflow of soil reservoir. Please report this bug.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);  
  }
  // Case 2 (underfilling due to simple numerical solution) --> Reduce all losses proportionally
  if (change < (-1. * wc)) {
    // f = Ratio of available water over losses due to runoff
    double f= (wc * soildepth / delta_t + input - et_real) /
      (r_surf + r_pref + r_inter + r_base);
    r_surf=  f * r_surf;
    r_pref=  f * r_pref;
    r_inter= f * r_inter;
    r_base=  f * r_base;
  }

  #if DEBUGMODE
    if (!isfinite(r_surf)) cout << "adj r_surf" << r_surf << endl;
    if (!isfinite(r_pref)) cout << "adj r_pref" << r_pref << endl;
    if (!isfinite(r_inter)) cout << "adj r_inter" << r_inter << endl;
    if (!isfinite(r_base)) cout << "adj r_base" << r_base << endl;
  #endif
}

#endif

