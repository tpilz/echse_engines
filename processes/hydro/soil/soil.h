
#ifndef SOIL_H
#define SOIL_H

#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// Sorptivity (ms^{-1/2})
// Quantity appearing in Philip's equation (Philip, 1957) with physical meaning:
// "a measure of the capillary uptake or removal of water", i.e. processes of
// adsorption and desorption causing horizontal movement of water in soil.
// Approach for calculation used here already appeared in Green and Ampt (1911)
// but term was only defined by Philip (1957).
// So far only the simple approach by Green-Ampt implemented (as shown in Stewart
// et al., 2013). Further revised approach in Stewart et al. (2013) but much more
// complicated.
// <tobias.pilz@uni-potsdam.de>, FEB 2016
////////////////////////////////////////////////////////////////////////////////
double sorptivity (
	const double wc,						// Actual volumetric water content at start of infiltration (m3/m3)
	const double wc_sat,				// Volumetric water content at saturation (m3/m3)
	const double ksat,					// Saturated hydraulic conductivity (ms-1)
	const double suc						// Suction at wetting front (m)
) {
	// eq. 4 from Stewart et al. (2013) citing Green and Ampt (1911)
	return( sqrt(2. * ksat * (wc_sat - wc) * suc) );
}


////////////////////////////////////////////////////////////////////////////////
// Matric potential / capillary suction (m of water) or (100 hPa)
// Different soil retention models are implemented.
// ATTENTION: Currently hysteresis is NOT accounted for!
// Choices have to be compliant with k_unsat() below!
// Choices:
//  1: Van Genuchten
//  2: Brooks and Corey
//  3: Campbell
// <tobias.pilz@uni-potsdam.de>, MAR 2016
////////////////////////////////////////////////////////////////////////////////
double matric_pot (
	const unsigned int choice,	// Choice flag for selection of method
	double wc,									// Actual volumetric water content at start of infiltration (m3/m3)
	const double wc_sat,				// Volumetric water content at saturation (m3/m3)
	const double wc_res,				// Residual volumetric water content (m3/m3)
	const double pores_ind,			// Pore-size-index (-)
	const double bubble,				// Bubbling pressure (m) or (100 hPa)
	const double na_val					// NA value
) {
	// check water content values
	if( (wc_res < 0.) || (wc_sat > 1.) || (wc_res > wc_sat) ) {
		stringstream errmsg;
		errmsg << "Cannot calculate matric potential, unreasonable soil parameters!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e);
	}
	
	// limit wc to wc_sat/wc_res as even small exceedance lets function crash
	wc = max(wc_res, min(wc, wc_sat));
	
	// initialise output
	double suction = na_val;
	
	
// VAN GENUCHTEN
	if(choice == 1) {
		// relative saturation (-)
		double sat_rel = (wc - wc_res) / (wc_sat - wc_res);
		
		// limit relativ saturation as suction goes to infinity when approaching zero
		sat_rel = max(sat_rel, 0.001);
		
		// van-Genuchten parameter m from pore size index (-)
		double m = pores_ind / (pores_ind + 1.);	
		
		// water potential
		suction = (pow(1./pow(sat_rel, 1./m) -1., 1./(pores_ind+1.))) * bubble;
	}
	
	
// BROOKS AND COREY
	if(choice == 2) {
		// relative saturation (-)
		double sat_rel = (wc - wc_res) / (wc_sat - wc_res);
		
		// limit relativ saturation as suction goes to infinity when approaching zero
		sat_rel = max(sat_rel, 0.001);
		
		// water potential
		suction = bubble / pow(sat_rel, 1./pores_ind);
	}
	
	
// CAMPBELL
	if(choice == 3) {
		// relative saturation (-)
		double sat_rel = wc / wc_sat;
		
		// water potential
		suction = bubble / pow(sat_rel, 1./pores_ind);
	}
	
	
// ERROR if inf_pot_sum still na_val
	if( abs(suction - na_val) < 0.01 ) {
		stringstream errmsg;
		errmsg << "Could not calculate matric potential: Value of choice parameter not supported, must be one of {1,2,3}!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	return(suction);
}


////////////////////////////////////////////////////////////////////////////////
// Unsaturated hydraulic conductivity (m/s)
// Different hydraulic models are implemented.
// ATTENTION: Currently hysteresis is NOT accounted for!
// Choices have to be compliant with matric_pot() above!
// Choices:
//  1: Van Genuchten
//  2: Brooks and Corey
//  3: Campbell
// <tobias.pilz@uni-potsdam.de>, MAR 2016
////////////////////////////////////////////////////////////////////////////////
double k_unsat (
	const unsigned int choice,	// Choice flag for selection of method
	double wc,									// Actual volumetric water content at start of infiltration (m3/m3)
	const double wc_sat,				// Volumetric water content at saturation (m3/m3)
	const double wc_res,				// Residual volumetric water content (m3/m3)
	const double pores_ind,			// Pore-size-index (-)
	const double ksat,					// Saturated hydraulic conductivity (ms-1)
	const double na_val					// NA value
) {
	// check water content values
	if( (wc_res < 0.) || (wc_sat > 1.) || (wc_res > wc_sat) ) {
		stringstream errmsg;
		errmsg << "Cannot calculate matric potential, unreasonable soil parameters!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e);
	}
	
	// set wc to wc_sat as even small exceedance (considered as insignificant during check) results in NaN output
	wc = max(wc_res, min(wc, wc_sat));
	
	// initialise output
	double k_u = na_val;
	
	
// VAN GENUCHTEN
	if(choice == 1) {
		// relative saturation (-)
		double sat_rel = (wc - wc_res) / (wc_sat - wc_res);
		
		// van-Genuchten parameter m from pore size index (-)
		double m = pores_ind / (pores_ind + 1.);	
		
		// unsaturated conductivity (according to van Genuchten), m/s
		k_u = ksat * pow(sat_rel,0.5) * pow(1. - pow(1. - pow(sat_rel, 1./m), m), 2.);
	}
	
	
// BROOKS AND COREY
	if(choice == 2) {
		// relative saturation (-)
		double sat_rel = (wc - wc_res) / (wc_sat - wc_res);
		
		// Brooy and Corey parameter n from pore size index; same as for Campbell (-)
		double n = 3. + 2. / pores_ind;	
		
		// unsaturated conductivity, m/s
		k_u = ksat * pow(sat_rel,n);
	}
	
	
// CAMPBELL
	if(choice == 3) {
		// relative saturation (-)
		double sat_rel = wc / wc_sat;
		
		// Campbell parameter n from pore size index; same as for Brooks and Corey (-)
		double n = 3. + 2. / pores_ind;	
		
		// unsaturated conductivity, m/s
		k_u = ksat * pow(sat_rel,n);
	}
	
	
// ERROR if inf_pot_sum still na_val
	if( abs(k_u - na_val) < 0.01 ) {
		stringstream errmsg;
		errmsg << "Could not calculate hydraulic conductivity: Value of choice parameter not supported, must be one of {1,2,3}!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	return(k_u);
}


////////////////////////////////////////////////////////////////////////////////
// Infiltration flux (ms-1).
// ATTENTION: approach only for homogeneous soil without macro pores and preferential
// flow!
// Choices:
//  1: Horton equation (parameters have to be given).
//  2: Philip equation (parameters given or can be estimated from soil properties)
//  3: Green-Ampt two-stage approach for layered soil (Green and Ampt (1911) and
//     Peschke (1977, 1987) or Mein and Larson (1971, 1973); also implemeted in
//     WASA: Guentner (2002), chap. 4.2.3 and WaSiM-ETH: Schulla (1997) chap. 2.3.8)
// <tobias.pilz@uni-potsdam.de>, AUG 2015, MAR 2016
////////////////////////////////////////////////////////////////////////////////

double infiltration(
// General input and parameters
	const unsigned int choice,	// Choice flag for selection of method
	const double input,					// Water hitting top of soil to be infiltrated (ms-1)
	const double wc,						// Actual volumetric water content at start of infiltration (m3/m3)
	const double wc_sat,				// Volumetric water content at saturation (m3/m3)
	const double ksat,					// Saturated hydraulic conductivity (ms-1)
	const double delta_t,				// time step length (s)
	const double na_val,				// NA value
// Horton-specific parameters
	const double Hort_ini,			// Horton parameter: initial infiltration rate (m/s)
	const double Hort_end,			// Horton parameter: final infiltration rate (m/s)
	const double Hort_k,				// Horton parameter: decay constant (1/s)
// Philip specific parameters
	double Phil_s,							// Philip parameter: Sorptivity (ms^{-1/2}); calculated internally if set to NA
	double Phil_a,							// Philip parameter: second term parameter (m/s); calculated internally if set to NA
	double Phil_cal,						// Philip parameter: calibration paremter for Phil_a if this is NA (fraction of ksat), should be within [0.2..1.0]
// Green-Ampt after Peschke specific parameters
	const double suc						// Suction at wetting front (m)
) {
	// catch significantly unrealistic input
// 	if (wc > (wc_sat + 1e-6) ) {
// 		stringstream errmsg;
// 		errmsg << "Could not compute infiltration, actual water content is greater than saturated water content (wc = " << wc << ", wc_sat = " << wc_sat << ").";
// 		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 		throw(e); 
// 	}
	
	// calculate refillable porosity
	double na = max(wc_sat - wc, 0.);
	
	// if soil is saturated nothing can infiltrate
	if( na < 0.001)
		return(0.);
	
	// if incoming water flux is lower than ksat everything infiltrates at input rate
	if( input <= ksat )
		return(input);
	
	
	
// POTENTIAL infiltration sum (m)
	double inf_pot_sum = na_val;
	
	
	// HORTON
	if(choice == 1) {
		// check Horton parameters
		if( (abs(Hort_ini - na_val) < 0.01) || (abs(Hort_end - na_val) < 0.01) || (abs(Hort_k - na_val) < 0.01)) {
			stringstream errmsg;
			errmsg << "Could not compute infiltration after Horton as not all Horton-specific parameters are given!";
			except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
			throw(e);
		}
		// calculate potential infiltration sum (m)
		inf_pot_sum = Hort_end * delta_t + 1./Hort_k * (Hort_ini - Hort_end) * (1. - exp(-1. * Hort_k * delta_t) );
	}
	
	
	// PHILIP
	if(choice == 2) {
		// calculate parameters if not available
		if(abs(Phil_a - na_val) < 0.01) {
			Phil_a = Phil_cal * ksat;
		}
		if(abs(Phil_s - na_val) < 0.01) {
			Phil_s = sorptivity(wc,wc_sat,ksat,suc);
		}
		
		// Philip equation
		inf_pot_sum = Phil_s * pow(delta_t, 0.5) + Phil_a * delta_t;
	}
	
	
	// GREEN-AMPT after PESCHKE (two-stage approach for homogeneous layered soil)
	if(choice == 3) {
		// calculate depth of wetting front at saturation (m)
		double d_wet = suc / (input/ksat - 1.); //
		
		// calculate potentially infiltrated water at time of saturation (m)
		double vol_inf_sat = d_wet * na; //
		
		// calculate time until saturation (s)
		double t_sat = vol_inf_sat / input; //
		
		// if t_sat is larger than current time step length soil will not be saturated wihin this time step and all input can infiltrate
		if(t_sat > delta_t) {
			
			inf_pot_sum = input * delta_t;
			
		} else {
			
		// otherwise calculate potential cumulated infiltration (vol_inf_sat + infiltration after saturation till end of time step); find solution interatively (recursive equation) (m)
			double inf_old = vol_inf_sat + input * (delta_t - t_sat) / 2.; // start value: vol_inf_sat + half of input which is still left
			double inf = 0.;
			double err = 1.;
			int n = 1;
			while ( (n<100) && (err>0.0001) ) {
				inf = vol_inf_sat + ksat * (delta_t-t_sat) + na * suc * log( (inf_old + na * suc) / (vol_inf_sat + na * suc) );
				err = abs(inf - inf_old);
				inf_old = inf;
				n = n+1;
			}
			
			inf_pot_sum = inf;
			
		}

	}
	
	
// ERROR if inf_pot_sum still na_val
	if( abs(inf_pot_sum - na_val) < 0.01 ) {
		stringstream errmsg;
		errmsg << "Could not calculate infiltration: Value of choice parameter not supported, must be one of {1,2,3}!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	
// Return ACTUAL infiltration rate as average over time step length (m/s)
	double inf_act = min(input, inf_pot_sum/delta_t);
	if( inf_act < 0. ) {
		stringstream errmsg;
		errmsg << "Calculated infiltration rate is negative!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}

	return( inf_act );
	
}


////////////////////////////////////////////////////////////////////////////////
// Percolation into next horizon or groundwater storage (m/s).
// Calculation based on Guentner (2002), eqs. 4.42, 4.43 citing Arnold et al. (1990)
// and Arnold and Williams (1995) (SWRRB model).
// Makes use of van Genuchten model to calculate unsaturated hydraulic conductivity,
// see (e.g.) Maidment (1993), chap. 5.6.
// <tobias.pilz@uni-potsdam.de>, AUG 2015, MAR 2016
////////////////////////////////////////////////////////////////////////////////

double percolation(
// General input and parameters
	const unsigned int choice,	// Choice flag for selection of method
	const double hor_depth,			// Depth of soil column (m)
	const double ku,						// Unsaturated hydraulic conductivity; can be calculated internally (ms-1)
// SWAT (storage approach) specific
	double wc,									// Actual volumetric soil water content (m3/m3)
	const double wc_fc,					// Volumetric water content at field capacity (for suction 316hPa (pF=2.6) in WASA) (m3/m3)
	const unsigned int delta_t,	// time step length (s)
// Gradient based (Richards' eq.) specific
	const double ku_n,					// Unsaturated hydraulic conductivity of next (deeper) horizon (or lower boundary condition); can be calculated internally (ms-1)
	const double mat_pot,				// Matric potential / capillary suction; can be calculated internally (m of water) or (100 hPa)
	const double mat_pot_n,			// Matric potential / capillary suction of next (deeper) horizon (or lower boundary condition); can be calculated internally (m of water) or (100 hPa)
	const double hor_depth_n		// Depth of soil column of next (deeper) horizon (or half of actual layer as boundary condition) (m)
) {
	// initialise output
	double perc_out = -9999.;
	
	
// SWAT storage approach
	if (choice == 1) {
		// no percolation if water content is below field capacity
		if(wc <= wc_fc)
			return(0.);
		
		// no percolation if underlying layer is close to saturation
// 		if( wc_n > (wc_fc_n + 0.5 * (wc_sat_n - wc_fc_n)) )
// 			return(0.);
		
		// excess water for percolation (m)
		double v_perc = hor_depth * (wc - wc_fc);
		
		// calculate travel time of excess water in horizon (Guentner 2002, eq. 4.43), s
		double ttime = v_perc / ku;


		// calculate percolation (m/s)
		perc_out = v_perc * (1. - exp(-1. * delta_t / ttime)) / delta_t;
	}
	
	
// GRADIENT of matric potential, based on Richards' equation (strongly dependent on soil water retention and concudctivity model)
	if (choice == 2) {
		// gradient depth (distance between midpoints of the two horizons) (m)
		double hor_dist = (hor_depth+hor_depth_n)/2.;
		
		// depth weights (weight conductivity and 
		double w1 = 0.5*hor_depth / hor_dist;
		double w2 = 0.5*hor_depth_n / hor_dist;
		
		// conductivity (weighted mean of current and lower soil layer) (m/s)
		double ku_m = w1*ku + w2*ku_n;
		//double ku_m = 2.*ku*ku_n / (ku + ku_n);
		
		// calculate flow rate (m/s)
		perc_out = ku_m * ( (mat_pot_n - mat_pot) / hor_dist + 1. );
		
// 		if( perc_out < 0. ) {
// 			stringstream errmsg;
// 			errmsg << "Calculated percolation rate is negative! mat_pot = " << mat_pot << ", mat_pot_n = " << mat_pot_n << ", ku_m = " << ku_m;
// 			except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 			throw(e); 
// 		}
	}
	
	
// ERROR if inf_pot_sum still na_val
	if( abs(perc_out + 9999.) < 0.01 ) {
		stringstream errmsg;
		errmsg << "Could not calculate percolation: Value of choice parameter not supported, must be one of {1,2}!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	
	
	return( perc_out );
}


////////////////////////////////////////////////////////////////////////////////
// Lateral subsurface outflow (m/s).
// Calculation based on Guentner (2002), chap. 4.2.4.
// Based on soil moisture state, slope length and slope gradient.
// NOTE: In contrast to descriptions in Guentner (2002) this is a generic implementation,
//			 i.e. slope and slopelength must be given specifically for the present 
//			 soil profile (in WASA these are TC and LU specific and equations include
//			 conversion factors).
// <tobias.pilz@uni-potsdam.de>, AUG 2015, MAR 2016
////////////////////////////////////////////////////////////////////////////////

double latflow(
	const double depth,					// Depth of horizon (m)
	const double wc,						// Actual volumetric soil water content (m3/m3)
	const double wc_fc,					// Volumetric water content at field capacity (for suction 316hPa (pF=2.6) in WASA) (m3/m3)
	const double wc_sat,				// Volumetric water content at saturation (m3/m3)
	const double ksat,					// Saturated hydraulic conductivity (ms-1)
	const double slopelength,		// Average slope length (m)
	const double slope					// Average slope gradient (-)
) {
	// catch significantly unrealistic input
// 	if (wc > (wc_sat + 1e-6) ) {
// 		stringstream errmsg;
// 		errmsg << "Could not compute lateral outflow, actual water content is greater than saturated water content (wc = " << wc << ", wc_sat = " << wc_sat << ").";
// 		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
// 		throw(e); 
// 	}
	
	// no lateral flow if water content is below field capacity
	if(wc <= wc_fc)
		return(0.);
	
	// saturated depth of soil horizon (m), Guentner (2002) eq. 4.46
	double d_sat = depth * min( (wc - wc_fc) / (wc_sat - wc_fc), 1.);
	
	// calculate lateral flow (m/s), Guentner (2002) eqs. 4.44, 4.45 modified
	double latfl_out = slope * ksat * d_sat / slopelength;
	if( latfl_out < 0. ) {
		stringstream errmsg;
		errmsg << "Calculated lateral outflow rate is negative!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	return( latfl_out );
}


////////////////////////////////////////////////////////////////////////////////
// Fraction of saturation of a spatial element with representative soil profile.
// Based on soil distribution function of Guentner (2002) fig. 4.3, WASA code rev. 130
// soildistr.f90 used in soilwat.f90 e.g. lines 519 ff. (ATTENTION: code might be erroneous?).
// Fraction of saturation depends on soil storage capacity whereas values are
// interpolated between 5 node points (cf. fig. 4.3 of Guentner, 2002):
// <tobias.pilz@uni-potsdam.de>, AUG 2015, MAR 2016
////////////////////////////////////////////////////////////////////////////////

double f_saturation(
	const double wc,						// Actual volumetric water content of soil profile (m3/m3)
	const double wc_sat,				// Saturated volumetric water content of soil profile (m3/m3)
	const double var1,					// Fraction of wc_sat (first node, area starts becoming saturated; should be greater than field capacity) (-)
	const double var2,					// Fraction of wc_sat (second node) (-)
	const double var3,					// Fraction of wc_sat when about half the soil is saturated (-)
	const double var4,					// Fraction of wc_sat (fourth node) (-)
	const double var5,					// Fraction of wc_sat when soil is saturated, var5 >= var4 >= var3 >= var2 >= var1 (-)
	const double frac1,					// Area fraction for first node point (-) - in WASA set to 0
	const double frac2,					// Area fraction for second node point (-) - in WASA set to 0.1
	const double frac3,					// Area fraction for third node point (-) - in WASA set to 0.5
	const double frac4,					// Area fraction for fourth node point (-) - in WASA set to 0.9
	const double frac5					// Area fraction for fifth node point (-) - in WASA set to 1
) {
	
	// check parameters
	if( (var2 < var1) || (var3 < var2) || (var4 < var3) || (var5 < var4) ) {
		stringstream errmsg;
		errmsg << "Could not compute fraction of saturation. Check parameters varN (fractions of saturated water content)!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	if( (frac2 < frac1) || (frac3 < frac2) || (frac4 < frac3) || (frac5 < frac4) ) {
		stringstream errmsg;
		errmsg << "Could not compute fraction of saturation. Check parameters fracN (values of saturated areal fraction)!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	
	double node [5];
	node[0] = wc_sat * var1;
	node[1] = wc_sat * var2;
	node[2] = wc_sat * var3;
	node[3] = wc_sat * var4;
	node[4] = wc_sat * var5;
	
	// no saturation at all
	if( wc <= node[0] )
		return(0.);
	
	// completely saturated
	if( wc > node[4] )
		return(1.);
	
	double frac [5] = { frac1, frac2, frac3, frac4, frac5 };
	
	// in all other cases: linear interpolation between node points
	for (int i = 0; i <= 4; i++) {
		if ( (wc > node[i]) && (wc <= node[i+1]) ) // find node points framing wc
			return( frac[i] + (frac[i+1] - frac[i]) * (1. - (node[i+1] - wc) / (node[i+1] - node[i])) );
	}
	
	// all other cases: something must have gone wrong
	stringstream errmsg;
	errmsg << "Could not compute fraction of saturation. Check your input data! wc = " << wc << ", wc_sat = " << wc_sat << ".";
	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
	throw(e); 
}


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

