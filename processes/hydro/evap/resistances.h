#ifndef RESISTANCES_H
#define RESISTANCES_H

// include libraries (header file) of meteorological functions and constants needed within some functions
#include "meteo/meteo.h"
#include "meteo/meteo_const.h"

#include <cmath>

// This file contains functions for calculation of resistances needed for calculation 
// of evapotranspiration in equations of Penman-Monteith type

// At bottom you find the resistance functions, above you find additional functions to 
// calculate plant stress factors, roughness length etc.


////////////////////////////////////////////////////////////////////////////////
// Roughness length for momentum transfer (m).
//
// Choice 1:
// SWAT manual (2011) eqs. 2:2.2.4, 2:2.2.5 citing Brutsaert (1975).
// This is what in Shuttleworth & Gurney (1990) is called the 'preferred value'.
//
// Choice 2:
// Shuttleworth & Gurney (1990) eq. 43.
// In contrast to choice 1 this method is superior as it describes the 
// transition between a closed canopy and a bare substrate via LAI.
// Used in original WASA code.
//
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double roughLen_mom (
	const double choice,						// Choose a method
	const double cano_height,				// Mean height of plant canopy (m)
	const double lai,								// Leaf Area Index (m2m-2)
	const double h_plantDispl,			// Displacement height of vegetation after Shuttleworth & Gurney (1990) (m)
	const double drag_coef,					// Effective value of mean drag coef. of vegetative elements (-) - in SG and WASA a value of 0.07
	const double rough_bare					// Roughness length of bare substrate (m) - in SW, SG and WASA a value of 0.01
) {
	
	double res = -9999.;
	
	switch((int)(choice+.5)) {
		case 1 : // SWAT manual (2011) eqs. 2:2.2.4, 2:2.2.5 citing Brutsaert (1975)
		{
			// set to 0.01 for bare soil, i.e. if cano_height is zero (Shuttleworth & Wallace (1985) citing Van Bavel & Hillel (1976))
			if( cano_height < 0.001 )
				res = rough_bare;
			
			// unit conversion
			double hc = cano_height * 100.;
			
			if (cano_height <= 2.0)
				res = 0.123 * hc / 100.;
			else
				res = 0.058 * pow(hc, 1.19) / 100.;
			
			break;
			
		}
		case 2 : // Shuttleworth & Gurney (1990) eq. 43
			
			if( (cano_height - h_plantDispl) < 0. ) {
				stringstream errmsg;
				errmsg << "Cannot calculate roughness length (SG model) as height of canopy is smaller than displacement height of vegetation!";
				except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
				throw(e);
			}
				
			if ( ((drag_coef * lai) >= 0.) && ((drag_coef * lai) < 0.2) )
				res = rough_bare + 0.3 * cano_height * pow(drag_coef * lai, 0.5);
			else if ( ((drag_coef * lai) >= 0.2) && ((drag_coef * lai) < 1.5) )
				res = 0.3 * cano_height * (1. - h_plantDispl/cano_height);
			else {
				stringstream errmsg;
				errmsg << "Cannot calculate roughness length (SG model) as 'lai' * 'drag_coef' out of range of definition!";
				except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
				throw(e);
			}
			break;
			
		default :
			stringstream errmsg;
      errmsg << "Invalid choice to calculate roughness length for momentum transfer! Currently supported is one of {1,2}.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e); 
	}
	
	// check result
	if( res < -999. ) {
		stringstream errmsg;
		errmsg << "Something during the calculation of roughness length went wrong!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e);
	}
	
	return(res);
}


////////////////////////////////////////////////////////////////////////////////
// Roughness length for vapor transfer (m).
// SWAT manual (2011) eq. 2:2.2.6 citing Stricker and Brutsaert (1978).
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double roughLen_vap (
	const double roughlen_mom				// Roughness length for momentum transfer (m)
) {

	return(0.1 * roughlen_mom);
}


////////////////////////////////////////////////////////////////////////////////
// Displacement heigt for a plant (m).
//
// Choice 1:
// SWAT manual (2011) eq. 2:2.2.7 citing Monteith (1981) and Plate (1971).
// This is what in Shuttleworth & Gurney (1990) is called the 'preferred value'.
//
// Choice 2:
// Shuttleworth & Gurney (1990) eq. 42.
// In contrast to choice 1 this method is superior as it describes the 
// transition between a closed canopy and a bare substrate via LAI.
// Used in original WASA code.
//
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double height_plantDisplace (
	const double choice,						// method to use
	const double cano_height,				// Mean height of plant canopy (m)
	const double lai,								// Leaf Area Index (m2m-2)
	const double drag_coef					// Effective value of mean drag coef. of vegetative elements (-) - in SG and WASA a value of 0.07
) {

	double res = -9999.;
	
	switch((int)(choice+.5)) {
		
		case 1 : // SWAT manual (2011) eq. 2:2.2.7
			res = 2./3. * cano_height;
			break;
		
		case 2 : // Shuttleworth & Gurney (1990) eq. 42
			res = 1.1 * cano_height * log(1 + pow(drag_coef * lai, 0.25));
			break;
	
		default:
			stringstream errmsg;
      errmsg << "Invalid choice to calculate displacement height for a plant! Currently supported is one of {1,2}.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e); 
	}
	
	// check result
	if( res < -999. ) {
		stringstream errmsg;
		errmsg << "Something during the calculation of displacement height went wrong!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e);
	}
	
	return(res);
}


////////////////////////////////////////////////////////////////////////////////
// Actual stomatal resistance of a single leaf under stress (sm-1).
// Model of Jarvis (1976)
//		Stomatal conductance (1/resistance) expressed as multiplicative relation to
//		five environmental factors: quantum flux density, co2 concentration, leaf-air
//		vapour pressure difference, leaf temperature and leaf water status. Each
//		factor is determined independently.
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double res_stom_leaf (
	const double res_min,				// Minimum plant-specific stomatal resistance of a single leaf without any stress (sm-1)
	const double cond_rad,			// Radiation stress factor limiting stomatal conductance [0..1] (-)
	const double cond_co2,			// CO2 stress factor limiting stomatal conductance [0..1] (-)
	const double cond_temp,			// Temperature stress factor limiting stomatal conductance [0..1] (-)
	const double cond_vap,			// Vapour pressure stress factor limiting stomatal conductance [0..1] (-)
	const double cond_water			// Soil water availability stress factor limiting stomatal conductance [0..1] (-)
) {
	//maximum conductance (ms-1)
	double cond_max = 1. / res_min;
	
	// actual conductance (ms-1)
	double cond = cond_max * cond_rad * cond_co2 * cond_temp * cond_vap * cond_water;
	
	// return resistance (sm-1)
	if (cond < 1e-10)
		return( 1.e10 );
	else
		return( 1. / cond );
	
}


////////////////////////////////////////////////////////////////////////////////
// Soil water stress factor for water uptake by plants (-).
// Guentner (2002), eq. 4.29 (limits are erroneous; in Code there is a different solution).
// ATTENTION:
// Reference to Hanan and Prince (1997), eq. 5. But this relation is slightly different 
// to the approach used in WASA.
// Formula for water potential from code, not shown in Guentner (2002)
// but reference to Van Genuchten (1980). However, therein I cannot find the formula.
// Furthermore, due to hysteresis water potential should not be estimated via soil
// water content. Look for alternative approaches.
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double stress_soilwater(
	const double wc,						// Actual volumetric soil water content (m3/m3)
	const double wc_sat,				// Volumetric water content at saturation (m3/m3)
	const double wc_res,				// Residual volumetric water content (m3/m3)
	const double bubble,				// Bubbling pressure (cm) or (hPa)
	const double pores_ind,			// Pore-size-index (-)
	const double wstressmin,		// Threshold for water stress effect on resistance (begin of stomata closure) (cm) or (hPa)
	const double wstressmax			// Threshold for water stress effect on resistance (total stomata closure, wilting point) (cm) or (hPa)
) {
	// check water content values
	if( (wc_res < 0.) || (wc_sat > 1.) || (wc_res > wc_sat) || (wc < wc_res) || (wc > wc_sat) ) {
		stringstream errmsg;
		errmsg << "Cannot calculate soilwater plant stress factor, unreasonable water content values! wc = " << wc << ", wc_sat = " << wc_sat << ", wc_res = " << wc_res << ".";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e);
	}
	
	// relative saturation (-)
	double sat_rel = (wc - wc_res) / (wc_sat - wc_res);
	
	// van-Genuchten parameter m from pore size index (-)
	double m = pores_ind / (pores_ind + 1.);	
	
	// water potential
	double suction = -9999.;
	if(sat_rel >= 0.999)
		suction = 0.;
	else
		suction = (pow(1./pow(sat_rel, 1./m) -1., 1./(pores_ind+1.))) * bubble; 
	
	// soil water stress factor
	if(suction < wstressmin)	// no water stress
		return(1.);
	else if(suction >= wstressmax)	// maximum water stress
		return(0.01); // minimum stress factor (at maximum water stress) as defined in WASA
	else	// between begin of and total stomatal closure
		return( 1. - (suction - wstressmin) / (wstressmax - wstressmin) );
}


////////////////////////////////////////////////////////////////////////////////
// Vapour pressure deficit stress factor for stomatal conductance of plants (-).
// Guentner (2002), eq. 4.29 (limits are erroneous; in Code there is a different solution).
// Hanan and Prince (1997) eq. 4 based on Jarvis (1976), used in WASA (Guentner, 2002, eq. 4.30).
// Empirical parameter:
//		0.26 kPa-1 for Picea sitchensis (Jarvis, 1976)
//		0.12 kPa-1 for Pseudotsuga menziesii (Jarvis, 1976)
//		0.03 hPa-1 WASA hard-coded, Guentner (2002)
//		0.01 to 0.07 mbar-1 for savanna vegetation in Sahel region (Hanan and Prince, 1997)
//		0.06 kg/g (uses specific humidity deficit) for Pinus forest in UK
// ATTENTION: Studies mentioned above use water vapour deficit in different units
// and slightly different equations (1-coef*vap_deficit vs. 1/(1+coef*vap_deficit).
// This affects the parameter values. Herein unit hPa is used. Convert the parameter unit
// during pre-processing if needed.
// Deficit within canopy is needed but usually measurements are available for 
// above canopy only which should be considered!
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double stress_humidity(
	const double vap_deficit,				// Water vapour pressure deficit within canopy (hPa)
	const double par								// Empirical parameter, see notes above.
) {
	return( 1. / (1. + par * vap_deficit) );
}





/////////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////////////
// Aerodynamic resistance to sensible heat and vapor transfer (sm-1).
// SWAT manual (2011) eq. 2:2.2.3.
// For original 'big-leaf' approach of Penman-Monteith.
// Also called "Thom's equation" (Thom, 1975).
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double res_aero (
	const double choice_plantDisplace, // mehod to calculate displacement height for a plant (see function description)
	const double choice_roughLenMom, // method to calculate roughness length (see function description)
	const double h_windMeas,			// height of windspeed measurement (m)
	const double h_humMeas,				// height of humidity measurement (psychrometer) (m)
	const double h_tempMeas,			// height of temperature measurement (m)
	const double windspeed,				// Windspeed (ms-1)
	const double cano_height,			// Mean height of vegetation canopy (m)
	const double rough_bare,			// Roughness length of bare substrate (m) - in SW, SG and WASA a value of 0.01
	const double lai,							// Leaf area index (m2/m2)
	const double drag_coef				// Effective value of mean drag coef. of vegetative elements (-) - in SG and WASA a value of 0.07
) {
	// conversions
	double z_w = h_windMeas;
	double z_p = (h_tempMeas + h_humMeas) / 2.;
	
	// Displacement heigt for a plant (m)
	double d = height_plantDisplace(choice_plantDisplace, cano_height, lai, drag_coef);
	
	// Roughness length for momentum transfer (m)
	double z_om = roughLen_mom(choice_roughLenMom, cano_height, lai, d, drag_coef, rough_bare);
	
	// Roughness length for vapor transfer (m)
	double z_ov = roughLen_vap(z_om);
	
	return( (log((z_w - d) / z_om) * log((z_p - d) / z_ov)) / (pow(KARMAN,2) * windspeed) );
}


////////////////////////////////////////////////////////////////////////////////
// Aerodynamic resistance between canopy height and reference height (sm-1), i.e.
// integration of eddy diffusion from (h_plantDispl+rough_len) to h_windMeas.
// For two-layer approach (in contrast to 'big-leaf' of Penman-Monteith).
// Shuttleworth & Wallace (1985) eq. 27
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double res_aa (
	const double windspeed,				// Windspeed at 'height' (ms-1)
	const double h_windMeas,			// height of windspeed measurement (m)
	const double cano_height,			// Mean height of vegetation canopy (m)
	const double h_plantDispl,		// Displacement height of vegetation (m)
	const double rough_len,				// Roughness length of vegetation (m)
	const double eddy_decay				// Eddy diffusivity decay constant (-) - by SW set to 2.5 for agricultural crops
) {
	if( (h_windMeas - cano_height) < 0.) {
		stringstream errmsg;
		errmsg << "Cannot calculate aerodynamic resistance between canopy height and reference height as height of windspeed measurement is smaller than height of vegetation!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e);
	}

	// Shuttleworth & Wallace (1985) eq. 27
	return( ( log((h_windMeas-h_plantDispl) / rough_len) / (pow(KARMAN,2) * windspeed) ) * 
					( ( log((h_windMeas-h_plantDispl) / (cano_height-h_plantDispl)) ) +
						( cano_height / (eddy_decay * (cano_height-h_plantDispl)) ) *
						( exp(eddy_decay * (1. - (h_plantDispl+rough_len) / cano_height)) - 1. )
					)
	);
}


////////////////////////////////////////////////////////////////////////////////
// Bulk boundary layer aerodynamic resistance of vegetative elements in the canopy (sm-1).
// For two-layer approach (in contrast to 'big-leaf' of Penman-Monteith).
// Shuttlewort & Wallace (1985) eq. 20
// Used in original WASA code (res_b = 25 hard-coded).
// TODO: maybe there are newer / other relationships?
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double res_ca (
	double lai,						// Leaf Area Idex (m2m-2)
	const double res_b					// Mean boundary layer resistance (sm-1) - in SW and WASA a value of 25
) {
	// when there is no vegetation (i.e. LAI -> 0) resistance goes to infinity
	if (lai < 1e-3)
		lai = 1e-3;
	return( res_b / (2 * lai) );
}


////////////////////////////////////////////////////////////////////////////////
// Aerodynamic resistance between soil surface and canopy height (sm-1), i.e.
// integration of eddy diffusion from soil surface to (h_plantDispl+rough_len).
// For two-layer approach (in contrast to 'big-leaf' of Penman-Monteith).
// Shuttleworth & Wallace (1985) eq. 26 but taking roughness length of bare soil
// into account (cf. Shuttleworth & Gurney (1985), eq. 45).
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double res_sa (
	const double windspeed,				// Windspeed at 'height' (ms-1)
	const double h_windMeas,			// height of windspeed measurement (m)
	const double cano_height,			// Mean height of vegetation canopy (m)
	const double h_plantDispl,		// Displacement height of vegetation (m)
	const double rough_len,				// Roughness length of vegetation (m)
	const double eddy_decay,			// Eddy diffusivity decay constant (-) - by SW and SG set to 2.5 for agricultural crops
	const double rough_bare				// Roughness length of bare substrate (m) - in SW, SG and WASA a value of 0.01
) {
	if( (h_windMeas - cano_height) < 0.) {
		stringstream errmsg;
		errmsg << "Cannot calculate aerodynamic resistance between soil surface and canopy height as height of windspeed measurement is smaller than height of vegetation!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e);
	}
	
	// Shuttlewort & Wallace (1985) eq. 26
	return( ( log((h_windMeas-h_plantDispl) / rough_len) / (pow(KARMAN,2) * windspeed) ) * 
					( cano_height / (eddy_decay * (cano_height-h_plantDispl)) ) *
					( exp(eddy_decay * (1. - rough_bare / cano_height)) - exp(eddy_decay * (1. - (h_plantDispl+rough_len) / cano_height)) )
	);
}


////////////////////////////////////////////////////////////////////////////////
// Soil surface resistance (sm-1).
// For two-layer approach (in contrast to 'big-leaf' of Penman-Monteith).
// Empirical relationship by Domingo et al. (1999) under clumped shrub vegetation in Spain.
// For bare soil: a = 15.4, b = -0.76
// Under shrub vegetation vegetation: a = 37.5, b = -1.23
// Used in original WASA code. For the empirical coefficients the mean value from
// Domingo et al. (1999) was taken (a = 26, b = -1).
// TODO: maybe there are newer / other relationships?
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double res_ss (
	double soilwat_grav,				// Gravimetric soil water content (kg/kg)
	double const a,										// Empirical constant.
	double const b										// Empirical constant.
) {
	// soil water content may not be exactly zero (res_ss goes to infinity) -> set minimum to 1e-6 (arbitrary)
	if (soilwat_grav < 1e-6)
		soilwat_grav = 1e-6;
		
	return( a * pow(soilwat_grav, b) );
}


////////////////////////////////////////////////////////////////////////////////
// Bulk stomatal resistance of canopy (sm-1).
// For Shuttlewort & Wallace two-layer approach and 'big-leaf' of Penman-Monteith.
// Different methods implemented:
// 1: Shuttlewort & Wallace (1985) eq. 19
//		Integrate single leaf resistance to the whole canopy using a simple relationship
//		with LAI.
//		TODO:
// 		Technically this should be the same as in SWAT but maybe there are differences
// 		in empirical parameter estimations due to different underlying assumptions or just a mistake?!
//		In SWAT (and FAO) it is res_ST / (0.5*LAI) in contrast SW approach: res_ST / (2*LAI)
// 2: Saugier and Katerji (1991) eq. 4, used in WASA (Guentner, 2002 eq. 4.32).
//		Includes light dependence in integration of resistance from leaf to whole canopy
//		in a physically more meaningful way than just using LAI.
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double res_cs (
	const double choice,				// Wich method shall be used?
	double lai,									// Leaf Area Idex (m2m-2)
	const double res_ST,				// Actual stomatal resistance of leaves (sm-1) - in SW a value of 400
	const double ext,						// Canopy extinction coefficient (-)
	const double glorad,				// Downward short-wave radiation (above canopy) (W/m2)
	const double glo_half				// Solar radiation at which stomatal conductance is half of its maximum (W/m2) - in WASA a value of 100
) {
	
	// if res_ST set to zero (e.g. to calculate potential ET) return zero
	if(res_ST < 1e-6)
		return(0.);
	
	switch((int)(choice+.5)) {
		case 1 : // Shuttlewort & Wallace (1985) eq. 19
			
			// when there is no vegetation (i.e. LAI -> 0) resistance goes to infinity
			if (lai < 1e-3)
				lai = 1e-3;
			return( res_ST / (2.0 * lai) );
			
		case 2 : // Saugier and Katerji (1991) eq. 4, Guentner (2002) eq. 4.32
		{
			if ( glo_half < -90. ) {
				stringstream errmsg;
				errmsg << "Solar radiation at which stomatal conductance is half of its maximum is missing!";
				except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
				throw(e); 
			}
			// leaf conductance (ms-1)
			double cond_l = 1. / res_ST;
			
			// canopy conductance
			double cond_cano = cond_l / ext * log( (glo_half + ext * glorad) / (glo_half + ext * glorad * exp(-1. * ext * lai)) );
			
			// return canopy resistance
			return( 1. / cond_cano );
		}
		default :
			stringstream errmsg;
      errmsg << "Invalid choice to calculate stomatal resistance of canopy! Currently supported is one of {1,2}.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e); 
	}
}



#endif
