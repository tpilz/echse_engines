#ifndef ET_APPROACHES_H
#define ET_APPROACHES_H


// include libraries (header file) of meteorological functions and constants needed within some functions
#include "meteo/meteo.h"
#include "meteo/meteo_const.h" 


///////////////////////////////////////////////////////////////////////////////
// Potential reference crop ET after Makkink corrected by crop factor.
// Returns: Rate of potential evapotranspiration (m/s)
///////////////////////////////////////////////////////////////////////////////

inline double et_pot_makkink (
  const double glorad,       // Downward short-wave radiation (W/m2)
  const double temper,       // Air temperature (°C)
  const double apress,       // Air pressure (hPa)
  const double cropfactor    // Crop-factor after Makkink (-)
) {
  double s= slopeSatVapPress(temper);
  return (
    // ET_pot after Makking (m/s)
    cropfactor                                // (-)
    * 0.65                                    // (-) Empirical constant
    * s/(s + psychroConst(temper, apress))    // (-)
    * glorad                                  // W/m2 = J/m2/s
    / (latentHeatEvap(temper) * 1000)         // J/kg
    / 1000                                    // kg/m3
  );
}

// Lake evaporation after Makking
// Returns: Rate of evaporation (m/s)
double lakeEvap_makkink(
  const double t,  // Average temperature (°C)
  const double g   // Short-wave downward radiation (W/m2)
) {
  return( 0.61 * (0.439 + 0.01124 * t) * g /
  1.e+06 / (2501. - 2.375 * t) - 0.012/100./86400. );
}


///////////////////////////////////////////////////////////////////////////////
// Evapotranspiration after Penman-Monteith (m/s). 
// Function works for crops with uniform and closed canopy only (big leaf approach)!
// Theoretically not applicable to forests and bare soil or sparsely vegetated surfaces.
// To obtain POTENTIAL evapotranspiration, surface canopy resistance r_c must be set
// to its minimum value (type-specific minimum of resistance; not zero!). 
// For ACTUAL et calculate resistance with respect to plant stress.
// <tobias.pilz@uni-potsdam.de>, FEB 2016
///////////////////////////////////////////////////////////////////////////////

double et_penmon (
	const double lambda,				// Latent heat of water evaporation (J/kg)
	const double delta,					// slope of saturation vapor pressure curve (hPa/K)
	const double H_net,					// net incoming (short-wave + long-wave) radiation (Wm-2)
	const double G,							// soil heat flux (Wm-2)
	const double rho_air,				// air density (kgm-3)
	const double ez_0,					// saturation vapor pressure of air (hPa)
	const double ez,						// vapor pressure of air (hPa)
	const double gamma,					// psychrometric constant (hPa/K)
	const double r_c,						// bulk canopy surface resistance (sm-1)
	const double r_a						// aerodynamic resistance (sm-1)
) {
	// calculate evapotranspiration (mm/s)
	double et = ( delta * (H_net - G) + rho_air * SPHEATMOIST * (ez_0 - ez) / r_a ) / 
							( (delta + gamma * (1. + (r_c/r_a))) * lambda );
							
	// return evapotranspiration in (m/s)
	return(et/1000.);
	
}


///////////////////////////////////////////////////////////////////////////////
// Potential evapotranspiration after Penman-Monteith for FAO reference grass surface (m/s). 
// See http://www.fao.org/docrep/X0490E/x0490e06.htm#chapter 2   fao penman monteith equation
// and chapter 4 for discussion of different temporal resolutions (currently implemented are the
// daily and hourly approach).
//
// Multiply result with a crop factor depending on land-cover to get site-specific et_pot.
//
// Assumptions:
// - surface resistance for a reference grass set to 70 sm-1
// - aerodynamic resistance assumes wind measurement at 2 m for a 0.12 m tall grass
// - albedo of 0.23
// <tobias.pilz@uni-potsdam.de>, SEP 2015, FEB 2016
///////////////////////////////////////////////////////////////////////////////

double et_penmon_ref (
	const double temper,			// Air temperatur (°C)
	double windspeed,					// Windspeed (ms-1)
	double airpress,					// Air pressure (hPa)
	double delta,							// slope of saturation vapor pressure curve (hPa/K)
	double H_net,							// net incoming (short-wave + long-wave) radiation (Wm-2)
	double ez_0,							// saturation vapor pressure of air (hPa)
	double ez,								// vapor pressure of air (hPa)
	const double h_windMeas,	// height of windspeed measurement (m)
	const unsigned int delta_t					// time step length (s)
) {
	
	// calculate windspeed at 2 m reference height if needed
	if(abs(h_windMeas - 2.) < 0.01)
		windspeed = windspeed * (4.87 / log(67.8 * h_windMeas - 5.42));
	
	airpress = airpress/10.; // hPa -> kPa
	
	// saturation vapor pressure of air (kPa)
	ez_0 = ez_0/10.;
	
	// vapor pressure of air (kPa)
	ez = ez/10.;
	
	// psychrometric constant (kPa/K), implicitly includes latent heat of vaporization
	double gamma = 0.665e-03 * airpress;
	
	// slope of saturation vapor pressure curve in (kPa/K)
	delta = delta/10.;
	
	// compute reference evapotranspiration (all assumptions implicitly included) (mm)
	double etp = -9999.;
	if(delta_t == 86400) {
		// net incoming radiation (Wm-2)
		H_net = H_net * 86400. / 1e6; // Wm-2 -> MJ/m2/day
		etp = ( 0.408 * delta * H_net + gamma * 900. / (temper+273.) * windspeed * (ez_0 - ez) ) /
								( delta + gamma * (1. + 0.34 * windspeed) );
		etp = etp / 86400.;
		
	} else if(delta_t == 3600) {
		// net incoming radiation (Wm-2)
		H_net = H_net * 3600. / 1e6; // Wm-2 -> MJ/m2/hour
		etp = ( 0.408 * delta * H_net + gamma * 37. / (temper+273.) * windspeed * (ez_0 - ez) ) /
								( delta + gamma * (1. + 0.34 * windspeed) );
		etp = etp / 3600.;
		
	} else {
		stringstream errmsg;
		errmsg << "Calculation of FAO reference evaporation currently only supported for daily (86400s) and hourly (3600s) time steps!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	// return in (m/s)
	return(etp/1000.);
}


///////////////////////////////////////////////////////////////////////////////
// Evapotranspiration after Suttleworth & Wallace (1985) (m/s). 
// Based on Penman-Monteith but relaxation of the big leaf approach: SW model
// can be applied to sparsely vegetated surfaces showing bare soil. Not applicable
// to forests.
// To obtain POTENTIAL evapotranspiration, surface resistances r_cs, r_ss must be set
// to its minimum values (not zero!). 
// For ACTUAL et calculate resistance with respect to plant stress and dry soil.
//
// Transpiration and soil evaporation can be determined separately if desired:
// 1. apply et_sw() to calculate total et, 2. apply vapPressDeficit_canopy()
// to determine D_0, 3. use D_0 and apply et_sw_soil() and/or et_sw_cano();
// do NOT use ez_0-ez as D_0!
// <tobias.pilz@uni-potsdam.de>, FEB 2016
///////////////////////////////////////////////////////////////////////////////

// sub-function for vapour pressure deficit at canopy source height (hPa)
// Suttleworth & Wallace (1985), eq. 8
double vapPressDeficit_canopy(
	const double H_net,						// net incoming (short-wave + long-wave) radiation (at reference/measurement height) (Wm-2)
	const double totalheat,			// heat conduction into soil AND plants (Wm-2)
	const double gamma,						// psychrometric constant (hPa/K)
	const double delta,						// slope of saturation vapor pressure curve (hPa/K)
	const double lambda,					// Latent heat of water evaporation (J/kg)
	const double vapPressDeficit,	// Vapour pressure deficit at reference (i.e. measurement) height (hPa)
	const double r_aa,							// Aerodynamic resistance between canopy height and reference height (sm-1)
	const double et_total,				// Total evapotranspiration (soil + canopy) (m/s)
	const double rho_air					// Actual (moist) air density (kgm-3)
) {
	// evapotranspiration flux from m/s to mm/s (= kg/m2/s)
	double et = et_total * 1000.;
	
	return( vapPressDeficit + (delta * (H_net-totalheat) - (delta + gamma) * lambda * et) * r_aa / (rho_air * SPHEATMOIST) );
}

// sub-function for soil evaporation (m/s)
// Suttleworth & Wallace (1985), eq. 9
double et_sw_soil (
	const double lambda,				// Latent heat of water evaporation (J/kg)
	const double delta,					// slope of saturation vapor pressure curve (hPa/K)
	const double H_soil,				// net incoming (short-wave + long-wave) radiation hitting the soil surface (Wm-2)
	const double soilheat, 		// soil heat flux (Wm-2)
	const double rho_air,				// air density (kgm-3)
	const double D_0,						// vapour pressure deficit at canopy source height (hPa)
	const double gamma,					// psychrometric constant (hPa/K)
	const double r_ss,					// soil surface resistance (sm-1)
	const double r_sa						// aerodynamic resistance between soil surface and canopy source height (sm-1)
) {
	// soil evaporation
	double ets = (delta * (H_soil - soilheat) + rho_air * SPHEATMOIST * D_0 / r_sa) / (delta + gamma * (1. + r_ss/r_sa)) / lambda;
	
	return(ets/1000.);
}

// sub-function for canopy evaporation (m/s)
// Suttleworth & Wallace (1985), eq. 10
double et_sw_cano (
	const double lambda,				// Latent heat of water evaporation (J/kg)
	const double delta,					// slope of saturation vapor pressure curve (hPa/K)
	const double H_net,					// net incoming (short-wave + long-wave) radiation (at reference/measurement height) (Wm-2)
	const double H_soil,				// net incoming (short-wave + long-wave) radiation hitting the soil surface (Wm-2)
	const double totalheat,		// heat conduction into soil AND plants (Wm-2)
	const double soilheat, 		// soil heat flux (Wm-2)
	const double rho_air,				// air density (kgm-3)
	const double D_0,						// vapour pressure deficit at canopy source height (hPa)
	const double gamma,					// psychrometric constant (hPa/K)
	const double r_cs,					// Bulk stomatal resistance of the canopy (sm-1)
	const double r_ca						// Bulk boundary layer resistance of the vegetative elements in the canopy (sm-1)
) {
	// Total available energy (Wm-2), Suttleworth & Wallace (1985), eq. 3
	double A_total = H_net - totalheat;
	// Energy available at substrat (Wm-2), Suttleworth & Wallace (1985), eq. 5
	double A_s = H_soil - soilheat;
	
	// canopy evaporation
	double etc = (delta * (A_total - A_s) + rho_air * SPHEATMOIST * D_0 / r_ca) / (delta + gamma * (1. + r_cs/r_ca)) / lambda;
	
	return(etc/1000.)
}

// sub-function total evapotranspiration (m/s)
// Suttleworth & Wallace (1985), eqs. 11-18
double et_sw (
	const double lambda,				// Latent heat of water evaporation (J/kg)
	const double delta,					// slope of saturation vapor pressure curve (hPa/K)
	const double H_net,					// net incoming (short-wave + long-wave) radiation (at reference/measurement height) (Wm-2)
	const double H_soil,				// net incoming (short-wave + long-wave) radiation hitting the soil surface (Wm-2)
	const double totalheat,		// heat conduction into soil AND plants (Wm-2)
	const double soilheat, 		// soil heat flux (Wm-2)
	const double rho_air,				// air density (kgm-3)
	const double ez_0,					// saturation vapor pressure of air (hPa)
	const double ez,						// vapor pressure of air (hPa)
	const double gamma,					// psychrometric constant (hPa/K)
	const double r_cs,					// Bulk stomatal resistance of the canopy (sm-1)
	const double r_ca,					// Bulk boundary layer resistance of the vegetative elements in the canopy (sm-1)
	const double r_ss,					// soil surface resistance (sm-1)
	const double r_sa,					// aerodynamic resistance between soil surface and canopy source height (sm-1)
	const double r_aa						// Aerodynamic resistance between canopy source height and reference/measurement height (sm-1)
) {
	// Total available energy (Wm-2), Suttleworth & Wallace (1985), eq. 3
	double A_total = H_net - totalheat;
	// Energy available at substrat (Wm-2), Suttleworth & Wallace (1985), eq. 5
	double A_s = H_soil - soilheat;
	
	// calculate vapor pressure deficit at reference/measurement height (hPa)
	double D = ez_0 - ez;
	
	// calculate term of canopy transpiration (SW eq. 12) (Wm-2)
	double PM_c = ( delta * A_total + ( (rho_air * SPHEATMOIST * D) - (delta * r_ca * A_s) ) / (r_aa + r_ca) ) /
								( delta + ( gamma * (1. + r_cs / (r_aa + r_ca)) ) );
	
	// calculate term of soil evaporation (SW eq. 13) (Wm-2)
	double PM_s = ( delta * A_total + ( (rho_air * SPHEATMOIST * D) - (delta * r_sa * (A_total-A_s)) ) / (r_aa + r_sa) ) /
								( delta + ( gamma * (1. + r_ss / (r_aa + r_sa)) ) );
	
	// SW eqs. 16-18
	double R_a = (delta + gamma) * r_aa;
	double R_s = (delta + gamma) * r_sa + gamma * r_ss;
	double R_c = (delta + gamma) * r_ca + gamma * r_cs;
	
	// coefficients, SW eqs. 14, 15 (-)
	double C_c = 1. / ( 1. + R_c * R_a / (R_s * (R_c + R_a)) );
	double C_s = 1. / ( 1. + R_s * R_a / (R_c * (R_s + R_a)) );
	
	// compute evapotranspiration rate, eq. 11 (mm/s)
	double et = (C_c * PM_c + C_s * PM_s) / lambda;
	
	return(et/1000.)
}

#endif
