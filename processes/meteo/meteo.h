
#ifndef METEO_H
#define METEO_H

#include <cmath>

// include library (header file) of constants needed within some functions
#include "meteo/meteo_const.h"


////////////////////////////////////////////////////////////////////////////////
// Atmospheric pressure (hPa) as function of elevation above sea level.
// Strongly simplified formula assuming: standard pressure (1013.25 hPa) at sea
// level, 20°C air temperature at sea level, and a temperature gradient of 0.65 K per 100 m
// elevation. This formula can be used to roughly estimate atmospheric pressure
// if the output of your target variable (e.g. evapotranspiration) is not very
// sensitive to pressure values.
//
// Source: e.g. http://www.fao.org/docrep/X0490E/x0490e07.htm#atmospheric pressure (p)
// or Wikipedia or textbooks.
////////////////////////////////////////////////////////////////////////////////
double apress_simple (
	const double elev				// elevation above sea level (m)
) {
	return( 1013.25 * pow( 1. - (0.0065 * elev / 293.), 5.255) );
}


////////////////////////////////////////////////////////////////////////////////
// Saturation vapor pressure as a function of temperature
//
// Source: Dyck & Peschke (1995), p. 58, Eq. 4.6
//
// NOTE T. Pilz: "overWater" Magnus eq. gives same results as (tested with R): 
// - SWAT manual (2011) eq. 1:2.3.2: exp( (16.87 * temp - 116.9) / (temp + 237.3) ) * 10.
// - Dyck & Peschke (1995) eq. 11.13: 6.11 * exp(17.3 * temp / (237.3 + temp))
//
// Input:  Temperature (degree celsius)
// Output: Result in (hPa)
double satVapPress_overWater (const double temp) {
 return( 6.11 * pow(10., 7.5*temp/(237.3+temp)) );  // Magnus equ. over water
}
double satVapPress_overIce (const double temp) {
  return( 6.11 * pow(10., 9.5*temp/(265.5+temp)) );  // Magnus equ. over ice
}


////////////////////////////////////////////////////////////////////////////////
// Vapor pressure as a function of temperature and relative humidity
//
// Source: Dyck & Peschke (1995), p. 61, Eq. 4.12 together with the equation
//         to compute the saturation vapor pressure

// Input:  Temperature (degree celsius)
// Input:  Relative humidity (%)
// Output: Result in (hPa)
double vapPress_overWater (const double temp, const double relhum) {
  return( satVapPress_overWater(temp) * relhum / 100. );
}
double vapPress_overIce (const double temp, const double relhum) {
  return( satVapPress_overIce(temp) * relhum / 100. );
}

////////////////////////////////////////////////////////////////////////////////
// Dew point temperature as a function of temperature and relative humidity
//
// Source: See explanation in documentation of the snow model

// Input:  Vapor Pressure (hPa)
// Output: Result in (hPa)
double dewTemp_overWater (const double vapPress) {
  return( 237.3 * log10(vapPress / 6.11) / (7.5 - log10(vapPress / 6.11)) );
}
double dewTemp_overIce (const double vapPress) {
  return( 265.5 * log10(vapPress / 6.11) / (9.5 - log10(vapPress / 6.11)) );
}

////////////////////////////////////////////////////////////////////////////////
// Vapor pressure deficit (E(T) - e) in hPa as a function of temperature T (°C)
// and relative humidity U (%)
//
// Source: Just a combination of the equations for the actual vapor pressure and
//         the saturation vapor pressure (see there).
//           
// Input:  Temperature (degree celsius)
// Input:  Relative humidity (%)
// Output: Result in (hPa)
double vapPressDeficit_overWater (const double temp, const double relhum) {
  return( satVapPress_overWater(temp) * (1. - relhum / 100.) );
}
double vapPressDeficit_overIce (const double temp, const double relhum) {
  return( satVapPress_overIce(temp) * (1. - relhum / 100.) );
}

////////////////////////////////////////////////////////////////////////////////
// Slope of the saturation vapor pressure curve, i.e. derivative of saturation
// vapor pressure with resp. to temperature  s = d E(T) / d T
//
// Source: Dyck & Peschke (1995), p. 188, Eq. 11.13
//
// TODO T. Pilz: In Dyck & Peschke (1995) factor 4098 instead of 4105.3 is used
// -> de-activate this formula util reason for using 4105.3 is known! (deviations
//    are negligible)
//
// Input:  Temperature (degree celsius)
// Output: Result in (hPa/K)
// double slopeSatVapPress (const double temp) {
//   return(
//     6.11 * exp(17.3 * temp / (237.3 + temp)) *
//     4105.3 / pow(237.3 + temp, 2)
//   );
// }

////////////////////////////////////////////////////////////////////////////////
// Slope of the saturation vapor pressure curve, i.e. derivative of saturation
// vapor pressure with resp. to temperature  s = d E(T) / d T
//
// Source: Dyck & Peschke (1995), p. 188, Eq. 11.14
//
// Input:  Temperature (degree celsius)
// Output: Result in (hPa/K)
double slopeSatVapPress (const double temp) {
  return(
    satVapPress_overWater(temp) *
    4098. / pow(237.3 + temp, 2)
  );
}

////////////////////////////////////////////////////////////////////////////////
// Specific humidity (-)
//
// Source: Dyck & Peschke (1995), p. 60, Eq. 4.10
//
// Input:  Air pressure (hPa)
// Input:  Vapor pressure (hPa)
// Output: Result is dimensionless
double specificHumidity (const double pressAir, const double vapPress) {
  return( 0.622 * vapPress / (pressAir - 0.378 * vapPress) );
}

////////////////////////////////////////////////////////////////////////////////
// Latent heat of water evaporation as a function of temperature
//
// Source: Dyck & Peschke (1995), p. 28, Eq. 3.2
//
// Input:  Temperature (degree celsius)
// Output: Result in kJ/kg
double latentHeatEvap (const double temp) {
  return( 2501. - 2.37 * temp );
}

////////////////////////////////////////////////////////////////////////////////
// Psychrometric constant
//
// Source: Dyck & Peschke (1995), p. 188, Eq. 11.15  (with p. 28, Eq. 3.2)
//
// Input:  Temperature (degree celsius)
// Input:  Air pressure (hPa)
// Output: Result in hPa/K
double psychroConst (const double temp, const double airpress) {
  return( 0.016286 * airpress / latentHeatEvap(temp) );
}

////////////////////////////////////////////////////////////////////////////////
// Density of dry air
// Note: The error when using this equation for moist air is low as under normal
//       atmospheric conditions the vapor pressure is small when compared to the
//       air pressure.
//
// Source: Dyck & Peschke (1995), p. 60, Eq. 4.9
// 
// Input:  Temperature (degree celsius)
// Input:  Air pressure (hPa)
// Output: Result in kg/m3
double densityDryAir (const double temp, const double airpress) {
  return( airpress * 0.1 / (0.287 * (273.15 + temp)) );
}


////////////////////////////////////////////////////////////////////////////////
// Mole fraction of water vapor (-).
// CIPM-81/91, eqs. A1.2+A1.3 in Picard et al. (2008), see evap/infos/penman-monteith/References/
// NOTE: In Picard et al. (2008) there are some unit inaccuracies.
//       I am not sure if the formula below is correct in terms of units but
//       results are reasonable.
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double moleFrac_waterVap (
	const double airpress,				// air pressure (hPa)
	const double temp,						// air temperature (°C)
	const double relhum						// relative humidity (%)
) {
		// unit conversions
	double R_h = relhum / 100.;											// in (-)
	double e0_z = satVapPress_overWater(temp)*100.;	// in (Pa)
	double p = airpress *100.;											// in (Pa)
	
	return( R_h * (1.00062 + 3.14e-08 * p + 5.6e-07 * pow(temp,2)) * e0_z / p );
}


////////////////////////////////////////////////////////////////////////////////
// Compressibility factor (-).
// CIPM-81/91, eqs. A1.3+A1.4 in Picard et al. (2008), see evap/infos/penman-monteith/References/
// Range of validity: 600 hPa <= airpress <= 1100 hPa; 15°C <= temp <= 27°C
// NOTE: In Picard et al. (2008) there are some unit inaccuracies.
//       I am not sure if the formula below is correct in terms of units but
//       results are reasonable.
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double f_compress (
	const double airpress,				// air pressure (hPa)
	const double temp,						// air temperature (°C)
	const double relhum						// relative humidity (%)
) {
	// unit conversions
	double p = airpress *100.;											// in (Pa)
	
	// calculate mole fraction of water vapor
	double x_v = moleFrac_waterVap(airpress,temp,relhum);
	
	return(
		1. - p / (temp + T_DEG_K) * 
		(1.58123e-06 - 2.9331e-08 * temp + 1.1043e-10 * pow(temp,2) + (5.707e-06 - 2.051e-08 * temp) * x_v + (1.9898e-04 - 2.376e-06 * temp) * pow(x_v,2)) + 
		pow(p,2) / pow(temp + T_DEG_K, 2) * (1.83e-11 - 0.765e-08 * pow(x_v,2))
	);
}


////////////////////////////////////////////////////////////////////////////////
// Density of moist air (kgm-3).
// CIPM-2007 formula introduced by Picard et al. (2008), see evap/infos/penman-monteith/References/
// NOTE: Tested against David's function densityDryAir() -> values deviate with
//       increasing humidity and temperature (for 20°C, 100% humidity and normal
//			 pressure the value of densityMoistAir() is 8.5% lower than for densityDryAir(),
//			 for 0% humidity values are almost identical as is to be expected).
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double densityMoistAir (
	const double airpress,				// air pressure (hPa)
	const double temp,						// air temperature (°C)
	const double relhum						// relative humidity (%)	
) {
	// unit conversion
	double p = airpress * 100.;
	
	// calculate compressibility factor
	double Z = f_compress(airpress,temp,relhum);
	
	// calculate mole fraction of water vapor
	double x_v = moleFrac_waterVap(airpress,temp,relhum);
	
	return( p * M_DA / (Z * R * (temp + T_DEG_K)) * (1. - x_v * (1. - M_W / M_DA)) );
}


// Precipitation correction after Sevruk (1989) for wind error.
//
// Source: Bremicker (2000), Freiburger Schriften zur Hydrologie No. 11, p. 43, 48
//
// Notes: The function returns a factor for multiplicative correction of
//        precipitation. Only the error due to wind is corrected. Evaporation
//        and interception losses are neglected. The equations derived by
//        Sevruk require the wind speed at the height of the rain gage as input.
//        Here, the wind speed at rain gage height is estimated from speed at
//        measuring height assuming a logarithmic wind profile with zero speed
//        at a selected level z0. The basic functional structure is just
//          u(z) = a * log(z) + b
//        with the two conditions
//          (1) 0 =  a * log(z0) + b
//          (2) observation = a * log(sensor height) + b
//        For two sensors at levels z1 and z2, this leads to
//          u(z1) = u(z2) * log(z1/z0) / log(z2/z0)
//
//        Warning: Using height=0 for the wind sensor will lead to division by
//        zero but this does not make sense anyway. Setting z0=zero is invalid
//        and will be corrected. 
//
// Input: Windspeed (m/s)
// Input: Air temperature (°C)
// Input: Height of rain gage over ground (m), usually 1 m
// Input: Height of wind sensor over ground (m), usually 10 m
// Input: Height over ground where wind speed becomes zero (m), must be > 0,
//        otherwise corrected to a small value.
// Output: Dimensionless correction factor
double factPrecipCorr (const double wind, const double temp,
  const double sensHeight_prec, const double sensHeight_wind, const double height_zeroWind
) {
  double wind_corr= wind * log10(sensHeight_prec/fmax(0.001,height_zeroWind)) /
                           log10(sensHeight_wind/fmax(0.001,height_zeroWind));
  if (temp >= 0.) {
    return ( 1. + (0.015 * wind_corr) );
  } else if (temp >= -8.) {
    return( 1. + (0.150 * pow(wind_corr, 1.18)) );
  } else if (temp >= -27.) {
    return( 1. + (0.280 * pow(wind_corr, 1.30)) );
  } else {
    return( 1. + (0.550 * pow(wind_corr, 1.40)) );
  }
}


////////////////////////////////////////////////////////////////////////////////
// Solar declination (rad) depending on current day of year.
// SWAT manual eq. 1:1.1.2, cites Perrin de Brichambaut (1975).
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double sol_decl (
	const double doy
) {
	return( asin(0.4 * sin(2. * PI / 365. * (doy - 82.))) );
}


////////////////////////////////////////////////////////////////////////////////
// Eccentricity correction factor (-) of the earth's orbit.
// SWAT manual eq. 1:1.1.1, cites Duffie and Beckman (1980).
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double e_corr (
	double doy					// Current day of year
) {
	// omit leap years; this limits the accuracy during leap years
	if (doy > 365.5) {
		doy = 365.;
	}
 
	return( 1. + 0.033 * cos(2. * PI * doy / 365.) );
}


////////////////////////////////////////////////////////////////////////////////
// Daylight time factor, i.e. approximate hour angle of sunset/sunrise (-) needed
// to compute extraterrestrial radiation.
// FAO: http://www.fao.org/docrep/X0490E/x0490e07.htm#calculation%20procedures eq. 25
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double dayTime_fac (
	const double doy,					// Current day of year
	const double lat					// Latitude (decimal degree)
) {
	// latitude in (rad)
	double phi = abs(lat * PI / 180.);
	
	double sunrise = -1. * tan(sol_decl(doy)) * tan(phi);
	
	if (sunrise > 1.)
		sunrise = 1.;
	
	if (sunrise < -1.)
		sunrise = -1;
  
  return(acos(sunrise));
}


////////////////////////////////////////////////////////////////////////////////
// Extraterrestrial (i.e. top of atmosphere) radiation (Wm-2). To calculate daily
// mean values only.
// SWAT manual sec. 1:1.2.1
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double rad_extraterr_daily (
	const double doy,					// Current day of year
	const double lat					// Latitude (decimal degree)
) {

	// solar declination (rad)
	double delta = sol_decl(doy);
	
	// eccentricity correction factor (-)
	double E_0 = e_corr(doy);
	
	// latitude in (rad)
	double phi = lat * PI / 180.;
	
	// daylight time factor (-)
	double t_sr = dayTime_fac(doy,lat);
	
	return( 1./PI * SOLAR_C * E_0 * (t_sr * sin(delta) * sin(phi) + cos(delta) * cos(phi) * sin(t_sr)) );
}


////////////////////////////////////////////////////////////////////////////////
// Extraterrestrial (i.e. top of atmosphere) radiation (Wm-2). Hourly values.
//
// Calculation following http://www.fao.org/docrep/X0490E/x0490e07.htm#calculation%20procedures
// eqs. 28-33.
//
// NOTE: To account for daylight saving time 'utc_add' has to be given. It is
// suggested to supply it as time series forcing (as 'doy') to account for
// changing values over the year due to daylight saving time.
// 
// <tobias.pilz@uni-potsdam.de>, JAN 2016
////////////////////////////////////////////////////////////////////////////////

double rad_extraterr_hourly (
	const double doy,					// Current day of year
	const double lat,					// Latitude (decimal degree)
	const int hour,						// hour of day in local time (including daylight daving time), range of [0..23]
	const int utc_add,				// Deviation of local time zone from UTC; may vary over the year due to daylight saving time; range of [-12..14] (hours)
	const double L_m					// Longitude of the location of interest (decimal degrees west of Greenwich, e.g. Greenwich: 0°, Berlin: 345°, New York: 75°)
) {
	// check input
	if ( (hour < 0) || (hour > 23) ) {
		stringstream errmsg;
		errmsg << "Computation of extraterrestrial radiation (hourly): Required hour of day not in the accepted range of [0..23]!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	if ( (utc_add < 0) || (utc_add > 23) ) {
		stringstream errmsg;
		errmsg << "Computation of extraterrestrial radiation (hourly): Required deviation from UTC not in the accepted range of [-12..14]!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	if ( (L_m < 0.) || (L_m > 360.) ) {
		stringstream errmsg;
		errmsg << "Computation of extraterrestrial radiation (hourly): Longitude in degrees west of Greenwich required, i.e. [0..360], e.g. Greenwich: 0°, New York: 75°, Berlin: 345°!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	// compute centre of current time zone (decimal degrees west of Greenwich)
	int L_z = (360 - utc_add * 15) % 360;
	
	// seasonal correction factor, eqs. 32, 33
	double b = 2. * PI * (doy - 81.) / 364.;
	double S_c = 0.1645 * sin(2.*b) - 0.1255 * cos(b) - 0.025 * sin(b);
	
	// solar time angle at the midpoint of the hour of day (rad), eq. 31
	double hour_mid = hour + 0.5;
	double w_mid = PI / 12. * ( (hour_mid + 0.06667 * (L_z-L_m) + S_c) - 12. );
	
	// solar time angles at begin and end of hour, eqs. 29, 30
	double w_ini = w_mid - PI / 24.;
	double w_end = w_mid + PI / 24.;
	
	// sunset/sunrise hour angle (rad)
	double w_s = dayTime_fac(doy,lat);
	
	// if no sun is shining (+/- half an hour accuracy) return zero
	if ( (w_mid < -1.*w_s) || (w_mid > w_s) )
		return(0.);
	
	// solar declination (rad)
	double delta = sol_decl(doy);
	
	// eccentricity correction factor (-)
	double E_0 = e_corr(doy);
	
	// latitude in (rad)
	double phi = lat * PI / 180.;
	
	return( 1./PI * SOLAR_C * E_0 * ( (w_end-w_ini) * sin(delta) * sin(phi) + cos(delta) * cos(phi) * (sin(w_end)-sin(w_ini)) ) );
}


////////////////////////////////////////////////////////////////////////////////
// Global radiation from sunshine duration or cloudiness, extraterrestrial 
// radiation and emprical parameters (Angström  coefficients) (Wm-2).
// Dyck & Peschke (1995), p. 30, eq. 3.10
// TODO: Account for observed cloudiness as input (nn NOT equal to 1 - cloud).
// ATTENTION: For daily resolution only! Resturns daily average value.
// Parameters: 
//			- Dyck & Peschke (1995): radex_a = 0.25, radex_b = 0.5
//			- WASA (Guentner, 2002): radex_a = 0.18, radex_b = 0.55
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double calc_glorad (
	const double radex,					// Incoming extraterrestrial radiation (i.e. at top of atmosphere) (W/m2)
	const double sundur,				// Sunshine duration of current day (h)
	const double cloud,					// Cloudiness (%)
	const double lat,						// Latitude (decimal degree)
	const double doy,						// Current day of year
	const double radex_a,				// Angstrom coefficient (share of radex on glorad under clouds) (-)
	const double radex_b				// Angstrom coefficient (radex_a+radex_b = share of radex on glorad under clear sky) (-)
) {
	double sundur_max = -9999.;
	double nn = -9999.;
	
	if( (radex_a + radex_b) > 1.0 ) {
		stringstream errmsg;
		errmsg << "Calculation of actual incoming shortwave radiation: Parameters (radex_a + radex_b) > 1.0 which is physically not possible!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	// use cloudiness or sunshine duration
	if (sundur < -90.) {
		stringstream errmsg;
		errmsg << "Calculation of global radiation: 'sundur' has to be given as input (cloudiness not yet fully supported)!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	} else {
		// calculate maximum possible sunshine duration for current day of year
		sundur_max = dayTime_fac(doy,lat) * 24./PI;
		nn = sundur / sundur_max;
	}
	
	return( radex * (radex_a + radex_b * nn) );
}


////////////////////////////////////////////////////////////////////////////////
// Maximum global radiation in case of clear sky (Wm-2).
// Dyck & Peschke (1995), p. 30, eq. 3.10 (method 1)
// OR Allen (2005), ASCE standard etp, eq. 19
// Parameters: 
//			- Dyck & Peschke (1995): radex_a = 0.25, radex_b = 0.5
//			- WASA (Guentner, 2002): radex_a = 0.18, radex_b = 0.55
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double calc_glorad_max (
	const int choice,						// Flag which method to use (see above)
	const double radex,					// Incoming extraterrestrial radiation (i.e. at top of atmosphere) (W/m2)
	const double radex_a,				// Angstrom coefficient (share of radex on glorad under clouds) (-)
	const double radex_b,				// Angstrom coefficient (radex_a+radex_b = share of radex on glorad under clear sky) (-)
	const double elev						// Elevation above sea level (m)
) {
	double res = -9999.;
	
	if( (radex_a + radex_b) > 1.0 ) {
		stringstream errmsg;
		errmsg << "Computation of maximum incoming short-wave radiation: Parameters (radex_a + radex_b) > 1.0 which is physically not possible!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
			
			
	switch(choice) {
		case 1: // Angstroem
			res = radex * (radex_a + radex_b);
			break;
			
		case 2: // Allen (2005), ASCE standard etp, eq. 19
			res = radex * (0.75 + 2e-5 * elev);
			break;
			
		default :
			stringstream errmsg;
      errmsg << "Invalid choice to calculate maximum incoming short-wave radiation! Currently supported is one of {1,2}.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e); 
	}
	
	return(res);
}


////////////////////////////////////////////////////////////////////////////////
// Net emissivity between land surface and atmosphere (-).
// Dyck & Peschke (1995), p. 31, eq. 3.12
// Maidment (1993) eqs. 4.2.8, 4.2.9
// Parameters:
//		- SWAT manual (2011): a = 0.34, b = -0.139
//		- Dyck & Peschke (1995): 0.34 < a < 0.44, -0.25 < b < -0.14; frequently used: a = 0.34, b = -0.14
//		- WASA (Guentner, 2002): a = 0.52, b = -0.065 (value of model code, in manual it is + 0.065)
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double net_emiss (
	const double temp,			// Air temperature (°C)
	const double relhum,		// Relative humidity (%)
	const double a,					// coefficient
	const double b					// coefficient
) {
	// two variants depending on availability of relhum measurements
	if(relhum < -90.)
		return( -0.02 + 0.261 * exp(-7.77e-04 * pow(temp,2)) );
	else {
		// calculate vapor pressure from temperature and relative humidity (kPa)
		double e = vapPress_overWater(temp, relhum) / 10.;
		
		return( a + b * sqrt(e) );
	}
}


////////////////////////////////////////////////////////////////////////////////
// Cloudiness correction factor (-).
// Dyck & Peschke (1995), p. 31, combining eqs. 3.10 and 3.13
// ATTENTION: Works for daily (or sub-daily values over daytime) only, i.e. when glorad_max > 0.
// To avoid this problem assume a constant cloudiness over night, i.e. use the last
// values of radiation before sunset. In a model engine this can be realized by using
// glorad_max and glorad as a state variables, i.e. update the state variables over daytime on
// time steps where glorad_max > 0 and don't update it otherwise (i.e. over nighttime)
// to ensure the last value before sunset is used over the all night time steps.
// Coefficients:
//		- SWAT manual (2011): a = 0.9, b = 0.1
//		- Dyck & Peschke (1995): arid areas: a = 1.35, b = -0.35; humid areas: a = 1, b = 0
//		- WASA (Guentner, 2002): a = 1.35, b = -0.35 (indirectly by using eq. 3.15 from Dyck & Peschke)
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double f_cloud (
	const double glorad,		// Downward short-wave radiation (W/m2)
	const double glorad_max, // Downward short-wave radiation under clear (cloudless) sky (Wm-2)
	const double a,					// Parameter (radiation coefficient a for cloudless sky)
	const double b					// Parameter (radiation coefficient b for cloudless sky)
) {
	
	if(glorad_max < 1e-6) {
		stringstream errmsg;
		errmsg << "Computation of cloudiness correction factor: Input 'glorad_max' is zero which is not allowed!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	if( (a + b) > 1.0 ) {
		stringstream errmsg;
		errmsg << "Computation of cloudiness correction factor: Parameters (a + b) > 1.0 which is physically not possible!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	return( a * glorad / glorad_max + b);
}


////////////////////////////////////////////////////////////////////////////////
// Incoming net long-wave radiation from Stefan-Boltzmann law (Wm-2).
// Negative values indicate net loss, positive values net gain of energy at surface.
// SWAT manual sec. 1:1.2.5.2, Dyck & Peschke (1995), pp. 30-31
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double net_longrad (
	const double temp,			// Air temperature (°C)
	const double relhum,		// Relative humidity (%)
	const double glorad,		// Downward short-wave radiation (W/m2)
	const double glorad_max, // Downward short-wave radiation under clear (cloudless) sky (Wm-2)
	const double emis_a,		// Coefficient a for calculating net emissivity (-)
	const double emis_b,		// Coefficient b for calculating net emissivity (-)
	const double fcorr_a,		// Coefficient a for calculating cloudiness correction factor (-)
	const double fcorr_b		// Coefficient b for calculating cloudiness correction factor (-)
) {
	// calculate cloud correction factor
	double f = f_cloud(glorad, glorad_max, fcorr_a, fcorr_b);

	
	// calculate net emissivity (-)
	double em = net_emiss(temp, relhum, emis_a, emis_b);
	
	return( -1. * f * em * SIGMA * pow(temp + T_DEG_K,4) );
}


////////////////////////////////////////////////////////////////////////////////
// Incoming net radiation from incoming shortwave radiation, albedo and net
// incoming longwave radiation (Wm-2).
// SWAT manual eq. 1:1.2.12, Dyck & Peschke (1995), p. 30, eq. 3.8
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

double net_rad (
	const double glorad,  	// Downward short-wave radiation (Wm-2)
	const double glorad_max, // Downward short-wave radiation under clear (cloudless) sky (Wm-2)
	const double temp,			// Air temperature (°C)
	const double relhum,		// Relative humidity (%)
	const double alb,				// albedo (-)
	const double emis_a,		// Coefficient a for calculating net emissivity (-)
	const double emis_b,		// Coefficient b for calculating net emissivity (-)
	const double fcorr_a,		// Coefficient a for calculating cloudiness correction factor (-)
	const double fcorr_b		// Coefficient b for calculating cloudiness correction factor (-)
) {
	
	// estimate net incoming long-wave radiation (Wm-2)
	double netlongrad = net_longrad(temp,relhum,glorad,glorad_max,emis_a,emis_b,fcorr_a,fcorr_b);
	
	double res = (1. - alb) * glorad + netlongrad;
	
	return (res);
}


////////////////////////////////////////////////////////////////////////////////
// Soil heat flux (Wm-2).
// Simplified assumption according to FAO: during daytime 10% (WASA: 20%) of net radiation,
// during nighttime 50% (WASA: 70%) of netradiation. For sub-daily time-steps only!
// Over a whole day it is assumed to be zero!
// http://www.fao.org/docrep/X0490E/x0490e07.htm#calculation%20procedures eq. 45, 46.
// <tobias.pilz@uni-potsdam.de>, FEB 2015
////////////////////////////////////////////////////////////////////////////////

double soil_heatflux (
	const double net_rad,		// Downward net (long- and short-wave) radiation (Wm-2)
	const double f_day,			// Fraction of net_rad over daytime (in case of sub-daily application) (-)
	const double f_night,		// Fraction of net_rad over nighttime (in case of sub-daily application) (-)
	const int daynight,			// flag: 1: calculation for day time; 0: calculation for night time
	const unsigned int delta_t				// time step length (s)
) {
	double res = -9999.;
	
	// calculate only for sub-daily timesteps
	if (delta_t < 86400) {
		// if global radiation is zero -> nighttime
		if (daynight == 0)
			res = f_night * net_rad;
		
		// if global radiation > zero -> daytime
		if (daynight == 1)
			res = f_day * net_rad;
		
	} else { // return zero for daily time steps
		res = 0.;
	}
	
	return(res);
}


#endif

