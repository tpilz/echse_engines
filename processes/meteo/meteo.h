
#ifndef METEO_H
#define METEO_H

#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// Saturation vapor pressure as a function of temperature
//
// Source: Dyck & Peschke (1995), p. 58, Eq. 4.6
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
// Input:  Temperature (degree celsius)
// Output: Result in (hPa/K)
double slopeSatVapPress (const double temp) {
  return(
    6.11 * exp(17.3 * temp / (237.3 + temp)) *
    4105.3 / pow(237.3 + temp, 2)
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

#endif

