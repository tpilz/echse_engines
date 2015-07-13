
#ifndef SNOW_H
#define SNOW_H

#include <cmath>

#include "../../meteo/meteo.h"

////////////////////////////////////////////////////////////////////////////////
// Definition of stoichiometry factors
////////////////////////////////////////////////////////////////////////////////

// Conversion of precipitation mass flux (m/s) to energy flux (kJ/m2/s)
// Unit of result: kJ/m3
double f_prec (
  const double tempAir,       // Air temperature (°C)
  const double tempAir_crit   // Threshold temp. for rain-/snowfall (°C)
) {
  if (tempAir > tempAir_crit) {
    // 4180. = 1000. * 4.18 = densityWater [kg/m3] * specHeatWater [kJ/kg/K]
    // 333.5e+03. = 1000. * 333.5 = densityWater [kg/m3] * meltHeatIce [kJ/kg]
    return( 4180. * fmax(tempAir, 0.) + 333.5e+03 );
  } else {
    // 2090. = 1000. * 2.09 = densityWater [kg/m3] * specHeatIce [kJ/kg/K]
    return( 2090. * fmin(tempAir, 0.) );
  }
}

// Conversion of sublimation mass flux (m/s) to energy flux (kJ/m2/s)
// Unit of result: kJ/m3
// 2837.e+03 = 1000. * 2837. = densityWater [kg/m3] * sublHeatIce [kJ/kg]
double f_subl() { return( 2837.e+03 ); }

// Conversion of meltwater loss mass flux (m/s) to energy flux (kJ/m2/s)
// Unit of result: kJ/m3
// 333.5e+03. = 1000. * 333.5 = densityWater [kg/m3] * meltHeatIce [kJ/kg]
double f_flow() { return ( 333.5e+03 ); }

////////////////////////////////////////////////////////////////////////////////
// General relations for estimation of derived variables
////////////////////////////////////////////////////////////////////////////////

// Fraction of liquid water (mass water / (mass water + mass ice))
// Unit of result: Dimensionless, range 0...1
double snowLiquidFrac (
  const double snowEnergyCont,  // Snow energy content (kJ/m2)
  const double snowWaterEquiv   // Snow water equivalent (m)
) {
  // 333500. = 333.5 * 1000. = meltHeatOfIce (kJ/kg) * densWater (kg/m3)
  // The result is bounded to the range 0...1.
  // If there is no snow, a fraction of 1 is returned.
  if (snowWaterEquiv > 0.) {
    return( fmin(1., fmax(0., snowEnergyCont / (snowWaterEquiv * 333500.))) );
  } else {
    return( 1. );
  }
}

// Mean temperature of the snow pack
// Unit of result: °C  (Range: -Inf ... 0)
double snowTemp_mean (
  const double snowEnergyCont,  // Snow energy content (kJ/m2)
  const double snowWaterEquiv,  // Snow water equivalent (m)
  const double soilDepth,       // Depth of interacting soil layer (m)
  const double soilDens,        // Density of soil (kg/m3)
  const double soilSpecHeat     // Spec. heat capacity of soil (kJ/kg/K)
) {
  if (snowWaterEquiv > 0.) {
    // If the snow pack is free of liquid water
    if (snowEnergyCont < 0.) {
      // 2090. = 1000. * 2.09 = WaterDensity (kg/m3) * specHeatCapIce (kJ/kg/K)
      return ( snowEnergyCont / (snowWaterEquiv * 2090. +
        soilDepth * soilDens * soilSpecHeat) );
    // If the snow pack contains some liquid water
    } else {
      return ( 0. );
    }
    // Note: Temperature for the case where all water is liquid is not computed.
    // Note: If there is no snow, a temperature of zero is returned.
  } else {
    return ( 0. );
  }
}

// Snow surface temperature
// Unit of result: °C  (Range: -Inf ... 0)
double snowTemp_surf (
  const double tempSnow_mean,  // Mean temperature of the snow pack (°C)
  const double tempAir,        // Air temperature (°C)
  const double weightAirTemp   // Weighting param. for air temp. (-) in 0...1
) {
  if (tempSnow_mean < 0.) {
    return( fmin(0., (1.-weightAirTemp) * tempSnow_mean +
      weightAirTemp * tempAir) );
  } else {
    return( 0. );
  }
}

////////////////////////////////////////////////////////////////////////////////
// Definition of rate expressions
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// Part 1: Radiation flux rates
//------------------------------------------------------------------------------

// Short-wave radiation balance
// Unit of result: W/m2
double R_netS (
  const double shortRad,      // Incoming short wave radiation (W/m2)
  const double albedo,        // Albedo (-)
  const double snowWaterEquiv // Snow water equivalent (m)
) {
  if (snowWaterEquiv > 0.) {
    return( shortRad * (1. - albedo) );
  } else {
    return( 0. );
  }
}

// Long-wave radiation balance
// Unit of result: W/m2
double R_netL (
  const double tempSnow_surf,     // Temperature of the snow surface (°C)
  const double emissivitySnowMin, // Minimum snow emissivity used for old snow (-)
  const double emissivitySnowMax, // Maximum snow emissivity used for new snow (-)
  const double tempAir,           // Air temperature (°C)
  const double relHumid,          // Relative humidity (%)
  const double cloudCoverage,     // Cloud cover (0 = clear sky, 1 = fully covered)
  const double albedo,            // Albedo
  const double albedoMin,         // Minimum albedo used for old snow (-)
  const double albedoMax,         // Maximum albedo used for new snow (-)
  const double snowWaterEquiv     // Snow water equivalent (m)
) {
  if (snowWaterEquiv > 0.) {
    // Outgoing part
    // Note: The snow emissivity decreases with the age of the snow cover (Dyck &
    //       Peschke, 1995). We use the dynamically computed albedo as an indicator
    //       for the age of the snow pack. The constant 0.001 is to prevent a Div0
    //       in case of a constant albedo (i.e. albedoMax = albedoMin).
    // 5.67e-08 = Stefan-Boltzmann constant (W/m2/K4)
    double R_outL = (emissivitySnowMin + (emissivitySnowMax - emissivitySnowMin) *
      (albedo-albedoMin+0.001) / (albedoMax-albedoMin+0.001)) *
      5.67e-08 * pow((tempSnow_surf + 273.15), 4.);
    // Incoming part 1 - Clear sky emission
    // This is the Stefan-Bolzmann equation with clear sky emissivity estimated
    // from vaporPressure using the empirical Brunt-formula.
    double R_inL_clear = (0.51 + 0.066 * sqrt(vapPress_overIce(tempAir,relHumid))) * 
      5.67e-08 * pow((tempAir + 273.15), 4.);
    // Incoming part 2 - Cloud emission
    // This is the stefan-Boltzmann equation with the emissivity set to 1 (clouds
    // treated as black body). The cloud temperature is approximated by the dew-
    // point temperature.
    double R_inL_cloud = 1.0 * 5.67e-08 *
      pow((dewTemp_overIce(vapPress_overIce(tempAir,relHumid)) + 273.15), 4.);
    // Radiation balance (net radiation)
    return( (1.-cloudCoverage)*R_inL_clear + cloudCoverage*R_inL_cloud - R_outL);    
  } else {
    return ( 0. );
  }
}

// Soil heat flux
// Unit of result: W/m2
double R_soil() { return ( 0. ); }

// Sensible heat flux
// Unit of result: W/m2
double R_sens (
  const double tempSnow_surf,  // Temperature of the snow surface (°C)
  const double tempAir,        // Air temperature (°C)
  const double pressAir,       // Air pressure (hPa)
  const double windSpeed,      // Wind speed (m/s)
  const double a0,             // Empirical coeff. (m/s)
  const double a1,              // Empirical coeff. (-)
  const double snowWaterEquiv     // Snow water equivalent (m)
) {
  if (snowWaterEquiv > 0.) {
    // (a0 + a1 * windSpeed) = Empirical estimation of transfer coefficient (m/s)
    // 1005. = Specific heat capacity of air in J/kg/K
    return ( (a0 + a1 * windSpeed) * densityDryAir(tempAir,pressAir) *
      1005. * (tempAir - tempSnow_surf));
  } else {
    return ( 0. );
  }
}

//------------------------------------------------------------------------------
// Part 2: Mass flux rates
//------------------------------------------------------------------------------

// Precipitation mass flux
// Unit or result: m/s
double M_prec (
  const double precipSumMM,   // Precipitation sum (mm / referenceInterval)
  const double precipSeconds  // Length of referenceInterval (seconds)
) {
  return ( precipSumMM / 1000. / precipSeconds );
}

// Sublimation mass flux
// Unit or result: m/s
double M_subl (
  const double tempSnow_surf,  // Temperature of the snow surface (°C)
  const double tempAir,        // Air temperature (°C)
  const double pressAir,       // Air pressure (hPa)
  const double relHumid,       // Relative humidity (%)
  const double windSpeed,      // Wind speed (m/s)
  const double a0,             // Empirical coeff. (m/s)
  const double a1,             // Empirical coeff. (-)
  const double snowWaterEquiv  // Snow water equivalent (m)
) {
  if (snowWaterEquiv > 0.) {
    // (a0 + a1 * windSpeed) = Empirical estimation of transfer coefficient (m/s)
    // 1000. = Density of water (kg/m3)
    // 100. = Relative humidity at saturation (%)
    return (
      fmax(0., // Negative flux rates are set to zero (no freezing air moisture)
      (a0 + a1 * windSpeed) * densityDryAir(tempAir,pressAir) / 1000. *
      (specificHumidity(pressAir, vapPress_overIce(tempSnow_surf,100.)) -
       specificHumidity(pressAir, vapPress_overIce(tempAir,relHumid)))
      )
    );
  } else {
    return ( 0. );
  }
}

// Meltwater flux
// Unit or result: m/s
double M_flow (
  const double snowLiquidFrac,// Fraction of liquid water (-)
  const double kSatSnow,      // Saturated hydraulic conductivity of snow (m/s)
  const double densDrySnow,   // Density of dry snow (kg/m3)
  const double specCapRet,    // Capill. retent. vol. as fract. of solid SWE (-)
  const double snowWaterEquiv // Snow water equivalent (m)
) {
  if (snowWaterEquiv > 0.) {
    // Relative saturation (-)
    // Don't allow 100% liquid water as this will cause a division by zero
    double rss= (fmin(0.99, snowLiquidFrac)/(1.-fmin(0.99, snowLiquidFrac)) - specCapRet) /
      (1000./densDrySnow - 1000./922. - specCapRet);
    // 1000. = Density of water (kg/m3)
    // 922. = Density of ice (kg/m3)
    // Fix negative values of the relative saturation (-)
    return ( kSatSnow * pow(fmax(0., rss), 3.) );
  } else {
    return ( 0. );
  }
}

//------------------------------------------------------------------------------
// Part 3: Other rates
//------------------------------------------------------------------------------

// Change rate of albedo
// Unit of result: 1/s
double G_alb (
  const double albedo,           // Current albedo (-)
  const double precipSumMM,      // Precipitation sum (mm / referenceInterval)
  const double precipSeconds,    // Length of referenceInterval (seconds)
  const double tempAir,          // Air temperature (°C)
  const double tempAir_crit,     // Threshold temperature for rain/snow (°C)
  const double albedoMin,        // Albedo of oldest snow (-)
  const double albedoMax,        // Albedo of fresh snow (-)
  const double agingRate_tAirPos, // Aging rate for air temperatures > 0 (1/s)
  const double agingRate_tAirNeg, // Aging rate for air temperatures < 0 (1/s)
  const double snowWaterEquiv     // Snow water equivalent (m)
) {
  if (snowWaterEquiv > 0.) {
    // Surface renewal if snow falls
    // Time of renewal set to the reference interval to keep the change rate reasonably
    // small in comparison to other change rates (prevent stiff ODE system). Complete
    // renewal is assumed for a snowfall of 10 mm (as in Utah Energy Balance Model) and
    // a partly renewal is assumed for less precipitation.
    if ((precipSumMM > 0.) & (tempAir < tempAir_crit)) {
      return ( (albedoMax - albedo) / precipSeconds * fmin(1., precipSumMM/10.) );  // This is positive
    // Surface aging
    } else {
      if (tempAir >= 0.) {
        return ( agingRate_tAirPos * (albedoMin - albedo) );  // This is negative
      } else {
        return ( agingRate_tAirNeg * (albedoMin - albedo) );  // This is negative
      }
    }
  } else {
    return ( 0. );
  }
}

////////////////////////////////////////////////////////////////////////////////
// Derivatives of the snow model's state variables with respect to time
////////////////////////////////////////////////////////////////////////////////

void snowModel_derivs(
  // Inputs
  // Part 1: Forcings
  const double precipSumMM,
  const double shortRad,
  const double tempAir,
  const double pressAir,
  const double relHumid,
  const double windSpeed,
  const double cloudCoverage,
  // Part 2: Parameters
  const double precipSeconds,
  const double a0,
  const double a1,
  const double kSatSnow,
  const double densDrySnow,
  const double specCapRet,
  const double emissivitySnowMin,
  const double emissivitySnowMax,
  const double tempAir_crit,
  const double albedoMin,
  const double albedoMax, 
  const double agingRate_tAirPos,
  const double agingRate_tAirNeg,
  const double soilDepth,
  const double soilDens,
  const double soilSpecHeat,
  const double weightAirTemp,
  // Part 3: States
  const double snowEnergyCont,
  const double snowWaterEquiv,
  const double albedo,
  // Outputs
  double &ddt_sec, // Result unit: kJ/m2/s
  double &ddt_swe, // Result unit: m/s
  double &ddt_alb, // Result unit: 1/s
  double &flux_melt, // Result unit: m/s, Useful for water balance check
  double &flux_subl  // Result unit: m/s, Useful for water balance check
) {
  // Derived variables
  double TEMP_MEAN= snowTemp_mean(snowEnergyCont, snowWaterEquiv,
    soilDepth, soilDens, soilSpecHeat);
  double TEMP_SURF= snowTemp_surf(TEMP_MEAN, tempAir, weightAirTemp);
  double LIQU_FRAC= snowLiquidFrac(snowEnergyCont, snowWaterEquiv);
  // Rate expressions used multiple times
  double M_P= M_prec(precipSumMM, precipSeconds);
  double M_S= M_subl(TEMP_SURF, tempAir, pressAir, relHumid, windSpeed, a0, a1, snowWaterEquiv);
  double M_F= M_flow(LIQU_FRAC, kSatSnow, densDrySnow, specCapRet, snowWaterEquiv);
  // Computation of derivatives
  ddt_sec= 0.001 * (
    R_netS(shortRad, albedo, snowWaterEquiv) +
    R_netL(TEMP_SURF, emissivitySnowMin, emissivitySnowMax, tempAir, relHumid,
      cloudCoverage, albedo, albedoMin, albedoMax, snowWaterEquiv) +
    R_soil() +
    R_sens(TEMP_SURF, tempAir, pressAir, windSpeed, a0, a1, snowWaterEquiv)
    )
    + f_prec(tempAir, tempAir_crit) * M_P
    - f_subl() * M_S
    - f_flow() * M_F;
  ddt_swe= M_P - M_S - M_F;
  ddt_alb= G_alb(albedo, precipSumMM, precipSeconds, tempAir, tempAir_crit,
      albedoMin, albedoMax, agingRate_tAirPos, agingRate_tAirNeg, snowWaterEquiv);
  flux_melt= M_F;
  flux_subl= M_S;
}

////////////////////////////////////////////////////////////////////////////////
// Function to return the fundamental variables of the snow model
// (for debugging purposes and in-depth analyses of the behaviour)
////////////////////////////////////////////////////////////////////////////////

void snowModel_debug(
  // Inputs
  // Part 1: Forcings
  const double precipSumMM,
  const double shortRad,
  const double tempAir,
  const double pressAir,
  const double relHumid,
  const double windSpeed,
  const double cloudCoverage,
  // Part 2: Parameters
  const double precipSeconds,
  const double a0,
  const double a1,
  const double kSatSnow,
  const double densDrySnow,
  const double specCapRet,
  const double emissivitySnowMin,
  const double emissivitySnowMax,
  const double tempAir_crit,
  const double albedoMin,
  const double albedoMax, 
  const double agingRate_tAirPos,
  const double agingRate_tAirNeg,
  const double soilDepth,
  const double soilDens,
  const double soilSpecHeat,
  const double weightAirTemp,
  // Part 3: States
  const double snowEnergyCont,
  const double snowWaterEquiv,
  const double albedo,
  // Outputs
  double &TEMP_MEAN,
  double &TEMP_SURF,
  double &LIQU_FRAC,
  double &flux_R_netS,
  double &flux_R_netL,
  double &flux_R_soil,
  double &flux_R_sens,
  double &stoi_f_prec,
  double &stoi_f_subl,
  double &stoi_f_flow,
  double &flux_M_prec,
  double &flux_M_subl,
  double &flux_M_flow,
  double &rate_G_alb
) {
  // Derived variables
  TEMP_MEAN= snowTemp_mean(snowEnergyCont, snowWaterEquiv,
    soilDepth, soilDens, soilSpecHeat);
  TEMP_SURF= snowTemp_surf(TEMP_MEAN, tempAir, weightAirTemp);
  LIQU_FRAC= snowLiquidFrac(snowEnergyCont, snowWaterEquiv);
  // Mass fluxes
  flux_M_prec= M_prec(precipSumMM, precipSeconds);
  flux_M_subl= M_subl(TEMP_SURF, tempAir, pressAir, relHumid, windSpeed, a0, a1, snowWaterEquiv);
  flux_M_flow= M_flow(LIQU_FRAC, kSatSnow, densDrySnow, specCapRet, snowWaterEquiv);
  // Radiation fluxes
  flux_R_netS= R_netS(shortRad, albedo, snowWaterEquiv);
  flux_R_netL= R_netL(TEMP_SURF, emissivitySnowMin, emissivitySnowMax, tempAir,
    relHumid, cloudCoverage, albedo, albedoMin, albedoMax, snowWaterEquiv);
  flux_R_soil= R_soil();
  flux_R_sens= R_sens(TEMP_SURF, tempAir, pressAir, windSpeed, a0, a1, snowWaterEquiv);
  // Stoichiometry factors
  stoi_f_prec= f_prec(tempAir, tempAir_crit);
  stoi_f_subl= f_subl();
  stoi_f_flow= f_flow();
  // Other rates
  rate_G_alb= G_alb(albedo, precipSumMM, precipSeconds, tempAir, tempAir_crit,
      albedoMin, albedoMax, agingRate_tAirPos, agingRate_tAirNeg, snowWaterEquiv);
}

#endif

