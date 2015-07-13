
#include "../snow.h"

extern "C" {

  // The cpp-functions to be called from R go here. These are just wrappers for
  // the functions defined in "snow.h". The basic difference is that in the
  // wrapper versions all arguments are pointers.

  void snowModel_derivs(
    // Inputs
    // Part 1: Forcings
    const double* precipSumMM,
    const double* shortRad,
    const double* tempAir,
    const double* pressAir,
    const double* relHumid,
    const double* windSpeed,
    const double* cloudCoverage,
    // Part 2: Parameters
    const double* precipSeconds,
    const double* a0,
    const double* a1,
    const double* kSatSnow,
    const double* densDrySnow,
    const double* specCapRet,
    const double* emissivitySnowMin,
    const double* emissivitySnowMax,
    const double* tempAir_crit,
    const double* albedoMin,
    const double* albedoMax, 
    const double* agingRate_tAirPos,
    const double* agingRate_tAirNeg,
    const double* soilDepth,
    const double* soilDens,
    const double* soilSpecHeat,
    const double* weightAirTemp,
    // Part 3: States
    const double* snowEnergyCont,
    const double* snowWaterEquiv,
    const double* albedo,
    // Outputs
    double* ddt_sec, // Result unit: kJ/m2/s
    double* ddt_swe, // Result unit: m/s
    double* ddt_alb, // Result unit: 1/s
    double* flux_melt, // Result unit: m/s, Useful for water balance check
    double* flux_subl  // Result unit: m/s, Useful for water balance check
  ) {
    snowModel_derivs(
      // In
      *precipSumMM, *shortRad, *tempAir, *pressAir, *relHumid, *windSpeed,
      *cloudCoverage, *precipSeconds, *a0, *a1, *kSatSnow, *densDrySnow,
      *specCapRet, *emissivitySnowMin, *emissivitySnowMax, *tempAir_crit,
      *albedoMin, *albedoMax, *agingRate_tAirPos, *agingRate_tAirNeg,
      *soilDepth, *soilDens, *soilSpecHeat, *weightAirTemp, *snowEnergyCont,
      *snowWaterEquiv, *albedo,
      // Out
      *ddt_sec, *ddt_swe, *ddt_alb, *flux_melt, *flux_subl
    );
  }

  void snowModel_debug(
    // Inputs
    // Part 1: Forcings
    const double* precipSumMM,
    const double* shortRad,
    const double* tempAir,
    const double* pressAir,
    const double* relHumid,
    const double* windSpeed,
    const double* cloudCoverage,
    // Part 2: Parameters
    const double* precipSeconds,
    const double* a0,
    const double* a1,
    const double* kSatSnow,
    const double* densDrySnow,
    const double* specCapRet,
    const double* emissivitySnowMin,
    const double* emissivitySnowMax,
    const double* tempAir_crit,
    const double* albedoMin,
    const double* albedoMax, 
    const double* agingRate_tAirPos,
    const double* agingRate_tAirNeg,
    const double* soilDepth,
    const double* soilDens,
    const double* soilSpecHeat,
    const double* weightAirTemp,
    // Part 3: States
    const double* snowEnergyCont,
    const double* snowWaterEquiv,
    const double* albedo,
    // Outputs
    double* TEMP_MEAN,
    double* TEMP_SURF,
    double* LIQU_FRAC,
    double* flux_R_netS,
    double* flux_R_netL,
    double* flux_R_soil,
    double* flux_R_sens,
    double* stoi_f_prec,
    double* stoi_f_subl,
    double* stoi_f_flow,
    double* flux_M_prec,
    double* flux_M_subl,
    double* flux_M_flow,
    double* rate_G_alb
  ) {
    snowModel_debug(
      // In
      *precipSumMM, *shortRad, *tempAir, *pressAir, *relHumid, *windSpeed,
      *cloudCoverage, *precipSeconds, *a0, *a1, *kSatSnow, *densDrySnow,
      *specCapRet, *emissivitySnowMin, *emissivitySnowMax, *tempAir_crit,
      *albedoMin, *albedoMax, *agingRate_tAirPos, *agingRate_tAirNeg,
      *soilDepth, *soilDens, *soilSpecHeat, *weightAirTemp, *snowEnergyCont,
      *snowWaterEquiv, *albedo,
      // Out
      *TEMP_MEAN, *TEMP_SURF, *LIQU_FRAC, *flux_R_netS, *flux_R_netL,
      *flux_R_soil, *flux_R_sens, *stoi_f_prec, *stoi_f_subl, *stoi_f_flow,
      *flux_M_prec, *flux_M_subl, *flux_M_flow, *rate_G_alb
    );
  }

} // extern "C"

