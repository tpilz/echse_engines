
#ifndef EVAP_H
#define EVAP_H

// Potential reference crop ET after Makkink corrected by crop factor.
// Returns: Rate of potential evapotranspiration (m/s)

inline double et_pot_makkink (
  const double glorad,       // Downward short-wave radiation (W/m2)
  const double temper,       // Air temperature (°C)
  const double apress,       // Air pressure (°C)
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


#endif

