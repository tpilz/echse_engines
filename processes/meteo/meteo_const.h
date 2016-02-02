
#ifndef METEO_CONST_H
#define METEO_CONST_H

////////////////////////////////////////////////////////////////////////////////
// Define constants frequently used for meteorological calculations.
// <tobias.pilz@uni-potsdam.de>, JUL 2015
////////////////////////////////////////////////////////////////////////////////

// conversion term for temperature [°C] -> [K]
const double T_DEG_K = 273.15;

// Stefan Boltzmann constant [Wm-2K-4]: https://de.wikipedia.org/wiki/Stefan-Boltzmann-Gesetz, Dyck & Peschke (1995), p. 31
const double SIGMA = 5.670373e-8;

// Pi
const double PI = 3.141592653589793;

// solar constant in (Wm-2) as measured by Kopp & Lean (2011); revised value, lower than previous estimates!
const double SOLAR_C = 1360.8;

// molar mass of dry air (kgmol-1), Picard et al. (2008)
const double M_DA = 28.96546e-03;

// molar mass of water (kgmol-1), Picard et al. (2008)
const double M_W = 18.01528e-03;

// molar gas constant (Jmol-1K-1), Mohr & Taylor (2005)
const double R = 8.314;

// Karman constant, SWAT manual (2011) citing Jensen et al. (1990); range of 0.36 to 0.43
const double KARMAN = 0.41;

// specific heat of moist air (Jkg-1K-1) under typical room conditions
// https://en.wikipedia.org/wiki/Heat_capacity
const double SPHEATMOIST = 1012.;

#endif
