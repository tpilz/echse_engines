
#ifndef RESERV_CASCADE_H
#define RESERV_CASCADE_H

// Rate of fast runoff reaching a single cascade (m3/s)
double qi_cascade_single;

// Temporary vectors for data exchange with ODE solver
const unsigned int CASCADELENGHT= 6;
vector<double> y_ini(CASCADELENGHT);  // Initial volumes
vector<double> y_new(CASCADELENGHT);  // Final volumes
// Utility variables in derivative computations
vector<double> hgt(CASCADELENGHT);    // Heights of water surface
vector<double> qx(CASCADELENGHT);     // Outflow rates

// Poleni formula for free overflow
// Restrictions/Assumptions:
// - Free overflow, i.e. no downstream submergence
// - Fixed empirical loss coefficient
double q_over(
  const double h,       // Height of water surface above reference
  const double hCrest,  // Height of dam crest above reference
  const double wCrest   // Flow widht
) {
  // 0.65 : Empirical loss coefficient depending on the crest's shape
  return(0.65 * 2./3. * wCrest * sqrt(2.*9.81) * pow(max(0., h-hCrest), 1.5));
}
/*
### R test code: ###############################################################
rm(list=ls())
pow= function(x,e) {x^e}
q_over= function(h, hCrest, wCrest) {
  return(0.65 * 2./3. * wCrest * sqrt(2.*9.81) * pow(max(0., h-hCrest), 1.5));
}
hCrest= 2
h= seq(from=0, to=2*hCrest, by=hCrest/100)
q= rep(NA, length(h))
for (i in 1:length(h)) q[i]= q_over(h[i], hCrest, wCrest= 0.5)
plot(h, q, type="l")
### End: R test code ###########################################################
*/

// Flow through a hole at the bottom of a dam
// Restrictions/Assumptions:
// - Free flow, i.e. no downstream submergence
// - Fixed empirical loss coefficient
// - Whole is located at the dam's bottom
// - If the hole falls (partly) dry, the flow rate is reduced in a smooth,
//   non-linear way using a power model similar to the Poleni formula
double q_hole(
  const double h,       // Height of water surface above reference
  const double dimOut   // Dimension of the outlet (width == height)
) {
  // 0.58 : Empirical loss coefficient for a whole with width == height
  if (h > dimOut) {
    return( 0.58 * pow(dimOut, 2.) * sqrt(2*9.81*(h-dimOut/2.)) );
  } else if (h > 0.) {
    return( 0.58 * pow(dimOut, 2.) * sqrt(2*9.81*(dimOut/2.)) * pow(h/dimOut, 1.5));
  } else {
    return( 0. );
  }
}

/*
### R test code: ###############################################################
rm(list=ls())
pow= function(x,e) {x^e}
q_hole= function(h, dimOut) {
  if (h > dimOut) {
    return( 0.58 * pow(dimOut, 2.) * sqrt(2*9.81*(h-dimOut/2.)) );
  } else if (h > 0.) {
    return( 0.58 * pow(dimOut, 2.) * sqrt(2*9.81*(dimOut/2.)) * pow(h/dimOut, 1.5));
  } else {
    return( 0. );
  }
}
dimOut= 0.1
h= seq(from=0, to=5*dimOut, by=dimOut/100)
q= rep(NA, length(h))
for (i in 1:length(h)) q[i]= q_hole(h[i], dimOut)
plot(h, q, type="l")
### End: R test code ###########################################################
*/

// Height of water level (h) as a function of storage volume (v) for a reservoir
// with
// - X-section area as a function of water level: a(h) = shapeFac * h ^ shapeExp
// - Linear bottom slope
// Note:
// - The expression h(v) is obtained by integrating a(h) over the length (L) of
//   the reservoir. The length L is related to h by s= h/L.
// - For small storage volumes lower than v_crit, a linear relation is assumed.
//   This is to make dh/dv less steep in that region. Otherwise, ODE solvers fail.
//   The value of v_crit must be > 0. Larger values speed up the numerical
//   integration of the storage equation - i.e. one should use the largest value
//   that still allows for a more or less accurate representation of the
//   reservoir's geometry.
double h_reservoir(
  const double v,        // Storage volume (m^3)
  const double shapeFac, // Parameter in eqn. for x-section area, see above
  const double shapeExp, // Parameter in eqn. for x-section area, see above
  const double s,        // Bottom slope of the reservoir (m/m)
  const double v_crit    // Critical volume (m^3) below which h(v) is linear
) {
  if (v > v_crit) {
    return( pow(v*(shapeExp+1.)*s/shapeFac, 1./(shapeExp+1.)) );
  } else {
    return( pow(v_crit*(shapeExp+1.)*s/shapeFac, 1./(shapeExp+1.)) * v/v_crit);
  }
}
/*
### R test code: ###############################################################
rm(list=ls())
pow= function(x,e) {x^e}
h_reservoir= function(v, shapeFac, shapeExp, s, v_crit) {
  if (v > v_crit) {
    return( pow(v*(shapeExp+1.)*s/shapeFac, 1./(shapeExp+1.)) );
  } else {
    return( pow(v_crit*(shapeExp+1.)*s/shapeFac, 1./(shapeExp+1.)) * v/v_crit);
  }
}
v= seq(0, 10, 0.001)
h= rep(NA, length(v))
for (i in 1:length(v)) h[i]= h_reservoir(v[i], shapeFac=1.78, shapeExp=2.09,
  s=0.01, v_crit=1)
plot(v, h, type="l", xlab="v", ylab="h")
### End: R test code ###########################################################
*/

#endif

