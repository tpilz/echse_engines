
#ifndef RESERV_LINEAR_H
#define RESERV_LINEAR_H

////////////////////////////////////////////////////////////////////////////////
// Intregration of linear reservoir equation with constant inflow
// Returns: New storage volume (m3/s) after time step of length dt
inline double v_new (
  const double v_ini,  // Initial volume (m3)
  const double k,      // Storage constant (s)
  const double dt,     // Length of time step (s)
  const double q_in    // Inflow rate (m3/s)
) {
  return ( (v_ini - q_in * k) * exp(-dt / k) + q_in * k );
}

////////////////////////////////////////////////////////////////////////////////
// Analytical solution of linear reservoir with linearly variable inflow
// Returns the storage volume of linear reservoir at the end of a time step.
// The inflow is assumed to vary linearly with time (known initial and final
// inflow rate).
// This is the analytical solution of the underlying ODE:
//
//   dV / dt = qi_ini + (qi_end - qi_ini) * t - 1/k * V
//             \_____________ ______________/   \__ __/
//                           V                     V
//           Inflow as linear function of time   Outflow of linear reservoir

double volume_end(
  const double v_ini,  // Storage volume at beginning of time step
  const double dt,     // Length of time step
  const double k,      // Retention constant
  const double inflow_ini, // Inflow rate at beginning of time step
  const double inflow_end) // Inflow rate at end of time step
{
  double e= exp(-dt/k);
  double a= (inflow_end - inflow_ini) / dt;
  double b= -1/k;
  return(v_ini*e + a/pow(b,2)*(e-1) + inflow_ini/b*(e-1) - a/b*dt);
}

/*
# Analytical solution of linear reservoir with linearly variable inflow
#
# See http://www.tutorvista.com/math/separable-differential-equation
# for the integration method to be used.

#-----------------------------------------------------------
# Some test code in R
v_new= function(v_ini, dt, k, qi0, dqidt) {
  e1= exp(-dt/k)
  e2= e1 - 1
  a= dqidt
  b= -1/k
  c= qi0
  return(v_ini*e1 + a/b^2*e2 + c/b*e2 - a/b*dt)
}
v_ini= 0
k= 1.5
qi0=10
dqidt= -0.5
t= c(0.001,seq(from=1, to=10, by=1))
plot(t, v_new(v_ini, t, k, qi0, dqidt)/k, type="p")
# Solution used by Maniak und and in LARSIM (for outflow rate)
qx= function(qx0, dt, k, qi0, dqidt) {
  c1= 1 - exp(-dt/k)
  c2= 1 - c1 * (k/dt)
  return(qx0 + (qi0 - qx0) * c1 + (qi0+dqidt*dt - qi0) * c2)
}
lines(t, qx(qx0=v_ini/k, dt=t, k=k, qi0=qi0, dqidt=dqidt), lty=3)
#-----------------------------------------------------------
*/

#endif

