
// Save initial volume for later use in mass balance
const double v_ini= stateScal(V);

// Numerical ODE integration
odesolve_nonstiff(
  stateScal_all(),    // Initial value(s) of state variable(s)
  delta_t,            // Length of time step
  1.e-08,             // Accuracy
  1000,               // Maximum number of sub-steps
  this,               // Pointer to this object
  set_stateScal_all() // New value(s) of state variable(s)
);

// Set outputs
//   (1) Average outflow rate from balance
set_output(Qex)= inputSim(Qin) - (stateScal(V)-v_ini) / delta_t;
//   (2) Water level at end of time step
set_output(dH)= paramFun(h, stateScal(V)) - paramNum(Hcrest);
