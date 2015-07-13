
// Save initial volume for later use in mass balance
const double v_ini= stateScal(v);

// Initialize volumes of precipitation and evaporation per time step to zero
set_stateScal(vp)= 0.;
set_stateScal(ve)= 0.;

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
// Average outflow rate (from water balance)
set_output(qx_avg)= inputSim(qi_avg) - (stateScal(v)-v_ini) / delta_t
  + stateScal(vp) / delta_t - stateScal(ve) / delta_t;

// Instantaneous outflow rate at end of time step
set_output(qx_end)= paramFun(h2q, paramFun(v2h,stateScal(v)));
// Water level at end of time step
set_output(h)= paramFun(v2h,stateScal(v));
