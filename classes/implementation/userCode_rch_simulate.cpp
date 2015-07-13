double inflow_ini;
double inflow_end;
double k;
double v_new;

// We estimate the inflow rate at the beginning of the time step "inflow_ini"
// using the inflow rate at the end of the time step "inflow_end" and the
// time-step average "inflow_avg". This allows us to approximate the inflow
// rate as a linear function of time while keeping a proper the mass balance.
// The equation follows from the equality qi_avg = (qi_ini + qi_end) / 2.
// Note that the estimate may be negative. A value of zero is used in that
// case. Therefore, "inflow_end" needs to be re-estimated as well, to satisfy
// the mass balance.

inflow_ini= max(0., 2. * inputSim(qi_avg) - inputSim(qi_end));
inflow_end= 2. * inputSim(qi_avg) - inflow_ini;

// We compute the storage volume at the end of the time step using the
// analytical solution of the linear reservoir with linearily varying
// inflow. The applicable retention constant is estimated taking into
// the initial filling and the current inflow (weighted as in LARSIM).

k= (paramFun(v2k, stateScal(vol)) + paramFun(q2k, inflow_ini) +
   paramFun(q2k, inflow_end)) / 3.;
v_new= volume_end(stateScal(vol), delta_t, k, inflow_ini, inflow_end);


// Set the output variables
// (1) Average outflow rate (--> from mass balance)
set_output(qx_avg)= (stateScal(vol) + (inflow_ini+inflow_end)/2. * delta_t -
  v_new) / delta_t;
// (2) Outflow rate at end of step (--> from final volume)
set_output(qx_end)= 1. / k * v_new;

// Update storage volume
set_stateScal(vol)= v_new;
