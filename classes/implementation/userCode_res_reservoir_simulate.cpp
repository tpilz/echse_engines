////////////////////////////////////////////////////////////////////////////
// Select the q(h) relation which is in effect

// Number of rules is 12 (fixed size; some rules may be unused)
const unsigned int funcindices_h2q [12]= {h2q_1.index, h2q_2.index,
  h2q_3.index, h2q_4.index, h2q_5.index, h2q_6.index, h2q_7.index,
  h2q_8.index, h2q_9.index, h2q_10.index, h2q_11.index, h2q_12.index};

// Get current rule (specified in external time series)
unsigned int rule= static_cast<unsigned int>(round(inputExt(ctrl_rule)));
if ((rule < 1) | (rule > 12)) {
  stringstream errmsg;
  errmsg << "Index of h2q function is out of range. Value is" << rule <<
    " (rounded from " << inputExt(ctrl_rule) << ") but should be in 1...12.";
  except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
  throw(e);
}

// Select index of the appropriate function
T_index_paramFun index= {funcindices_h2q[rule-1]};

////////////////////////////////////////////////////////////////////////////
// Estimate the inflow rate at the beginning and end of a time step  
// "inflow_ini" and "inflow_end" taking into account the mass balance. See comments
// in the method for the reach class. For 4th order RK, we also need the
// rate at 1/2 the time step.

double inflow_ini= max(0., 2 * inputSim(qi_avg) - inputSim(qi_end));
double inflow_end= 2 * inputSim(qi_avg) - inflow_ini;
double inflow_mid= (inflow_ini + inflow_end) / 2;

////////////////////////////////////////////////////////////////////////////
// Compute new volume by 4th order Runge-Kutta
// The derivative is: ddt_volume = inflow - outflow

double s1,s2,s3,s4;
double y0,y1,y2,y3;
double vol_new;

y0= stateScal(vol);
s1= inflow_ini - paramFun(index,paramFun(v2h,y0)); // 1st temp. derivative
y1= y0 + s1 * delta_t/2.;
s2= inflow_mid - paramFun(index,paramFun(v2h,y1)); // 2nd temp. derivative
y2= y0 + s2 * delta_t/2.;
s3= inflow_mid - paramFun(index,paramFun(v2h,y2)); // 3rd temp. derivative
y3= y0 + s3 * delta_t;
s4= inflow_end - paramFun(index,paramFun(v2h,y3)); // 4th temp. derivative
vol_new= y0 + delta_t * (s1 + 2.*s2 + 2.*s3 + s4) / 6.;

////////////////////////////////////////////////////////////////////////////
// Set output variables

// (1) Average outflow rate (--> from mass balance)
set_output(qx_avg)= (stateScal(vol) + (inflow_ini+inflow_end)/2 * delta_t - vol_new) / delta_t;

// (2) Outflow rate at end of step (--> from final volume)
set_output(qx_end)= paramFun(index,paramFun(v2h,vol_new));

// Update volume
set_stateScal(vol)= vol_new;

