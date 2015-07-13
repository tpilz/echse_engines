
// Note: qx_avg and qx_end hold either simulated or observed values, depending on
// whether observation data are read. These are variables whose values are to be
// passed to downstream objects.
// The variables qx_***_sim always hold the simulated values. They can be used,
// for example, to analyze the model error at a gage even in the case where the
// model imports observation data.

if ((inputExt(qobs_avg) >= paramNum(obs_lbound)) &
    (inputExt(qobs_avg) <= paramNum(obs_ubound))) {
  // Use observed flow (identical values for avg and end)
  set_output(qx_avg)= inputExt(qobs_avg);
  set_output(qx_end)= inputExt(qobs_avg);
} else {
  // Use the simulated values
  set_output(qx_avg)= inputSim(qi_avg);
  set_output(qx_end)= inputSim(qi_end);
}

set_output(qx_avg_sim)= inputSim(qi_avg);
set_output(qx_end_sim)= inputSim(qi_end);
