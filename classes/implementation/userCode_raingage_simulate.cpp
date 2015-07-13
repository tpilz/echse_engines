
set_output(precip)= max(0., inputExt(precip_resid) + inputExt(precip_slope) *
  paramNum(elev) + inputExt(precip_inter));


