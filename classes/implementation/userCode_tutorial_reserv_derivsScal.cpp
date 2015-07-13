
const double ZERO= 0.;
double dH, Qex;

// Flow height above crest
dH= paramFun(h, u[INDEX_V]) - paramNum(Hcrest);

// Outflow rate
if (dH > ZERO) {
  Qex= paramNum(a) * pow(dH, paramNum(b));
} else {
  Qex= ZERO;
}

// Derivative of storage volume with respect to time
dudt[INDEX_V]= inputSim(Qin) - Qex;

