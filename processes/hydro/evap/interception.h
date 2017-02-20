
#ifndef INTERCEPTION_H
#define INTERCEPTION_H

////////////////////////////////////////////////////////////////////////////////
// Interception rate (m/s).
// Simple bucket model based on implementation in model WASA-SED rev. 94, 
// soilwat.f90 line 1073 ff.
// <tobias.pilz@uni-potsdam.de>, AUG 2015
////////////////////////////////////////////////////////////////////////////////

double intercept(
	const double prec,			// Precipitation rate (m/s)
	const double lai,				// Leaf Area Index (m2/m2)
	const double intcf,			// interception capacity per unit LAI (m) - calibration value
	const double intercept,	// interception storage filling state at beginning of timestep (m)
	const double delta_t		// time step length (s)
) {
	// check input
	if(intercept < -1e-12) {
		stringstream errmsg;
		errmsg << "Negative interception storage: intercept = " << intercept;
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e); 
	}
	
	// interception occurs only if there is precipitation
	if (prec > 0.) {
		// maximum possible interception
		double maxintc = lai * intcf;
		// calculate amount that is intercepted in this timesetp
		double inter_p = max(0., maxintc-intercept);
		return( min(prec, inter_p / delta_t ) );
	} else {
		return(0.);
	}
}

#endif
