
#ifndef ROUTE_H
#define ROUTE_H

#include <cmath>


///////////////////////////////////////////////////////////////////////////////
// Routing, i.e. translation and retention, of water along a river / channel reach.
// Output is the river reach outflow for the current time step.
// Choices:
// 		1: Unit hydrograph approach as in WASA (Guentner, 2002)
// <tobias.pilz@uni-potsdam.de>, FEB 2017
///////////////////////////////////////////////////////////////////////////////

double chan_route (
	const double qin,			// inflow to the river reach, (m3/s)
// Unit hydrograph (WASA) approach specific
	const vector<double> &uh_f,	// unit hydrograph ordinate values (-)
	vector<double> &uh_q,		// outflow storage for upcomming time steps (m3/s)
// Choice flags
	const unsigned int choice	// Choice flag for a routing method; see comments above for implementations
) {
	double res = -9999.;
	
	switch(choice) {
		case 1: // Unit hydrograph, WASA
			
			// calculate output for this time step
			res = qin * uh_f[0] + uh_q[0];
			
			// calculate output for the upcomming time steps according to unit hydrograph
			for (unsigned int i=1; i<uh_f.size(); i++) {
				uh_q[i-1] = qin * uh_f[i] + uh_q[i];
			}
			
			break;
			
			
		default:
			stringstream errmsg;
			errmsg << "Invalid choice to calculate channel routing! Currently supported is one of {1}.";
			except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
			throw(e); 
	}
	
	return(res);
}



///////////////////////////////////////////////////////////////////////////////
// River transmission losses.
// Choices:
// 		1: Evaporation without infiltration as in WASA (Guentner, 2002)
// <tobias.pilz@uni-potsdam.de>, FEB 2017
///////////////////////////////////////////////////////////////////////////////

double trans_loss (
	const double q,				// River flow experiencing transmission losses (m3/s)
	const double pet,			// potential evapotranspiration within the corresponding subbasin (m/s)
	const double chan_len,		// River reach length (m)
	const unsigned int choice	// Choice flag for a transmission losses method; see comments above for implementations
) {
	double res = -9999.;
	double evap = -9999.;
	double w = -9999.;
	
	switch(choice) {
		case 1: // Evaporation without infiltration, river width according to relationship by Leopold (1994), WASA

			// calculate river width for a given river flow; global relationship by Leopold (1994) (m)
			w = pow(10., log10(q)*0.494+1.031);
			// calculate river evaporation from river surface and subbasin potential evapotranspiration (m3/s)
			evap = w * chan_len * pet;
			// subtract evaporation from river flow (m3/s)
			res = min(q, evap);
			
			break;
			
			
		default:
			stringstream errmsg;
			errmsg << "Invalid choice to calculate river transmission losses! Currently supported is one of {1}.";
			except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
			throw(e); 
	}
	
	return(res);
}

#endif