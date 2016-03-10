#ifndef ET_ACT_H
#define ET_ACT_H

// include libraries (header file) of meteorological functions and constants needed within some functions
#include "meteo/meteo.h"
#include "meteo/meteo_const.h"
// include library (header file) of aerodynamic and surface resistances
#include "resistances.h"
// include header file with et functions
#include "et_approaches.h"


///////////////////////////////////////////////////////////////////////////////
// Master function for choosing a specific method to calculate actual 
// evapotranspiration by using a choice flag. Function returns the act. evapotransp.
// rate calculated by the method chosen in [m/s]. All possible input values have
// to be given when calling the function but dummy values can be used if not all
// input variables are needed for the desired method. The function internally
// checks the input and calculates other (hydro-)meteorological quantities needed
// for the method chosen. In contrast to et_pot() the approaches herein need
// information about the actual soil moisture state.
//
// Choices:
// 	Simple empirical relations:
// 		1: Makkink
// 	Process-based models:
// 		11: Penman-Monteith
// 		12: FAO reference evaporation
// 		13: Shuttleworth-Wallace
// <tobias.pilz@uni-potsdam.de>, FEB 2015
///////////////////////////////////////////////////////////////////////////////
double et_act (
// Parameters and variables needed specifically for et_act()
	const double wc_vol_top,							// Actual volumetric soil water content at topmost horizon (m3/m3)
	const double wc_vol_root,							// Actual volumetric soil water content of the root zone (m3/m3)
	const double wc_sat,									// Volumetric water content at saturation of the root zone (m3/m3)
	const double wc_pwp,									// Volumetric water content of the root zone at permanent wilting point (m3/m3)
	const double wc_res,									// Residual volumetric water content of the root zone (m3/m3)
	const double wc_etmax,								// Parameter giving the volumetric water content where et_act equals et_pot, typically wc_etmax / wc_fk = [0.5..0.8] (m3/m3)
	const double bubble,									// Bubbling pressure of the root zone (cm) or (hPa)
	const double pores_ind,								// Pore-size-index of the root zone (-)
	const double wstressmin,							// Threshold for water stress effect on resistance (begin of stomata closure) (cm) OR (hPa)
	const double wstressmax,							// Threshold for water stress effect on resistance (total stomata closure, wilting point) (cm) OR (hPa)
	const double par_stressHum,						// Parameter to calculate water vapour deficit stomatal conductance stress factor (hPa-1) - in WASA a value of 0.03
// Common meteorological variables
	double temper,												// Air temperature (°C)
	const double temp_min,								// Minimum temperature within time step (°C)
	const double temp_max,								// Maximum temperature within time step (°C)
  double glorad,       									// Downward short-wave radiation (W/m2)
	const double rhum,										// relative humidity (%)
	const double wind,										// wind speed (m/s)
	double apress,       									// Air pressure (hPa)
	const double sundur,									// Sunshine duration of current day (h)
	const double cloud,										// Cloudiness (%)
// Common site-secific parameters
	const double lat,											// Latitude (decimal degree)
	const double lon,											// Longitude of the location of interest (decimal degrees west of Greenwich, e.g. Greenwich: 0°, Berlin: 345°, New York: 75°)
	const double elev,										// Elevation above sea level (m)
// Specific meteorological quantities (commonly calculated internally)
	double radex,													// Incoming extraterrestrial radiation (i.e. at top of atmosphere) (W/m2)
	double glorad_max, 										// Downward short-wave radiation under clear (cloudless) sky (Wm-2)
	double H_net,													// net incoming ( (1-alb) * short-wave + long-wave) radiation (Wm-2)
	double H_soil,												// net incoming (short-wave + long-wave) radiation hitting the soil surface (Wm-2)
	double H_long,												// Net incoming long-wave radiation (Wm-2)
	double soilheat,											// Soil heat flux (Wm-2)
	double totalheat,											// heat conduction into soil AND plants due to physical and biochemical energy storage (Wm-2)
// Meteorological parameters
	const double h_tempMeas,							// height of temperature measurement above ground (m)
	const double h_humMeas,								// height of humidity measurement above ground (psychrometer) (m)
	const double h_windMeas,							// height of windspeed measurement above ground (m)
	const double emis_a,									// Coefficient a for calculating net emissivity (-)
	const double emis_b,									// Coefficient b for calculating net emissivity (-)
	const double fcorr_a,									// Coefficient a for calculating cloudiness correction factor (-)
	const double fcorr_b,									// Coefficient b for calculating cloudiness correction factor (-)
	const double radex_a,									// Angstrom coefficient (share of radex on glorad under clouds) (-)
	const double radex_b,									// Angstrom coefficient (radex_a+radex_b = share of radex on glorad under clear sky) (-)
	const double f_day,										// soil heat flux calculation: Fraction of net_rad over daytime (in case of sub-daily application) (-)
	const double f_night,									// soil heat flux calculation: Fraction of net_rad over nighttime (in case of sub-daily application) (-)
// Vegetation and land-cover parameters and variables
	const double crop_makk,     					// Crop-factor after Makkink (-)
	const double crop_faoref,     				// Crop-factor for FAO reference approach (-)
	const double cano_height,							// canopy height (m)
	const double lai,											// leaf area index (m2/m2)
	const double alb,											// albedo (-)
	const double ext,											// Canopy extinction coefficient, Beer's law (-) -  in original WASA model code set to 0.5
	const double res_leaf_min,						// Plant-specific minimum (i.e. no stress occurs) stomatal resistance of a single leaf (sm-1)
	const double soil_dens,								// bulk density of soil of topmost soil horizon (kg/m3)
	const double glo_half,								// Solar radiation at which stomatal conductance is half of its maximum (W/m2) - in WASA model a value of 100
	const double res_b,										// Mean boundary layer resistance (sm-1) - in SW and WASA a value of 25
	const double drag_coef,								// Effective value of mean drag coef. of vegetative elements (-) - in SG and WASA a value of 0.07
	const double rough_bare,							// Roughness length of bare substrate (m) - in SW, SG and WASA a value of 0.01
	const double eddy_decay,							// Eddy diffusivity decay constant (-) - by SW and SG set to 2.5 for agricultural crops
	const double rss_a,										// Empirical constant for calculation of soil surface resistance (-) - in original code set to 26
	const double rss_b,										// Empirical constant for calculation of soil surface resistance (-) - in original code set to -1
// Computational parameters
	const double doy,											// Current day of year
	const int hour,												// hour of day in local time (including daylight daving time), range of [0..23]
	const int utc_add,										// Deviation of local time zone from UTC; may vary over the year due to daylight saving time; range of [-12..14] (hours)
	const double na_val,									// NA value
	const unsigned int delta_t,						// time step length (s)
// Choice flags
	const int choice,											// flag which method to use
	const int ch_rcs,											// Flag: canopy stomatal resistance by up-scaling of single leaf resistance, 1: Shuttlewort & Wallace (1985) eq. 19, 2: Saugier and Katerji (1991) eq. 4 used in WASA
	const int ch_roughLen, 								// Flag: Roughness length for momentum transfer, 1: SWAT manual (2011) eqs. 2:2.2.4, 2:2.2.5, 2: Shuttleworth & Gurney (1990) eq. 43
	const int ch_plantDispl,							// Flag: Displacement height for a plant, 1: SWAT manual (2011) eq. 2:2.2.7, 2: Shuttleworth & Gurney (1990) eq. 42
	const int ch_gloradmax								// Flag: calculation of maximum incoming short-wave radiation (clear sky), 1: Angstroem, 2: Allen (2005), ASCE standard etp, eq. 19 (based on elevation)
) {
	
// PRE-PROCESSING FOR ALL APPROACHES
	// air temperature
	if (abs(temper - na_val) < 0.01) {
		// calculate temper from temp_min and temp_max (simply the average) if available
		if ( (abs(temp_max - na_val) < 0.01) || (abs(temp_min - na_val) < 0.01) ) {
			stringstream errmsg;
			errmsg << "Value for air temperature is missing!";
			except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
			throw(e);
		}
		
		temper = (temp_max + temp_min) / 2.;
	}
	
	// apress (calculate if not given)
	if (abs(apress - na_val) < 0.01) {
		apress = apress_simple(elev);
	}
	
	// initialize variable for etp
	double etp = -9999.;
	
	
	
	
// EMPIRICAL RELATIONS, choice = {1..10}
	if(choice <= 10) {
		
		// glorad (try to calculate if not available)
		if (abs(glorad - na_val) < 0.01) {
			// cannot be calculated for sub.daily time steps
			if (delta_t < 86400) {
				stringstream errmsg;
				errmsg << "Value for incoming short-wave radiation is missing and cannot be calculate because of sub-daily resolution!";
				except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
				throw(e);
			}
			
			// calculate radex if not given
			if (abs(radex - na_val) < 0.01) {
				radex = rad_extraterr_daily(doy,lat);
			}
			
			// calculate glorad
			glorad = calc_glorad(radex,sundur,cloud,lat,doy,radex_a,radex_b);
		}
		
		// CALC ET
		
		//Makkink
		if(choice == 1) {
			if (abs(crop_makk - na_val) < 0.01) {
				stringstream errmsg;
				errmsg << "Value for crop factor is missing!";
				except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
				throw(e); 
			}
			
			// calculate pot. et; eta calculated later on by soil moisture reduction function
			etp = et_pot_makkink(glorad,temper,apress,crop_makk);
			
		}
	}
	
	
	
	
// PHYSICALLY-BASED RELATIONS, choice = {11..20}
	if(choice > 10) {
		
		// PRE-PROCESSING
		
		// net incoming radiation (Wm-2)
		if (abs(H_net - na_val) < 0.01) {
			
			// glorad (try to calculate if not available)
			if (abs(glorad - na_val) < 0.01) {
				// cannot be calculated for sub.daily time steps
				if (delta_t < 86400) {
					stringstream errmsg;
					errmsg << "Value for incoming short-wave radiation is missing and cannot be calculate because of sub-daily resolution!";
					except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
					throw(e);
				}
				
				// calculate radex if not given
				if (abs(radex - na_val) < 0.01) {
					if(delta_t == 86400)
						radex = rad_extraterr_daily(doy,lat);
					else if(delta_t == 3600)
						radex = rad_extraterr_hourly(doy,lat,hour,utc_add,lon);
					else {
						stringstream errmsg;
						errmsg << "Extraterrestrial radiation needs to be calculated but this is currently only possible for daily (86400s) and hourly (3600s) time steps!";
						except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
						throw(e);
					}
				}
				
				// calculate glorad
				glorad = calc_glorad(radex,sundur,cloud,lat,doy,radex_a,radex_b);
			}
			
			// calculate long-wave radiation if not given
			if (abs(H_long - na_val) < 0.01) {
				// calc glorad_max if not given
				if (abs(glorad_max - na_val) < 0.01) {
					// calculate radex if not given
					if (abs(radex - na_val) < 0.01) {
						if(delta_t == 86400)
							radex = rad_extraterr_daily(doy,lat);
						else if(delta_t == 3600)
							radex = rad_extraterr_hourly(doy,lat,hour,utc_add,lon);
						else {
							stringstream errmsg;
							errmsg << "Extraterrestrial radiation needs to be calculated but this is currently only possible for daily (86400s) and hourly (3600s) time steps!";
							except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
							throw(e);
						}
					}
					glorad_max = calc_glorad_max(ch_gloradmax, radex, radex_a, radex_b, elev);
				}
				H_long = net_longrad(temper,rhum,glorad,glorad_max,emis_a,emis_b,fcorr_a,fcorr_b);
			}
			
			// for FAO reference et always assume an albedo of 0.23
			if (choice == 12)
				H_net = (1. - 0.23) * glorad + H_long;
			else
				H_net = (1. - alb) * glorad + H_long;
		}
		
		// radiation at soil surface (only for SW; set to H_net for other aproaches) (Wm-2)
		if (abs(H_soil - na_val) < 0.01) {
			if (choice == 13)
				H_soil = H_net * exp(-1. * ext * lai); // canopy extinction according to Beer's law -> net radiation at soil surface (Wm-2)
			else
				H_soil = H_net;
		}
		
		// soil heat flux (Wm-2)
		if (abs(soilheat - na_val) < 0.01) {
			double daynight = -9999.;
			if (delta_t < 86400) {
				// calc glorad_max if not given
				if (abs(glorad_max - na_val) < 0.01) {
					// calculate radex if not given
					if (abs(radex - na_val) < 0.01) {
						if(delta_t == 3600)
							radex = rad_extraterr_hourly(doy,lat,hour,utc_add,lon);
						else {
							stringstream errmsg;
							errmsg << "Extraterrestrial radiation needs to be calculated but this is currently only possible for daily (86400s) and hourly (3600s) time steps!";
							except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
							throw(e);
						}
					}
					glorad_max = calc_glorad_max(ch_gloradmax, radex, radex_a, radex_b, elev);
				}
				if (glorad_max > 1e-6)
					daynight = 1;
				else
					daynight = 0;
			}
			soilheat = soil_heatflux(H_soil, f_day, f_night, daynight, delta_t);
		} else {
			// for simplicity set totalheat equal to soilheat if the latter was given as input and the former is missing (currently hard to interpret and quantify 'totalheat')
			if (abs(totalheat - na_val) < 0.01)
				totalheat = soilheat;
		}
		
		// air density (kgm-3)
		double rho_air = densityMoistAir(apress, temper, rhum);
		
		// saturation vapor pressure of air (hPa)
		double ez_0 = satVapPress_overWater(temper);
		
		// vapor pressure of air (hPa)
		double ez = vapPress_overWater(temper, rhum);
		
		// slope of saturation vapor pressure curve (hPa/K)
		double delta = slopeSatVapPress(temper);
		
		// psychrometric constant (hPa/K)
		double gamma = psychroConst(temper, apress);
		
		// Latent heat of water evaporation (J/kg)
		double lambda = latentHeatEvap(temper) * 1000.;
		
		// aerodynamic resistance (sm-1)
		double r_a = res_aero(ch_plantDispl, ch_roughLen, h_windMeas,h_humMeas,h_tempMeas,wind,cano_height, rough_bare, lai, drag_coef);
		
		// calculate actual stomatal resistance (i.e. under stress) of leaves (sm-1) and stress factors (-)
		double cond_water = stress_soilwater(wc_vol_root, wc_sat, wc_res, bubble, pores_ind, wstressmin, wstressmax); // soil water stress factor
		double cond_vap = stress_humidity(ez_0-ez, par_stressHum); // vapour pressure deficit stress factor, TODO: deficit should be WITHIN canopy, this is above
		double cond_co2 = 1.; // CO2 stress factor limiting stomatal conductance, currently not implemented (as in WASA)
		double cond_temp = 1.; // Temperature stress factor limiting stomatal conductance, currently not implemented (as in WASA)
		double cond_rad = 1.; //  Radiation stress factor limiting stomatal conductance, currently not implemented (as in WASA)
		double res_leaf = res_stom_leaf(res_leaf_min, cond_rad, cond_co2, cond_temp, cond_vap, cond_water); // actual stomatal resistance
		
		// canopy resistance (sm-1)
		double r_c = res_cs(ch_rcs, lai, res_leaf, ext, glorad, glo_half);
		
		
		// CALC ET
		
		// Penman-Monteith
		if (choice == 11) {
			return(et_penmon(lambda,delta,H_net,soilheat,rho_air,ez_0,ez,gamma,r_c,r_a));
		}
		
		// FAO reference evapotranspiration for reference grass surface
		if (choice == 12) {
			// calculate etp for reference surface
			etp = et_penmon_ref(temper,wind,apress,delta,H_net,soilheat,ez_0,ez,h_windMeas,delta_t);
			
			// etp for current crop; eta calculated later on by soil moisture reduction function
			etp = etp * crop_faoref;
		}
		
		// Shuttleworth-Wallace
		if (choice == 13) {
			
			// check water content values
			if( (wc_vol_top < wc_res) || (wc_vol_top > wc_sat) ) {
				stringstream errmsg;
				errmsg << "Cannot calculate actual evapotranspiration, unreasonable value of water content at top soil!";
				except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
				throw(e);
			}
			
			// total heat flux (conduction into soil and energy storage in plants and atmophere below reference height; for SW model only)
			// for simplicity calculate similar to soilheatflux; TODO: better approaches?
			if (abs(totalheat - na_val) < 0.01) {
				double daynight = -9999.;
				if (delta_t < 86400) {
					// calc glorad_max if not given
					if (abs(glorad_max - na_val) < 0.01) {
						// calculate radex if not given
						if (abs(radex - na_val) < 0.01) {
							if(delta_t == 3600)
								radex = rad_extraterr_hourly(doy,lat,hour,utc_add,lon);
							else {
								stringstream errmsg;
								errmsg << "Extraterrestrial radiation needs to be calculated but this is currently only possible for daily (86400s) and hourly (3600s) time steps!";
								except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
								throw(e);
							}
						}
						glorad_max = calc_glorad_max(ch_gloradmax, radex, radex_a, radex_b, elev);
					}
					if (glorad_max > 1e-6)
						daynight = 1;
					else
						daynight = 0;
				}
				totalheat = soil_heatflux(H_net, f_day, f_night, daynight, delta_t);
			}
			
			// bulk boundary layer resistance of vegetative elements in the canopy (sm-1)
			double r_ca = res_ca(lai,res_b);
			
			// calculate displacement length of vegetation (m)
			double d = height_plantDisplace(ch_plantDispl, cano_height, lai, drag_coef);
			
			// roughness length of vegetation (m)
			double z = roughLen_mom(ch_roughLen, cano_height, lai, d, drag_coef, rough_bare);
			
			// Aerodynamic resistance between soil surface and canopy height (sm-1)
			double r_sa = res_sa(wind, h_windMeas, cano_height, d, z, eddy_decay, rough_bare);
			
			// Aerodynamic resistance between canopy height and reference/measurement height (sm-1)
			double r_aa = res_aa(wind, h_windMeas, cano_height, d, z, eddy_decay);
			
			// surface resistance of substrate (sm-1)
			double theta_grav = wc_vol_top * 1000. / soil_dens; // density of water (set one) divided by bulk density of soil = maximum possible gravimetric water content
			double r_ss = res_ss(theta_grav, rss_a, rss_b);
			
			return(et_sw(lambda,delta,H_net,H_soil,totalheat,soilheat,rho_air,ez_0,ez,gamma,r_c,r_ca,r_ss,r_sa,r_aa));
		}
		
	}
	
	
	// check water content values
	if( (wc_pwp < 0.) || (wc_etmax > 1.) || (wc_pwp > wc_etmax) || (wc_vol_root < 0.) || (wc_vol_root > 1.) ) {
		stringstream errmsg;
		errmsg << "Cannot calculate soil moisture factor for calculation of actual evapotranspiration, unreasonable values of water content at root zone and/or wilting point and/or parameter wc_etmax!";
		except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		throw(e);
	}
	
	// calculate soil-moisture factor for empirical and FAO reference approaches
	double f_moist = min(1., max(0., (wc_vol_root - wc_pwp)/(wc_etmax - wc_pwp) ));
	
	// return eta
	return(etp * f_moist);
	
	
	
// ERROR IF NOT YET RETURNED
	stringstream errmsg;
	errmsg << "Could not calculate et_pot: Value of choice parameter not supported!";
	except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
	throw(e); 
}

#endif