// Updated Gas Exchange Class - Soo Kim, Updated Dennis Timlin

/*!  Coupled model of photosynthesis-stomatal conductance-energy balance for a maize leaf this unit simulates Maize leaf gas-exchange characteristics
  including photosynthesis, \n traspiration, boundary and stomatal conductances,
  and leaf temperature based \n on von Caemmerer (2000) C4 model, BWB stomatal
conductance (1987) and \n Energy balance model as described in Campbell and Norman (1998)

Photosynthetic parameters were calibrated with PI3733 from
SPAR experiments at Beltsville, MD in 2002.

For potato, parameters come from Potato data in indoor chambers at Beltsville, MD in 2003
Stomatal conductance parameters were not calibrated

@authors Soo-Hyung Kim, Univ of Washington \n Dennis Timlin, USDA-ARS \n  David Fleisher, USDA-ARS \n
@version 1.0
@date August 2013

@note <b>-Bibliography </b>\n
Kim, S.-H., and J.H. Lieth. 2003. A coupled model of photosynthesis, stomatal conductance and transpiration for a rose leaf (Rosa hybrida L.). Ann. Bot. 91:771–781. \n
Kim, S.-H., D.C. Gitz, R.C. Sicher, J.T. Baker, D.J. Timlin, and V.R. Reddy. 2007. Temperature dependence of growth, development, and photosynthesis in maize under elevated CO2. Env. Exp. Bot. 61:224-236. \n
Kim, S.-H., R.C. Sicher, H. Bae, D.C. Gitz, J.T. Baker, D.J. Timlin, and V.R. Reddy. 2006. Canopy photosynthesis, evapotranspiration, leaf nitrogen, and transcription  \n
*/
// ws

//Z simplied by Z for cereal rye, 10/19/2023

#include "stdafx.h"
#include "gas_exchange.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>

// General fixed parameters
#define R 8.314			//ideal gas constant
#define maxiter 200		//maximum number of iterations
#define epsilon 0.97	//emissivity (Cambpell and Norman, 1998), pg 163
#define sbc 5.6697e-8	//stefan-boltzmann constant W m-2 k-4 - actually varies somewhat with temperatur
#define O 205.0			//gas units are mbar
#define Q10 2.0			//Q10 factor, photosynthesis activity increments for every 10 K
#define cPFD 4.6        //conversion factor from PAR (W/m2) to PFD (umol m-2 s) for solar radiation, see Campbell and Norman (1998) p 149
#define MAX_N_PCT 3.5
#define MIN_N_PCT 0.35

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

	
//Z the constructor is simplied from maizsim version (C4 plant) for a unified IO stream comparing to maizsim
//  the simplification occurs when substituting C4 functions by C3 plant functions, like the one for glysim
CGasExchange::CGasExchange(const TGasExSpeciesParam& photoparam)
{
	lfNContent = MIN_N_PCT;
	isCiConverged = false;
	errTolerance = 0.001;
	eqlTolerance = 1.0e-6;
	//Z initialize the photosynthesis parameters
	sParms = photoparam;
}


CGasExchange::~CGasExchange()
{
}

//Z set the weather condition to the photosynthesis model
//  include "TInitInfo" object
//  include the leaf nitrogen content as the input variable
//  remove leaf ET supply

/**Sets environment variables for a single execution of the module
	* Calls GasEx_psil() to calculate photosynthetic rate and stomatal conductance.
	* @param[in] PhotoFluxDensity	Photosynthetic Flux Density (umol Quanta m-2 s-1) (check)
	* @param[in] Tair	Air Temperature (C)
	* @param[in] CO2	CO2 concentration of the air (umol mol-1)
	* @param[in] RH	    Relative Humidity (%)
	* @param[in] wind	Windspeed at 2.5 m, m s-1
	* @param[in] Press	Atmospheric pressure (kpa m-2)
	* @param[in] ConstantTemperature boolian if true, leaf temperature=air temperature when calculating gas exchange
	\return nothing
*/
void CGasExchange::SetVal(double PhotoFluxDensity, const TInitInfo info, double Tair, double CO2, 
		                      double RH, double wind, double Press, double width,
		                      double psil, double lfNContent, double BiomassLeafOver)
{
	/***********************************/
	/* Sets initial parameters and runs routines
	/* for the coupled biochemical C3 gas exchange, stomatal conductance, energy balance routine
	/* orignally coded by Soo Kim
	/***********************************/

	//Z photosynthesis in both quantum flux PFD (umol m-2 s) or par energy PAR (W/m2)
	//  assign PFD value
	this->PhotoFluxDensity = PhotoFluxDensity;
	//  PAR is watts m-2
	double PAR = PhotoFluxDensity / cPFD; 
	//  If total solar radiation unavailable, assume NIR the same energy as PAR waveband from the prespective of "energy density"
	//  PAR - photosynthesis active radiation
	//  NIR - near infrared radiation
	double NIR = PAR; 
	
	//Z leaf long wave radiation
	//  times 2 for projected area basis, leaf has two surfaces
	this->R_abs = (1 - sParms.scatt) * PAR + 0.15 * NIR + 2.0 * (epsilon * sbc * pow(Tair + 273.15, 4.0));
    // shortwave radiation (PAR (=0.85) + NIR (=0.15) solar radiation absorptivity of leaves: =~ 0.5
	
	//Z transfer variables to local scope
	this->CO2 = CO2;
	this->RH = __min(100.0, __max(RH, 25.0)) / 100.0;
	this->Tair = Tair;
	this->width = width;			//leaf width in cm
	this->wind = wind;				//m s-1
	this->Press = Press;			//kPa
	this->psileaf = psil;			//Mpa or bar leaf water potential
	this->lfNContent = lfNContent;	//leaf N content in g N/m^2 leaf
	this->BiomassLeafOver = BiomassLeafOver;

	GasEx_psil(psil, info);   // Gas exchange calculations here
}

//Z carries out calculations for photosynthesis and stomatal conductance.
void CGasExchange::GasEx_psil(double psileaf, const TInitInfo info)
{
	/***********************/
	/* Main looping routine for coupled models
	/*	Incorporate water stress effect
	/***********************/

	//Z previous leaf temperture (for iteration)
	double Tleaf_old;
	int iter = 1;
	iter_total = 0;
	Tleaf = Tair; 
	Tleaf_old = 0.0;

	//Z a constant for initialize internal/intercellular CO2 partial pressure
	Ci = sParms.internalCO2Ratio * CO2;
	//Z bounday layer conductance
	gb = gbw();
	//Z stomatal conductance
	gs = gsw(psileaf, info);

	//DHF - Assume enzymatic activity negligable at this point, so no photosynthesis.  Need literature?  At what point does plant die?
	// Loomis and Conner indicate C3 leaves can go to 0C without permanent damage, not sure what temporary affect is on gas exchange rates is.
	if (Tair <= 2.0)
	{
		//Z low temperature for photosynthesis shutdown
		A_net = -Rd();
		A_gross = __max(Rd(), 0.0);
		ET = 0.0;
		return;
	}
	
	while ((abs(Tleaf_old -Tleaf)>0.01) && (iter < maxiter))
	{
		Tleaf_old=Tleaf;
		Ci=SearchCi(Ci, psileaf, info);
		gs=gsw(psileaf, info);
		EnergyBalance();
		iter2 =++iter; //iter=iter+1, iter2=iter; 
	}
}

//Z boundary layer conductance to vapor (mol H2O m^-2 s^-1)
double CGasExchange::gbw(void)
{
	//Z stomatal ratio: fraction of stomatal conductances of one side of the leaf to the other
	//  The stomatal ratio is the ratio of stomatal frequency on the adaxial surface to that on the abaxial surface.
	//  wheat, could be ~ 1.2: Fig 6. https://nph.onlinelibrary.wiley.com/doi/epdf/10.1111/nph.17563
	//const double stomaRatio = 0.5; // for soybean.

	//Z temporary holding variable for stomatal ratio calculations
	double ratio;
	//Z characteristic dimension of leaf
	double d;

	ratio = Square(sParms.stomaRatio + 1) / (Square(sParms.stomaRatio) + 1.0);

	//Z characteristic dimension of a leaf, leaf width is converted from cm to m using 100.0
	//Z note that we are computing leaf width based on the ryegrowth model, so this will be an input value
	d = width / 100.0 * 0.72;

	// wind is in m per second
	return (1.4 * 0.147 * sqrt(__max(0.1, wind) / d)) * ratio;
	// multiply by 1.4 for outdoor condition, Campbell and Norman (1998), p109, gva
	// multiply by ratio to get the effective blc (per projected area basis), licor 6400 manual p 1-9
}

//Z stomatal conductance for water vapor (mol H2O m^-2 s^-1)
//  calculates and returns stomatal conductance Uses Ball - Berry model.
double CGasExchange::gsw(double pressure, const TInitInfo info)
{
	double Pn,
		aa,		//aa, a value in quadratic equation 
		bb,		//bb, b value in quadratic equation 
		cc,		//cc, calcuation variable (x) in quadratic equation
		Ha,		//Ha, relative humidity
		Hs,		//hs, solution for relative humidity
		Cs,		//Cs, estimate of mole fraction of CO2 at the leaf surface
		Ca,		//Ca, air CO2 mole fraction
		gg,
		gamma,	//Gamma, CO2 compensation point in the absence of mitochondirial respiration, in ubar
		Ds;		//Ds, VPD at leaf surface 

	double temp = set_PSIleafeffect(pressure, info);				//Z temp = 1.00; no water stress
	gsModel myModel = BBW;											//Z need to put this in main program
	Ca = CO2;
	Ha = RH;

	//Z CO2 compensation point in the absence of Rd
	//Z need to check: Bykov O.D.; Koshkin V.A.; Catsky J Carbon dioxide compensation concentration of C3 and C4 plants: Dependence on temperature [wheat, bean, beet, sugar beet].
	//  A different equation: Fig.2 https://www.sciencedirect.com/science/article/pii/S0044328X83802037
	gamma = 36.9 + 1.88 * (Tleaf - 25.0) + 0.036 * Square(Tleaf - 25.0);			

	double P = Press / 100;
	Cs = Ca - (1.37 * A_net / gb) * P; // surface CO2 in mole fraction
	if (Cs == gamma) Cs = gamma + 1.0;
	if (Cs <= gamma) Cs = gamma + 1.0;

	
	gg = sParms.g0;
	if (A_net <= 0.0) { Pn = 0.00001; }
	else { Pn = A_net; }

	// Quadratic equation to obtain hs by combining StomatalConductance with diffusion equation
	//Z Hs=RH at leaf surface

	aa = temp * sParms.g1 * A_net / Cs;
	bb = gg + gb - aa;
	cc = (-Ha * gb) - gg;
	Hs = QuadSolnUpper(aa, bb, cc);		
	if (Hs > 1.0) { Hs = 1.0; }
	if (Hs < 0.0) { Hs = 0.0; }

	// VPD at leaf surface
	Ds = (1.0 - Hs) * Es(Tleaf);

	//Z compute gs from here
	if (A_net < 0) {
		return gg;
	}
	else
	{
		switch (myModel)//cccz: a good summary on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4242327/
		{
		case BBW:	//Ball-Berry-Woodrow model
			return (gg + temp * sParms.g1 * (A_net * Hs / Cs));
			break;
		case L:		//Leuning, 1995
			return gg + sParms.g1 * A_net / ((Cs - gamma) * (1 + Ds / sParms.g2));
			break;
		case H:
			return (gg + sParms.g1 * (A_net * Ha / Ca));
			break;
		default:
			return (gg + temp * sParms.g1 * (A_net * Hs / Cs));
		}
	}
}

//Z Reduction in stomatal conductance using hourly bulk leaf water potential in MPa
//  calculating effect of leaf water potential on stomatal conductance model of Tuzet et al. 2003 used  Yang 8/21/06
double CGasExchange::set_PSIleafeffect(double pressure, const TInitInfo info)   
{	
	if (pressure < -0.2)
	{
		this->psileaf_stress = __max((1.0 + exp(sParms.sf * sParms.phyf)) / (1.0 + exp(sParms.sf * (sParms.phyf - pressure))), 0.0);
	}
	else
	{
		this->psileaf_stress = 1.0;
	}
	
	double temp = this->psileaf_stress;
	return temp;
}

//Z secant search to find optimal internal CO2 concentration
//  will be a function of both CO2 and leaf water potential
double CGasExchange::SearchCi(double CO2i, double psileaf, const TInitInfo info)
{
	int iter = 0;
	double fprime, Ci1, Ci2, Ci_low, Ci_hi, Ci_m;
	double temp;
	Ci1 = CO2i;
	Ci2 = CO2i + 1.0;
	Ci_m = (Ci1 + Ci2) / 2.0;
	iter_Ci = 0;

	isCiConverged = true;

	do
	{
		iter++;
		//Secant search method
		if (abs(Ci1 - Ci2) <= errTolerance) { break; }
		if (iter >= maxiter)
		{
			isCiConverged = false;
			break;
		}
		//Z take derivative of EvalCi
		fprime = (EvalCi(Ci2, psileaf, info) - EvalCi(Ci1, psileaf, info)) / (Ci2 - Ci1);
		if (fprime != 0.0)
		{
			Ci_m = __max(errTolerance, Ci1 - EvalCi(Ci1, psileaf, info) / fprime);
		}
		else
		{
			Ci_m = Ci1;
		}
			
		Ci1 = Ci2;
		Ci2 = Ci_m;
		temp = EvalCi(Ci_m, psileaf, info);
//		double temp2 = maxiter;
	} while ((abs(EvalCi(Ci_m, psileaf, info)) >= errTolerance) || (iter < maxiter));

	// use bisectional search if above doesn't converge
	if (iter > maxiter)
	{
		Ci_low = 0.0;
		Ci_hi = 2.0 * CO2;
		isCiConverged = false;

		while (abs(Ci_hi - Ci_low) <= errTolerance || iter > (maxiter * 2))
		{
			Ci_m = (Ci_low + Ci_hi) / 2.0;
			if (abs(EvalCi(Ci_low, psileaf, info) * EvalCi(Ci_m, psileaf, info)) <= eqlTolerance) { break; }
			else if (EvalCi(Ci_low, psileaf, info) * EvalCi(Ci_m, psileaf, info) < 0.0) { Ci_hi = __max(Ci_m, errTolerance); }
			else if (EvalCi(Ci_m, psileaf, info) * EvalCi(Ci_hi, psileaf, info) < 0.0) { Ci_low = __max(Ci_m, errTolerance); }
			else { isCiConverged = false; break; }
		}
	}

	//Z make final assignment
	//Z it is a little bit weird to assign a value to argument of the function, i.e., CO2i
	//  but we can treat it as an temperoary holder for a double variable anyway
	CO2i = Ci_m;
	Ci_Ca = CO2i / CO2;
	iter_Ci = iter_Ci + iter;
	iter_total = iter_total + iter;
	return CO2i;
}

//Z Calculates a new value for Ci for the current values of photosynthesis and stomatal conductance
//  Determined from parameters from prior stem where energy balance was solved
double CGasExchange::EvalCi(double Ci, double psileaf, const TInitInfo info)
{
	/*
	* Called by SearchCi() to calculate a new value of Ci for the current values of photosynthesis and stomatal conductance
	* determined using parameters from a previous step where the energy balance was solved.
	* Estimate of internal CO2 concentration, umol mol-1
	* Return the difference between the passed value of Ci (old)and the new one.
	*/
	double newCi;
	Photosynthesis(Ci, psileaf, info);
	if (abs(gs) > eqlTolerance)
	{
		newCi = __max(1.0, CO2 - A_net * (sParms.SC_param / gs + sParms.BLC_param / gb) * (Press / 100.0));
	}
	else
	{
		newCi = __max(1.0, CO2 - A_net * (sParms.SC_param / eqlTolerance + sParms.BLC_param / gb) * (Press / 100.0));
	}
	return (newCi - Ci);
}

	

void CGasExchange::Photosynthesis(double Ci, double psileaf, const TInitInfo info)
//C3 photosynthesis
{
	const double curvature=0.999;		// curvature factor of Av and Aj colimitation
	const double theta = 0.7;			// curvature of electron transport in response to PAR
	double alpha, Kc, Ko, gamma, Ia, I2, Jmax, Vcmax, TPU, J, Av, Aj, Ap, Ac, Km, Ca, Cc, P, Tk;
	double f_age = 1.0;
	
	//Z Light response function parameters
	//  1. solar radiation => absorbed irradiance
	Ia = PhotoFluxDensity*(1-sParms.scatt);
	//  2. apparent quantum efficiency, params adjusted to get value 0.3 for average C3 leaf
	alpha = (1.0 - sParms.f) / 2.0;
	//  3. absorbed irradiance => useful light absorbed by PSII
	I2 = Ia * alpha;

	//Z leaf temperature and co2 compensation point
	Tk = Tleaf + 273.15;
	gamma = 36.9 + 1.88 * (Tleaf - 25) + 0.036 * Square(Tleaf - 25); //CO2 comp point in absense of mito respiration, in ubar

	//Z initialize net assimilation
	A_net = 0;

	//Z other input parameters and constants
	P = Press / 100.0;
	Ca = CO2 * P; // conversion to partial pressure Atmospheric partial pressure of CO2, kPa
	
	//Z Arrhenius temperature adjustment of the RuBisCO reaction speed, adding C atom or oxidiation
	Kc = sParms.Kc25 * exp(sParms.Eac * (Tleaf - 25.0) / (298.15 * R * Tk));
	Ko = sParms.Ko25 * exp(sParms.Eao * (Tleaf - 25.0) / (298.15 * R * Tk));
	Km = Kc * (1.0 + O / Ko); // effective M-M constant for Kc in the presence of O2

	//Z adding Nitrogen stress using exponential function 
	double CriticalNitrogen;
	CriticalNitrogen = __max(MIN_N_PCT, lfNContent);
	double N_effect = 2.0 / (1.0 + exp(-2.9 * (CriticalNitrogen - MIN_N_PCT))) - 1.0;
	if (CriticalNitrogen > 1.5750 / 3.5 * MAX_N_PCT) { N_effect = 1.0; }
	N_effect = __max(N_effect, 0.001);
	double Vcm25_L = sParms.Vcm25 * N_effect;
	double Jm25_L = sParms.Jm25 * N_effect;
	double TPU25_L = sParms.TPU25 * N_effect;

//	//Z testing non-N stress
//	Vcm25_L = sParms.Vcm25;
//	Jm25_L = sParms.Jm25;
//	TPU25_L = sParms.TPU25;

	//Z Arrhenius or Modified Arrhenius function for temperature effects
	// de Pury 1997
	Jmax = f_age * Jm25_L * __max(exp(((Tk - 298.15) * sParms.Eaj) / (R * Tk * 298.15)) *
		(1.0 + exp((sParms.Sj * 298.15 - sParms.Hj) / (R * 298.15))) /
		(1.0 + exp((sParms.Sj * Tk - sParms.Hj) / (R * Tk))), 0.5);
	// Used peaked response, DHF
	Vcmax = f_age * Vcm25_L * __max(exp(((Tk - 298.15) * sParms.EaVc) / (R * Tk * 298.15)) *
		(1.0 + exp((sParms.Sv * 298.15 - sParms.Hv) / (R * 298.15))) /
		(1.0 + exp((sParms.Sv * Tk - sParms.Hv) / (R * Tk))), 0.5);
	TPU = f_age * TPU25_L * __max(exp(sParms.EaVp * (Tleaf - 25) / (298.15 * R * Tk)), 0.5);

	//Z testing sugar effects
	Vcmax = Vcmax * exp(-10.0 * __max(this->BiomassLeafOver - 0.15, 0.0));
	Jmax = Jmax * exp(-10.0 * __max(this->BiomassLeafOver - 0.15, 0.0));

	//assume infinite gi (stomata is sooooo conductive)
	Cc = Ci;
	gs = gsw(psileaf, info);
	gb = gbw();

	//Z comes the three photosynthesis limiting speeds
	//  that is the critical assumption of Farquhar von Caemmerer Berry model
	Av = (Vcmax * (Cc - gamma)) / (Cc + Km);
	J = (((alpha * Ia + Jmax) - sqrt(Square(alpha * Ia + Jmax) - 4 * alpha * Ia * (Jmax)*theta)) / (2 * theta));
	Aj = J * (Cc - gamma) / (4 * (Cc + 2 * gamma));
	Ap = 3 * TPU;
	Ac = ((Av + Aj) - sqrt(Square(Av + Aj) - 4 * curvature * Av * Aj)) / (2 * curvature); //curvature account for collimitation between Av and Aj

	if (Cc > gamma)
	{
		A_net = __min(Ac, Ap) - Rd();
	}
	else
	{
		A_net = Av - Rd();
	}
	A_gross = __max(A_net + Rd(), 0.0);
	gs = gsw(psileaf, info);
}

//Z see Campbell and Norman (1998) pp 224-225
//  because Stefan-Boltzman constant is for unit surface area by denifition,
//  all terms including sbc are multilplied by 2 (i.e., gr, thermal radiation)
void CGasExchange::EnergyBalance()
	/* 
	Calculates Transpiration rate (T) and leaf temperature (Tleaf). Iterates by recalculating photosynthesis until leaf temperatures converge
	See Campbell and Norman (1998) pp 224-225
    Because Stefan-Boltzman constant is for unit surface area by denifition,
	all terms including sbc are multilplied by 2 (i.e., RadiativeConductance, thermal radiation)
		
	Return nothing but calculates transpiration (T) and leaf temperature (Tleaf)
	*/
{
	const long lambda = 44000;		//latent heat of vaporization of water J mol-1 at 25C
	const double psc = 6.66e-4;		//psycrometric constant units are C-1
	const double Cp = 29.3;			// thermodynamic psychrometer constant and specific hear of air, J mol-1 C-1
	double gha, gv, gr, ghr, psc1, Ea, thermal_air, Ti, Ta;
	double lastTi, newTi;
	int iter;

	Ta = Tair;
	Ti = Tleaf;
	//gha = gb*(0.135/0.147);  // heat conductance, gha = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Note: this was only true if stomatal ratio = 1
	gha = 1.4 * 0.135 * sqrt(__max(0.1, wind) / (width / 100.0 * 0.72));
	gv = gs * gb / (gs + gb);
	gr = (4.0 * epsilon * sbc * pow(Ta + 273.15, 3) / Cp) * 2.0;	// radiative conductance, 2 account for both sides
	ghr = gha + gr;
	thermal_air = epsilon * sbc * pow(Ta + 273.15, 4) * 2.0;		// emitted thermal radiation
	psc1 = psc * ghr / gv;											// apparent psychrometer constant
	VPD = Es(Ta) * (1 - RH);										// vapor pressure deficit
	Ea = Es(Ta) * RH;												// ambient vapor pressure

	//iterative version
	newTi = -10.0;
	iter = 0;
	lastTi = Tleaf;
	double Res, dRes;
	double thermal_leaf;
	while ((abs(lastTi - newTi) > 0.001) && (iter < maxiter))
	{
		lastTi = newTi;
		Tleaf = Ta + (R_abs - thermal_air - lambda * gv * VPD / Press) / (Cp * ghr + lambda * Slope(Tair) * gv);
		thermal_leaf = epsilon * sbc * pow(Tleaf + 273.15, 4) * 2.0;
		Res = R_abs - thermal_leaf - Cp * gha * (Tleaf - Ta) - lambda * gv * 0.5 * (Es(Tleaf) - Ea) / Press;
		dRes = -4 * epsilon * sbc * pow(273.15 + Tleaf, 3) * 2 - Cp * gha * Tleaf - lambda * gv * Slope(Tleaf);
		newTi = Tleaf + Res / dRes;
		iter++;
	}
	Tleaf = newTi;
	ET = __max(0.0, 1000.0 * gv * ((Es(Tleaf) - Ea) / Press) / (1 - (Es(Tleaf) + Ea) / (Press)));
	//1000 is to go from mol to mmol
}


//Z calculates and returns Saturation vapor pressure (kPa)
double CGasExchange::Es(double Temperature) 
{
	//Campbell and Norman (1998), p 41 Saturation vapor pressure in kPa
	return (0.611 * exp(17.502 * Temperature / (240.97 + Temperature)));
}

//Z Calculates the slope of the sat vapor pressure curve: 
//  first order derivative of Es with respect to T
double CGasExchange::Slope(double Temperature) 
{
	// slope of the sat vapor pressure curve: first order derivative of Es with respect to T
		// units of b and c are  degrees C
	const double b = 17.502; const double c = 240.97;
	return (Es(Temperature) * (b * c) / Square(c + Temperature) / Press);
}

//Z an isolated function for dark respiration
//  same theory as we did Vcm, J and TPU,
//  but for Rd, we just put another function for that
double CGasExchange::Rd()   //Should be expanded to include other env. and physiological factors
{
	//sParms.Ear: exponential rate of arrhenious function for mitochondrial respiration (J mol)
	return sParms.Rd25 * exp(sParms.Ear * (Tleaf - 25.0) / (298.15 * R * (Tleaf + 273.15)));
}
	

//Z These two functions solve the quadratic equation.
// NOT designed for linear equations
double CGasExchange::QuadSolnUpper(double a, double b, double c)
{
	/* solves the uppper part of the quadratic equation ax2+bx2=c

	@param[in] a
	@param[in] b
	@param[in] c 

	\return lower portion of x
	*/
	if (a == 0.0) { return 0.0; }
	else if ((b * b - 4.0 * a * c) < 0.0) { return -b / a; }   //imaginary roots
	else { return (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a); }
}

double CGasExchange::QuadSolnLower(double a, double b, double c)
{
	/** solves the lower part of the quadratic equation ax2+bx=c

	@param[in] a
	@param[in] b
	@param[in] c
	\return lower portion of x
	*/
	if (a == 0.0) { return 0.0; }
	else if ((b * b - 4.0 * a * c) < 0.0) { return -b / a; }   //imaginary roots
	else { return (-b - sqrt(b * b - 4.0 * a * c)) / (2.0 * a); }
}

