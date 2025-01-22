#pragma once
#include "gas_ex_species_param.h"
#include "initinfo.h"

/*! \class CGasExchange 
* \brief Class for gas exchange calculations\n
* This class simulates gas exchange in plant leaves for C3 and C4 plants. \n
* \par Usage
	- Use <b>SetParams</b> to initialize the GasExchange object with parameters for a specific variety of plant
	- Use <b>SetVal</b> to pass environmental variables for a simulation and return a structure with output. \n
* See \ref Interface for details on input and output	
    

	*/

 
class CGasExchange  
{
public:
	CGasExchange(const TGasExSpeciesParam& photoparam);
	~CGasExchange(void);

	//Z sets input values for calculations for a particular set of environmental variables
	void SetVal(double PhotoFluxDensity, const TInitInfo info, double Tair, double CO2,
		double RH, double wind, double Press, double width, double psil, double lfNContent, double BiomassLeafOver);
	
	//Z These functions return results of calculations 
	
	double get_ANet() { return A_net; }						//Z return net photosynthesis (umol CO2 m-2 s-1) 
	double get_AGross() { return A_gross; }					//Z return gross photosynthesis  (umol CO2 m-2 s-1)
	double get_Transpiration() { return ET; }				//Z return transpiration rate (umol H2O m-2 s-1)
	double get_LeafTemperature() { return Tleaf; }			//Z return leaf temperature (C)
	double get_Ci() { return Ci; }							//Z return internal CO2 concentration (umol mol-1)
	double get_StomatalConductance() { return gs; }			//Z return stomatal conductance to water vapor (mol m-2 s-1)
	double get_BoundaryLayerConductance() { return gb; }    //Z return boundary layer conductance (mol m-2 s-1)
	double get_VPD() { return VPD; }						//Z return vapor pressure deficit (kpa)
	double get_leafpEffect() { return this->psileaf_stress; }

	double A_gross; //gross photosynthesis (umol CO2 m-2 s-1) 
	double A_net;	//net photosynthesis (umol CO2 m-2 s-1) 
	double ET;		//transpiration rate (umol H2O m-2 s-1) 
	double Tleaf;	//leaf temperature (C) 
	double Ci;		//internal CO2 concentration (umol mol-1) 
	double gs;		//stomatal conductance to water vapor (mol m-2 s-1) 
	double gb;      //boundary layer conductance (mol m-2 s-1) 
	double Rdc;		//Plant dark respiration umol m-2 s-1. 
	double VPD;		//vapor pressure deficit (kpa) 
	double temp;

private:
	// variables passed as arguments to the constructor
	TGasExSpeciesParam sParms;
	enum gsModel { BBW, L, H }; //Z the model for stomata conductance

	double  PhotoFluxDensity;	// Photosynthetic Flux Density (umol photons m-2 s-1) 
	double 	R_abs;				// Absorbed incident radiation (watts m-2)        
	double 	Tair;				// Air temperature at 2m, (C) 
	double 	CO2;				// CO2 concentration (umol mol-1 air) 
	double 	RH;					// Relative Humidity (%, i.e., 80) 
	double 	wind;				// Windspeed at 2 meters (m s-1) 
	double 	width;				// Leaf width (cm), when using leaf width, can find width/100.0, converting cm to m 
	double 	Press;				// Air pressure (kPa) 
	double 	Theta;				// Initial slope of CO2 response(umol m2 s - 1) - de Pury(1997)
	double 	psileaf;			// leaf water potential, MPa
	double 	psileaf_stress;		// 0 to 1 factor for stomatal closure effect of leaf water potential, i.e., effect of leaf water pressure (Mpa) on stomatal conductance.
	double 	lfNContent;			//YY lfNContent is leaf nitrogen content in unit of g m-2(leaf)
	double 	leaf_age;
	double  BiomassLeafOver;

	double Ci_Ca;				//Z Ratio of internal to external CO2, unitless
	double errTolerance;		//Z error tolerance for iterations
	double eqlTolerance;		//Z equality tolerance
	int    iter_total;			//Z holds total number of iterations
	int    iter1, iter2;		//Z holds iteration counters
	int    iter_Ci;				//Z iteration value for Ci umol mol-1, internal CO2 concentration
	bool   isCiConverged;		//Z true if Ci (internal CO2 concentration) iterations have converged

	// Main module to calculate gas exchange rates
	// produce photosynthesis and leaf temperature
	// Incorporate water stress effect
	void GasEx_psil(double psileaf, const TInitInfo info);
	//Z C3 photosynthesis calculation
	void Photosynthesis(double Ci, double psileaf, const TInitInfo info);  
	//Z calculates leaf temperature and transpiration
	void EnergyBalance();     

	//Z called iterively to find optimal internal CO2 concentration returns optimal internal CO2 concentration (CO2i) 
	//  secant search to find optimal internal CO2 concentration
	double SearchCi(double CO2i, double psileaf, const TInitInfo info); 

	//Z Calls photosynthesis modules to evaluate Ci dependent calculations returns difference between old Ci and new Ci
	//  Calculates a new value for Ci for the current values of photosynthesis and stomatal conductance
	//  used in "SearchCi" to fulfull the Ci iteration
	double EvalCi(double Ci, double psileaf, const TInitInfo info);   

	//Z Stomatal conductance (mol m-2 s-1)
	//  stomatal conductance for water vapor in mol m-2 s-1
	//  need a parameter to convert to CO2
	double gsw(double pressure, const TInitInfo info);   
	
	//Z Conductance for turbulant vapor transfer in air - forced convection (mol m-2 s-1)
	//  boundary layer conductance to vapor
	double gbw();           

	//Z Saturated vapor pressure at given temperature. kPa
	double Es(double Temperature);      
	//Z Slope of the vapor pressure curve, first order derivative of Es with respect to T
	//  sometimes present as "gamma" 
	double Slope(double Temperature);  

	//Z mitochondrial respiration
	double Rd();

	//Z Reduction in stomatal conductance using hourly bulk leaf water potential in MPa
	double set_PSIleafeffect(double pressure, const TInitInfo info);

	inline double Square(double a) { return a * a; } /*!< Squares number */ 
	inline double Min(double a, double b, double c) {return (__min(__min(a,b),c));} /*!< Finds minimum of three numbers */
	double QuadSolnUpper(double a, double b, double c); //Z Upper part of quadratic equation solution
	double QuadSolnLower(double a, double b, double c); //Z Lower part of quadratic equation solution

};


 
