// Entry point of the RYESIM dll
// the structure is inherit from MAIZSIM

// myplant.cpp : Defines the entry point for the DLL application.
//#define MYPLANT_EXPORTS

#include "stdafx.h"
#include "RyeCrop.h"
#include "RyeController.h"

#include "time.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;
#include <cmath>
#define endl "\n"
#define comma ","

// note that we have to dereference the variable in order to assign a value that can be passed back to 2DSOIL.
// This is because the FORTRAN program expects a pointer rather than a value.
// I don't think this applies to structures as it does variables that may be in the arguments list.
// note use of lower case names. Upper and lower case conversions between fortran and C++ don't matter here because these are arguments and not a function name. 
// CROP must be upper case because it is a function name
int compare(const void* arg1, const void* arg2)
{
	/* Compare all of both strings: */
	if (*(double*)arg1 > *(double*)arg2) return 1;
	else if (*(double*)arg1 < *(double*)arg2) return -1;
	return 0;
};
#ifdef _WIN32
void _stdcall CROP(struct
#else
void crop(struct
#endif
	ShootCommon* SHOOTR,
	WeathCommon* Weather,
	GridCommon* grid_public,
	NodeCommon* node_public,
	BoundaryCommon* bound_public,
	TimeCommon* time_public,
	ModuleCommon* module_public,
	FileCommon* file_public
)

{
	//SK, declare as static to ensure only one copy is instantiated during 2DSOIL execution
	// varFile contains variety information, GraphicFile holds output,LeafFile holds individual leaves

	//Z weather and initial info will be fed into the model within controller
	//  other info will be fed or computed in the crop model

	static RyeController* pSC;

	// *************************** model initialization **************************************************
	
	//SK initialiing crop module
	//   This will only be called once
	if (time_public->lInput == 1) 
	{
		//Z this should be in the first call
		//  mark the index of the crop module in the whole (big) model structure 
		//  used to aligned the time steps.
		module_public->NumMod = module_public->NumMod + 1;
		ModNum = module_public->NumMod;

		//Parse the file names from the FORTRAN strings passed from 2dsoil
		//KY looks like GNU Fortran handle linebreak differently, making filename detection unusable
		//KY this new macro based on std::string should work on both platforms with smaller code
#define SETSTR(s, n) std::string s(n, sizeof(n)); s.erase(s.find_last_not_of(" \n\r\t")+1);
		SETSTR(varFile, file_public->VarietyFile);
		SETSTR(GraphicFile, file_public->GraphicsFile);
		SETSTR(LeafFile, file_public->LeafFileIn);

		//Z import data from the soil data structure
		initInfo.plantDensity = SHOOTR->PopArea;
		initInfo.latitude = Weather->LATUDE;
		initInfo.longitude = Weather->Longitude;
		initInfo.altitude = Weather->Altitude;
		initInfo.year = time_public->Year;
		initInfo.sowingDay = time_public->sowingDay;
		initInfo.beginDay = time_public->beginDay;
		initInfo.endDay = time_public->endDay;
		initInfo.timeStep = time_public->TimeStep; //unit in minute
		time_public->iTime = 1;
		ActualPopSlab_prev = 0;
		//Z to assign seed depth, should first know the y coordinate of the soil surface
		double maxY = 0.0;
		int GridNodeNum = grid_public->NumNP;
		for (int ii = 0; ii < GridNodeNum; ii++)
		{
			if (grid_public->y[ii] > maxY) maxY = grid_public->y[ii];
		}
		initInfo.seedDepth = maxY - SHOOTR->yBStem;

		SHOOTR->LCAI = 0.0;
		SHOOTR->LAREAT = 0.0;
		SHOOTR->Height = 0.0;
		//dt change for debugging purposes
		SHOOTR->Convr = 1.0; // was 0.1 or 0.38 should be 1.0 as dry matter is used in all cases
		SHOOTR->AWUPS = 0.0;  //initialize AWUPS, AWUPS_old and LeafWP in 2DSOIL Yang 8/15/06
		SHOOTR->LeafWP = -0.5;
		SHOOTR->PCRS = 0.0;
		SHOOTR->ET_demand = 0.0;
		SHOOTR->HourlyCarboUsed = 0;  //it is also zero'd upon initialization in 2dsoil
		Period = time_public->TimeStep / 60.0 / 24.0; // period should be in days, input in minutes
		PopSlab = SHOOTR->PopRow / 100.0 * SHOOTR->EOMult;
		SHOOTR->isEmergeStart = 0;
		/*These two lines show the relationships among some of the space variables.
			  --PopSlab=SHOOTR->PopRow/100*SHOOTR->RowSp*SHOOTR->EOMult;
			  --PlantDensity=SHOOTR->PopRow*100.0/SHOOTR->RowSp;
		*/
		
		//Z A new plant model object is created and initialized (calls initialize function) here
		//  plant model IO is in the controller
		//  plant model (the main model) is in RyePlant and its derived classes
		// 
		//Z note that the "varFile, GraphicFile, LeafFile" in CROP is a c++ string
		//  .c_str() converts c++ string to a pointer, pointing to a C char array

		pSC = new RyeController(varFile.c_str(), GraphicFile.c_str(), LeafFile.c_str(), initInfo);

		//initialize nitrogen uptake with what is already in the plant
		//SK 8/20/10: this is curious but OK
		//Z note that all plant N is in "mg/plant"
		// NitrogenUptake = pSC->get_plant()->get_plantTotalNitrogen() * PopSlab; 		
		NitrogenUptake = 0.0;
		NitrogenUptakeOld = 0.0;

		SHOOTR->NDemandError = 0;
		SHOOTR->CumulativeNDemandError = 0;
		SHOOTR->isGerminatedStart = 0;
		SHOOTR->isGerminatedEnd = 0;
		time_public->RunFlag = 1;
		time_public->tNext[ModNum - 1] = pSC->get_sowingDay();
	}
	//end initialization


	// *************************** model execution **************************************************
	//SK Running the crop module step by step
	//Z  this if-statement follows the 2DSOIL steps to cumulate root water and N uptake
	//   the crop model is not necessarily run
	if (module_public->NShoot > 0)
	{
		//Z g water per slab taken up in an hour
		//  in 2dsoil, it should actually be cm^3 water per slab, but it is fine since 1g water = 1 cm^3 water
		WaterUptake = WaterUptake + SHOOTR->AWUPS * time_public->Step;		

		//Z Cumulative N (mass, g (plant slab)-1) in this time step
		//  SIncrsink is in "ug" for uptake
		//  we uniformly use "mg" for N mass in plant model, therefore, we have "/0.001", now the N uptake is in "mg"
		//EXCEPTION: for photosynthesis, leaf N content is "g/m^2"
		NitrogenUptake = NitrogenUptake + SHOOTR->SIncrSink * 0.001;    
	}
	// Note that SIncrSink has been multiplied by time step in the solute uptake routing
	// the 1000 scales from ug to mg.

	//Z run the cereal rye model at this point 
	//  use sowingday to be a "NShoot" condition because if there exist a plant (even a seed), we need to run the simulations
	//  add "pSC != NULL" to make sure there is something to run
	if (fabs(time_public->Time - time_public->tNext[ModNum - 1]) < fabs(0.001 * time_public->Step) && pSC != NULL)
	{
		//If the sowing date has come and there is not plant, let the program know so other calculations are not done
		if ((module_public->NShoot == 0) && (fabs(time_public->Time - pSC->get_sowingDay())) < 0.001)
		{
			module_public->NShoot = 1;
		}

		//Z because plant (seed) density may change during germination, emergence and growth (die by any reason)
		//  use "ActualPopSlab" as the true plant density
		//  "plantLivingFrac" rnages 0 to 1, with default value 1
		// when the routine is first entered, plantlivingfraction is 1 and 
		// the defaults are at the max. maybe should represent 0 germination at first?
		double PlantLivingFraction = pSC->get_plant()->get_plantLivingFrac();
		double ActualPopSlab = __min(1.0, PopSlab * PlantLivingFraction);
		double ActualPlantDensity = __max(1.0, initInfo.plantDensity * PlantLivingFraction);
		SHOOTR->PlantLivingFraction = PlantLivingFraction;

		//Z saturated vapor pressurea at airT kPa
		double Es;
		// calculate error for demand and actual uptake, if negative, demand is greater then uptake
		CurrentNUptakeError = NitrogenUptake / ActualPopSlab - pSC->get_plant()->get_cumulativeNitrogenDemand();
		CumulativeNUptakeError += CurrentNUptakeError;

		TWeather wthr;
		{
			wthr.HourlyOutput = time_public->HourlyOutput;
			wthr.DailyOutput = time_public->DailyOutput;
			wthr.jday = Weather->JDAY;
			wthr.time = time_public->Time - Weather->JDAY;

			// ***** Z add year computation to this block **********
			//  winter rye simulation from fall to the next year, 
			//  so "year" (such as 1990, 1991, ..., 2020, 2021, ....) is NOT a constant and need to be updated
			int mm, dd, yy;
			RyeTimer.caldat(wthr.jday, mm, dd, yy);
			//Z even initInfo should be only initialized once, the year number should be updated
			wthr.year = yy;
			initInfo.year = yy;
			int YearDayZero = RyeTimer.julday(1, 1, yy);
			wthr.doy = wthr.jday - YearDayZero + 1;
			// *****************************************************


	//?????? need to choose one of these two, there is one CO2 column in the weater file
			wthr.CO2 = Weather->CO2;
			if (initInfo.CO2 > 0)
			{
				wthr.CO2 = initInfo.CO2;         //Can set CO2 in initials for specific simulations where CO2 is constant
			}
	//??????

			wthr.airT = Weather->TAIR[time_public->iTime - 1];
			wthr.PFD = Weather->par[time_public->iTime - 1] * 4.6;					//Z conversion from PAR in (W m-2) to (umol s-1 m-2)
			wthr.solRad = Weather->WATTSM[time_public->iTime - 1];					//Z hourly Total Radiation incident at soil surface (W m-2)
			Es = (0.611 * exp(17.502 * wthr.airT / (240.97 + wthr.airT)));			//Z saturated vapor pressure (kPa)
			wthr.RH = (1.0 - (Weather->VPD[time_public->iTime - 1] / Es)) * 100.0;	//Z relative humidity in percent
			wthr.rain = Weather->RINT[time_public->iTime - 1];

			//Z wind speed, converts from km hr-1 to m s-1
			wthr.wind = Weather->WIND / 3.6; 

			//Z daylength, hour
			wthr.dayLength = Weather->daylng;

			//Z leaf water potential Yang 8/15/06 MPa, SHOOTR->LeafWP = -0.5 bar
			//  since LeafWP in 2dsoil is in bar but in maizesim is in MPa, so, have to
			//  divide it by 10 to convert it into MPa before passing the value to Maizesim 1 bar = 100 kPa
			wthr.LeafWP = SHOOTR->LeafWP / 10.0;

			//Z basic and max allowed biomass rates for root growth, hourly rate
			//  (g plant^-1 hour^-1) = (g per slab per day^-1 ) / (plant per slab) / (24 hour day^-1)
			wthr.pcrl = SHOOTR->PCRL / ActualPopSlab / 24.0;
			wthr.pcrq = SHOOTR->PCRQ / ActualPopSlab / 24.0;
			//cout << wthr.time << "  PCRL:" << wthr.pcrl << "   PCRQ:" << wthr.pcrq << endl;

			// pass actual carbohydrate amount used in root (@DSOIL) back to the plant
							//ToDo - make pcrs a new variable (ActualRootCarboUsed) and make it a member of plant.
							//dt here I changed this temporarily for debugging
							//don't need to divide by 24 since the value has been integrated over an hour
							//dividing it by PopSlab converts it to g/day/plant;
							//ToDo: need to document this better, what is pcrs being used for.
			//Z the 1 hour time integration is done in soil module
			wthr.pcrs = SHOOTR->HourlyCarboUsed / ActualPopSlab;
			//cout << "PCRL:" << wthr.pcrl << "  PCRS:" << wthr.pcrs << endl;
			SHOOTR->HourlyCarboUsed = 0.0;

			//YY If time is at 5 am, then pass the leaf water potential (the predawn leaf water potential) from SHOOTR to the wthr object. 
			if (abs(wthr.time - 0.2083) < 0.0001) 
			{
				//YY. LeafWP is in bar. Since the LWPeffect in leaf.cpp uses leaf water potential in bar,
				//    so here PredawnLWP is in bar, instead of being scaled to MPa.
				wthr.PredawnLWP = SHOOTR->LeafWP; 
			}

			//Pass through nitrogen uptake (total mg per slab in the one hour) from 2DSOIL. 
			wthr.TotalRootWeight = SHOOTR->TotalRootWeight / ActualPopSlab;
			wthr.MaxRootDepth = SHOOTR->MaxRootDepth;
			// Available water is cm per profile - should be divided by PopSlab
			wthr.ThetaAvail = node_public->ThetaAvail / ActualPopSlab;

/* ??? Z comment: This, want to make hourly N input into the model system 		
       Z comment: More importantly, when rye emerge graduately, at the time of emergence, maybe only 5 plants out of 500 seeds emerge in the first hour
	              Therefore, "/ PopSlab" could be dangerous, where actual "PopSlab" in that hour can be as small as 0.001

			if (NitrogenUptake > 0)
			{
				double nuptake = NitrogenUptake / PopSlab; //uptake is now grams per plant
				double leafloss = pSC->getPlant()->get_droppedLfArea();
				double nloss = leafloss * SLNmin;

				//SK 8/20/10: Here seems to be the only place where totalN of the plant is set. NitrogenUptake is initiated from get_N at the begining of the timestep so OK. 


				pSC->getPlant()->set_N(NitrogenUptake / PopSlab);  // Units are converted from g slab-1 to g plant -1 YY
					// need to look at loss of N in the dropped leaf (plant N goes negative?)				                                                         
					//pSC->getPlant()->set_N(NitrogenUptake/PopSlab-pSC->getPlant()->get_droppedLfArea()*SLNmin);
			}
*/

			//SLNmin base Specific leaf nitrogen content; for now assume it's 0.5 YY
// debuging dt 1-30-2012

			if (SHOOTR->LAI == 0)
			{
				wthr.ET_supply = 0;
			}
			else
			{
				//Z WaterUptake is hourly summed in 2DSOIL and zeroed (at the end of this code)
				//  g water plant^-1 hour^-1
				wthr.ET_supply = WaterUptake / ActualPopSlab;
				//Z ET_diff to see possible water deficit
				ET_diff = wthr.ET_supply * 24.0 - SHOOTR->ET_demand;
			}

			//Z add daily max and min temperature for Vernalisation
			wthr.airTDmax = Weather->TAIR[0];
			wthr.airTDmin = Weather->TAIR[0];
			for (int ii = 1; ii < 24; ii++)
			{
				wthr.airTDmax = __max(wthr.airTDmax, Weather->TAIR[ii]);
				wthr.airTDmin = __min(wthr.airTDmin, Weather->TAIR[ii]);
			}
		}

		//Z loop to compute soil conditions, such as soil temperature and water condition
		int GridNodeNum = grid_public->NumNP;
		int count = 0;
		//Z cumulators for soil tmpr and soil degree of saturation
		double soilTmpr = 0.0;
		double soilWfps = 0.0;
		//Z elevation of soil surface
		double maxY = 0.0;
		//Z soil matric potential
		double soilMP[400];
		
		//Z Lowerboundary, average soil Tmpr and degree of saturation till that depth
		//  this loop means we only care about the top 5 cm soil layer
		double LowerBoundary = 5;
		for (int ii = 0; ii < GridNodeNum; ii++)
		{
			if (grid_public->y[ii] > maxY) maxY = grid_public->y[ii];
		}
		LowerBoundary = maxY - LowerBoundary;

		// average temperature in layer between surface and lower boundary
		// average degree of saturation (water content/porosity) in layer between surface and lower boundary
		for (int ii = 0; ii < GridNodeNum; ii++)
		{
			if (grid_public->y[ii] >= LowerBoundary)
			{
				soilTmpr = soilTmpr + node_public->Tmpr[ii];
				soilWfps = soilWfps + __min(1.0, (node_public->ThNew[ii] / node_public->DegreeSaturation[ii]));
				soilMP[ii] = node_public->hNew[ii];
				count++;
			}
		}

		wthr.soilT = soilTmpr / count;
		wthr.wfps = soilWfps / count;
		std::qsort((void*)soilMP, count, sizeof(double), compare);
		wthr.SoilMP_med = soilMP[count / 2];


		// The model code to simulate growth ect begins here when the plant object is called :
		//TODO add some error catching code here
		int ier = pSC->get_errStatus();
		if (ier == 0)
		{
			ier = pSC->run(wthr); //Pass weather  into the "run" function of the controller pSC YY
		}

		//Z after the plant model run once, 
		//  update living plant number
		//  i.e. plantLivingFraction
		PlantLivingFraction = pSC->get_plant()->get_plantLivingFrac();
		ActualPopSlab = __min(1.0, PopSlab * PlantLivingFraction);
		ActualPlantDensity = __max(1.0, initInfo.plantDensity * PlantLivingFraction);
		SHOOTR->PlantLivingFraction = PlantLivingFraction;

		//Z write some output to the screen 
		//  prepare some variables that return to 2DSOIL via the defined structures
		//  Assumes that germination takes place about halfway through the sowing date
		
		if (!pSC->get_plant()->get_develop()->is_germinationEnd() && (pSC->get_plant()->get_develop()->is_germinationStart()))
		{
			if (wthr.time >= 0.49 && wthr.time <= 0.51)
			{
				cout << "Germinating:" << wthr.jday << endl;
			}
		}
		// The remaining groups of code handle carbon and nitrogen exchange between 2dsoil and maizsim
		//Z germination finished but not emerge yet
		//  begin root growth at germination
		//cout << " actual pop slab: " << ActualPopSlab << endl;
		if (pSC->get_plant()->get_develop()->is_germinationStart())
		{
			//Z root carbon pools are stored in plant
			//  root mass that send to 2DSOIL is the real root mass
			// need to know germination boundaries in 2dsoil so we can
			// increase root mass as more seeds germinate
			if (!SHOOTR->isGerminatedStart)
			{
				SHOOTR->isGerminatedStart = 1;
				// get initial root mass to distribute over initial nodes
				// 2DSOIL computes based on slab while corp model computes based on plant
				SHOOTR->InitialRootCarbo = pSC->get_plant()->get_rootMass() * ActualPopSlab; 
				pSC->get_plant()->add_BiomassRootAllocationDuringGermination(SHOOTR->InitialRootCarbo);
			}
			if (!SHOOTR->isEmergeEnd)
			{
				if (pSC->get_plant()->get_develop()->is_emergenceEnd())
				{
					SHOOTR->isEmergeEnd = 1;
					cout << "emergence ended" << endl;
				}
				
				SHOOTR->InitialRootCarbo = pSC->get_plant()->get_rootMass() * (ActualPopSlab-ActualPopSlab_prev);
				pSC->get_plant()->add_BiomassRootAllocationDuringGermination(SHOOTR->InitialRootCarbo);
				double rm = pSC->get_plant()->get_rootMass();
				//cout << SHOOTR->InitialRootCarbo <<" " <<(ActualPopSlab-ActualPopSlab_prev) << " " << pSC->get_plant()->get_BiomassRootAllocationDuringGermination() <<
				//	" " << SHOOTR->TotalRootWeight << endl;
				ActualPopSlab_prev = ActualPopSlab;
			}
		}
		// pass appropriate data to 2DSOIL file structures 
		//Z this can be done by assigning values to the predefined structures
		if (pSC->get_plant()->get_develop()->is_emergenceStart())
		{
			if (!SHOOTR->isEmergeStart)
			{
				cout << "Emerging:" << wthr.jday << endl;
				SHOOTR->isEmergeStart = 1;
			}
			// ActualCarboIncrement is calculated from "assimilate", which in turn is calculated from photosynthsis_net in plant;
			// the unit of assimilate then is in g/plant/hour, thus, at this point, pcrl has unit g/plant/hour
			// Multiply by 24 to get g plant-1 day-1; multiply by popslab to get g Carbo slab-1 day-1
			
			//Z DT add a pool to hold carbon left over for previous root growth, 
			//  but Z thinks it should not be here but in the plant class during mass partitioning
			//  also Z will define biomasspools to hold roots

			SHOOTR->PCRL = (pSC->get_plant()->get_actualRootBiomassAssignment_PCRL()) * 24.0 * ActualPopSlab;

			if (!pSC->get_plant()->get_develop()->is_grainFillBegan())
			{
				SHOOTR->PCRQ = SHOOTR->PCRL +
					pSC->get_plant()->get_actualShootBiomassAssignment() * 24.0 * ActualPopSlab;
			}
			else
			{
				SHOOTR->PCRQ = SHOOTR->PCRL +
					0.75 * pSC->get_plant()->get_actualShootBiomassAssignment() * 24.0 * ActualPopSlab;
			}
			//Z after determine 
			wthr.pcrl = SHOOTR->PCRL / ActualPopSlab / 24.0;
			wthr.pcrq = SHOOTR->PCRQ / ActualPopSlab / 24.0;

			//Z plant density = plant number / m^2 but leaf area is in cm^2, thus need (100.0 * 100.0)
			SHOOTR->LCAI = pSC->get_plant()->get_greenLeafArea() * ActualPlantDensity / (100.0 * 100.0);  
			SHOOTR->LAI = (float)SHOOTR->LCAI;
			SHOOTR->Cover = 1.0 - exp(-0.79 * SHOOTR->LCAI);
			SHOOTR->Shade = (float)SHOOTR->Cover * SHOOTR->RowSp * SHOOTR->EOMult;
			SHOOTR->Height = __min(SHOOTR->Shade, SHOOTR->RowSp);
			SHOOTR->ET_demand = (pSC->get_plant()->get_ET() * 24.0); //Z this ET should be based on plant class, unit conversion should be in the plant class
			
			//YY: calculate total shoot mass per meter aquared and its increase
			//    The biomass returned by getPlant()->get_shootMass() is the weight of each single plant (g/plant), 
			//    use "pSC->getIniInfo().plantDensity" to convert it into (g/m-2)
			shoot_weightPerM2 = pSC->get_plant()->get_shootMass() * ActualPlantDensity;
			massIncrease = __max((shoot_weightPerM2 - old_shoot_weightPerM2), 0.0); //Calculated increase in above-ground biomass per m2 YY

			// optimal N ratio according to N Dilution ratio
			//Z if shoot weight < 100.0, N% = max allowed N%, i.e., no N stress
			//  if shoot weight > 100.0, will be N rate issue.
			double NitrogenRatio;     
			if (shoot_weightPerM2 < 10.0)
			{
				NitrogenRatio = N_min / 100.0; //when shoot weight is lower than 100 g m-2, then nitrogen concentration is assumed to by .0410 g/g
			}
			else
			{
				//concentration as a function of aboveground biomass (Greenwood et al., 1990; Lindquist et al., 2007) YY
				NitrogenRatio = N_min / 100.0 * pow(shoot_weightPerM2, (1.0 - N_shape)); //sqrt(shoot_weightPerM2); //Calcualte above ground potential N concentration 
			}
			pSC->get_plant()->set_NitrogenRatio(NitrogenRatio / 10.0);


			//YY U_N maximum observed N uptake rate (g N m-2 ground d-1) (Lindquist et al, 2007) 
			//	 The unit of U_N is g N m-2 ground d-1
			//	 d: shape coefficient in the logistic function to simulate cumulative N uptake (Equation 9 in Lindquist et al. 2007), double d=075
			U_N = 0.359 * d / 4.0; 
			
			//YY U_M maximum uptake rate as limited by maximum N fraction per unit (Equation 2 in Lindquist et al., 2007), g N m-2 ground d-1
			//   the unit of massIncrease is g m-2/step (here one hour), thus, need 24.0
			//   Here in the simulation, the default length of one step is an hour; so, we have to scale it up to one day by multiplying it by 24
			//	 q_n the maximum ratio of daily N uptake to measured daily growth rate (g N g-1) (Lindquist et al., 2007) double q_n = 0.032;
			U_M = q_n * massIncrease * 24.0; 
			
			//YY unit of U_P is also g N m-2 ground d-1; 
			//   however, the unit of massIncrease is g/step. Here in the simulation, 
			//   the default length of one step is an hour; so, we have to scale it up to one day by multiplying it by 24
			//
			if (shoot_weightPerM2 < 10) //if shoot weight<100 (g m-2) then U_P is calculated this way 
			{
				U_P = (N_min / 100.0) * massIncrease; // U_P potential rate of N accumulation (g N m-2 ground d-1) (Lindquist et al. 2007)
			}
			else //otherwise, it is calculated like this (Equation 6, Lindquist et al., 2007) YY
			{
				U_P = ((1.0 - N_shape) * N_min / 10.0) * pow(shoot_weightPerM2, -N_shape) * massIncrease;
			}
			
			//YY U_D U uptake rate (g N m-2 d-1) as limited by the difference between potential and actual amount of N in existing biomass, equation 3 in Lindquist et al. 2007)
			//   value from get_N() is in g N/plant (need unit adjustment from mg to g). 
			//   It has to be converted to g/m-2 ground that's why the actual n content is mulitpled by pSC->getIniInfo().plantDensity/(100*100)
			U_D = N_min * 10 / 100 * pow(shoot_weightPerM2, -N_shape) - pSC->get_plant()->get_plantTotalNitrogen() * (ActualPlantDensity / (100.0 * 100.0));
			

			// set up accounting of N here.
			// first hourly//Actual and needed N uptake in the last hour per plant per day
			double HourlyActualNFromSoil = (NitrogenUptake - NitrogenUptakeOld) / ActualPopSlab; //houly rate per hour
			double HourlyNitrogenDemand = __max(U_P, 0.0); // ActualPlantDensity / 24.0; //Determine the hourly nitrogen demand (equation 1 Lindquist et al. 2007) in grams plant-1
 			pSC->get_plant()->set_HourlyNitrogenSoilUptake(HourlyActualNFromSoil);
			pSC->get_plant()->set_HourlyNitrogenDemand(HourlyNitrogenDemand);

			//cout << "N_fromSoil:"<<HourlyActualNFromSoil << " Demand:" << HourlyNitrogenDemand << '\n';

			// now do cumulative amounts
			CumulativeNitrogenDemand += HourlyNitrogenDemand; // grams N plant-1 hour-1
			pSC->get_plant()->set_CumulativeNitrogenDemand(CumulativeNitrogenDemand); //g N plant hour-1  
			pSC->get_plant()->set_CumulativeNitrogenSoilUptake(NitrogenUptake / ActualPopSlab);


			SHOOTR->nitroDemand = (float)HourlyNitrogenDemand * (float)ActualPopSlab * 1e6 * 24.0; //Pass the nitrogen demand into 2dsoil YY
																		   //Units are ug slab-1
			old_shoot_weightPerM2 = shoot_weightPerM2; //Save the value of the above_ground biomass of this time-step
			NitrogenUptakeOld = NitrogenUptake; // save the cumulative N uptake from this time step;
			SHOOTR->NDemandError = (float)CurrentNUptakeError;
			SHOOTR->CumulativeNDemandError = (float)CumulativeNUptakeError;

		}
			//			if (pSC->getLastDayOfSim() <= pSC->getTime()->get_day_of_year()) 
		if (pSC->get_plant()->get_develop()->is_dead())
		{
			cout << "Completing crop simulation..." << endl;
			module_public->NShoot = 0; //tell 2dsoil that crops harvested
			time_public->tNext[ModNum - 1] = 1e12; // set the next time so the model
			if (pSC != NULL)
			{
				delete pSC;
				pSC = NULL;
			}
			time_public->RunFlag = 0;
		}
		else
		{
			time_public->tNext[ModNum - 1] = time_public->Time + Period;
			WaterUptake = 0;
			//NitrogenUptake=0;
		}
	}// end hourly calculations code
 // end if NShhot>0 section 
	return;
}