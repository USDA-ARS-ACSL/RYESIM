//
#include "stdafx.h"
#include "RyePlant.h"
#include "radtrans.h"
#include "timer.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>

#define MINUTESPERDAY 1440.0
#define DAYPERMINUTES 0.00069444444
#define CO2_MW 44.0098
#define C_MW 12.011
#define CH2O_MW 30.03
#define MAXLEAFNUM 20
#define MAX_N_PCT 4.0
#define MIN_N_PCT 0.5
#define MAX_N_PCT_ROOT 1.5
#define cPFD 4.60

using namespace std;

RyePlant::RyePlant(const TInitInfo& info, TGasExSpeciesParam& photoparam)
{
	//Z set pointer to NULL before allocate memory
	mainstem = NULL;
	roots = NULL;
	develop = NULL;
	sunlit = NULL;
	shaded = NULL;

	PlantLivingFraction = 1.0;

	//Z initialize leaf number, area, mass
	//  by default, this will be representative numbers, rather than "one single plant" number
	//  length in cm, area in cm^2, biomass in g, nitrogen mass in mg
	LeafNumPlt = 0.0;
	GreenLfNumPlt = 0.0;
	DropLfNumPt = 0.0;

	LeafAreaPlt = 0.0;
	GreenLfAreaPlt = 0.0;
	GreenLfLengthPlt = 0.0;
	GreenLfWidthPlt = 0.0;
	SenescentLfAreaPlt = 0.0;
	DropLfAreaPlt = 0.0;

	LeafMassPlt = 0.0;
	SheathMassPlt = 0.0;
	DropLeafMassPlt = 0.0;
	DropSheathMassPlt = 0.0;
	LeafNitrogenMassPlt = 0.0;
	SheathNitrogenMassPlt = 0.0;
	LeafNitrogenReleasePlt = 0.0;
	SheathNitrogenReleasePlt = 0.0;

	LeafNitrogenDeadMassPlt = 0.0;
	SheathNitrogenDeadMassPlt = 0.0;
	LeafNitrogenYellowMassPlt = 0.0;
	SheathNitrogenYellowMassPlt = 0.0;
	LeafNitrogenDroppedMassPlt = 0.0;
	SheathNitrogenDroppedMassPlt = 0.0;

	PtnLeafMassIncreasePlt = 0.0;
	PtnLeafNitrogenMassIncreasePlt = 0.0;

	//Z initialize tiller and internode number, mass, length
	//  length in cm, area in cm^2, biomass in g, nitrogen mass in mg
	TillerNumSingle = 0;
	TillerNumPlt = 0.0;
	IntrNodeLengthPlt = 0.0;
	IntrNodeMassPlt = 0.0;
	IntrNodeNitrogenMassPlt = 0.0;
	IntrNodeNitrogenReleasePlt = 0.0;

	IntrNodeNitrogenDeadMassPlt = 0.0;

	PtnIntrNodeMassIncreasePlt = 0.0;
	PtnIntrNodeNitrogenMassIncreasePlt = 0.0;

	//Z initialize plant based mass and N variables, % number
	LeafNitrogenContentPlt = MAX_N_PCT;

	//Z initialize plant part sizes
	SeedMass = 0.04; // g seed^-1
	SeedNitrogenMass = SeedMass * MAX_N_PCT * 10.0; // total weight of N per seed (mg N seed^-1, aka, mg N plant^-1)
	//  mass partitioning to plant organs
	shootPart = 0.60;
	rootPart = 0.40;
	shootPart_old = 0.0;
	rootPart_old = 0.0;
	leafPart = 0.90;
	internodePart = 0.10;

	RyeMass = SeedMass;
	ShootMass = SeedMass * shootPart;
	RootMass = SeedMass * rootPart;
	LeafMassPlt = leafPart * ShootMass;
	IntrNodeMassPlt = internodePart * ShootMass;
	NitrogenMassPlt = SeedNitrogenMass; // total weight of N per seed (mg N seed^-1, aka, mg N plant^-1)

	//Z initialize root mass, actual root growth and hence the N demand
	//  root biomass allocation has its own rule, 
	//  while root N is affliated to root biomass
	//  therefore, we only assume potential root N mass
	RootNitrogenMassPlt = 0.0;
	NitrogenRootGrowth_ptn = 0.0;
	
	//Z photosynthesis and gas exchange model
	sunlit_LAI = shaded_LAI = 0.0;
	sunlit_PFD = shaded_PFD = 0.0;

	sunlit_A_net = shaded_A_net = 0;
	sunlit_A_gross = shaded_A_gross = 0;
	sunlit_gs = shaded_gs = 0;
	sowingDay = 1.0;
	age = 0.0;
	initInfo = info;
	gasExparam = photoparam;

	//Z initialize the plant body
	develop = new RyeDevelopment(initInfo);
	roots = new RyeRoots();
	//Z we start the plant by initialize the first tiller, aka the mainstem
	//  the initial living fraction of the mainstem MUST be 1.0
	//Z the mainstem does not have "rank" and its order is 0, hence its cumuRank is also 0;
	mainstem = new RyeTiller(0, 0, 0, true, develop, 1.0);

	photosynthesis_net = 0.0;
	photosynthesis_gross = 0.0;
	transpiration = 0.0;
	transpirationOld = 0.0;
	assimilate = 0.0;
	VPD = 0;
	conductance = 0;
	emerge_gdd = 0;
	SunlitRatio = 0.0;
	temperature = develop->get_Tcur();

	//Z initialize the biomass partitioning factors
	//  pools and storage factors
	BiomassReserve = 0.00001;
	BiomassPool = SeedMass;
	NitrogenPool = SeedNitrogenMass;
	BiomassSupply = 0.0;
	BiomassDemand = 0.0;
	BiomassLeftover = 0.0;
	BiomassRoot = 0.00001;
	BiomassLeafGrowth_ptn = 0.0;
	BiomassIntrNodeGrowth_ptn = 0.0;
	BiomassShootGrowth_ptn = 0.0;
	BiomassPltGrowth_ptn = 0.0;
	
	//Z all the following nitrogen is based on mg plant^-1
	CumulativeNitrogenDemand = 0.0;
	HourlyNitrogenDemand = 0.0;
	CumulativeNitrogenSoilUptake = 0.0;
	HourlyNitrogenSoilUptake = 0.0;
	NitrogenShootGrowth_ptn = 0.0;
	ShootNitrogenAvailiableAllocation = 0.0;

	ActuralRootBiomassAssignment_PCRL = 0.0;
	ActuralShootBiomassAssignment = 0.0;
	C2_effect = 0.0;

	LeafBiomassAssignmentRate = 0.0;
	IntrNodeBiomassAssignmentRate = 0.0;
	LeafNitrogenAssignmentRate = 0.0;
	IntrNodeNitrogenAssignmentRate = 0.0;
	maintRespiration = 0.0;

}

//Z even plant died, this does not mean plant objective will be destructed.
//  destruction only occurs when simulation is finished
//  need to delete and release all the pointers
RyePlant::~RyePlant()
{
	if (develop != NULL) delete develop;
	develop = NULL;
	if (roots != NULL) delete roots;
	roots = NULL;
	if (mainstem != NULL) delete mainstem;
	mainstem = NULL;
}

void RyePlant::RyePlantUpdate(const TWeather& weather)
{
	//Z update the development stage
	develop->RyeDelpUpdate(weather);
	//Z this living fraction is NOT the same as the tiller, leaf, internode livingfrac
	//  this is for emergence percentage
	PlantLivingFraction = develop->get_plantLivingFraction();

	//Z update the plant morphology
	//Z first germination starts (not complete), just mark the time
	//  germination start, namely "germination initiation, germinInit", should be always earlier than anything
	if (!develop->is_germinInit())
	{
		//Z seed is not germinated yet, nothing need to be done
		temperature = develop->get_Tcur();
	}
	
	//Z germination marks 1/2 of the "active" seeds germinate, based on a germination ratio
	//  may complete after emergence
	//  so germination seems not that useful for plant.cpp, because in the module, we still tend to simulate one representative plant
	//  but maybe important for output
	if ((develop->is_germinInit()) && (!develop->is_germination()))
	{
		temperature = develop->get_Tcur();
	}

	//Z between "germination start" and "emergence start"
	//  representative plant, start to trace
	//  seed is the source of biomass and nitrogen
	if ((develop->is_germinInit()) && (!develop->is_emerge()))
	{
		//Z initialize the root
		//Z root initialization is not that critical at all
		//  the root mass will be assigned from plant class to the 2DSOIL model
		//  so the mass in root class will just be an "apparent number" that will never be used
		//Z I do not know why I have to put root class in this rysim model, can be simply removed in future
		temperature = develop->get_Tcur();
		if (!roots->get_Initialization())
		{
			roots->set_Initialization();
			roots->inputBiomass(RootMass);
		}
		//Z mainstem initialization, 
		//  apsim assumes two leaves grow on the mainstem
		//  therefore, initialize the two leaves before emergence
		if (!mainstem->get_mainstemInitilization())
		{
			mainstem->MainstemInitialize();
		}
		
		//Z update morph features, like leaf area, stem length, mass and N etc.
		//  this need to be done after we initialize 2 leaves and 2 internodes
		//Z after update the growth, also compute the biomass requirement under a very ideal condition
		//		+ potential leaf area increase * 0.0045 g/cm^2 (specific weight)
		//      + potential internode increase * 0.015 g/cm (specific weight)
		//Z after update the growth, also compute the nitrogen requirement under a very ideal condition
		//		+ potential leaf area increase * 0.0045 g/cm^2 (specific weight) * 0.035 (N mass fraction) * 1000 (convert to mg)
		//      + potential internode increase * 0.015 g/cm (specific weight) * 0.035 (N mass fraction) * 1000 (convert to mg)
		RyeRecursiveUpdateSummary();
		
		//Z potential growth can be consider as the demand
		//  biomass g
		//  nitrogen mass mg
		if (develop->is_accel())
		{
			BiomassLeafGrowth_ptn = PtnLeafMassIncreasePlt * 2.0;
			BiomassIntrNodeGrowth_ptn = PtnIntrNodeMassIncreasePlt;
			BiomassShootGrowth_ptn = BiomassLeafGrowth_ptn + BiomassIntrNodeGrowth_ptn;
			NitrogenLeafGrowth_ptn = PtnLeafNitrogenMassIncreasePlt * 2.0;
			NitrogenIntrNodeGrowth_ptn = PtnIntrNodeNitrogenMassIncreasePlt;
			NitrogenShootGrowth_ptn = NitrogenLeafGrowth_ptn + NitrogenIntrNodeGrowth_ptn;
		}
		else
		{
			BiomassLeafGrowth_ptn = PtnLeafMassIncreasePlt * 1.5;
			BiomassIntrNodeGrowth_ptn = PtnIntrNodeMassIncreasePlt;
			BiomassShootGrowth_ptn = BiomassLeafGrowth_ptn + BiomassIntrNodeGrowth_ptn;
			NitrogenLeafGrowth_ptn = PtnLeafNitrogenMassIncreasePlt * 1.5;
			NitrogenIntrNodeGrowth_ptn = PtnIntrNodeNitrogenMassIncreasePlt;
			NitrogenShootGrowth_ptn = NitrogenLeafGrowth_ptn + NitrogenIntrNodeGrowth_ptn;
		}
		RootMass += weather.pcrs;
		NitrogenRootGrowth_ptn = __max((RootMass * MAX_N_PCT_ROOT * 10.0 - RootNitrogenMassPlt), 0.0);
		
		//Z sum the mass component at this time
		//  including biomass and nitrogen mass
		//  note seed mass is changing, but the nitrogen fraction in the seed mass is a constant
		RyeMass = SeedMass + LeafMassPlt + SheathMassPlt + IntrNodeMassPlt + RootMass;
		ShootMass = LeafMassPlt + SheathMassPlt + IntrNodeMassPlt;
		NitrogenMassPlt = SeedNitrogenMass + LeafNitrogenMassPlt + SheathNitrogenMassPlt + IntrNodeNitrogenMassPlt
			+ LeafNitrogenDeadMassPlt + SheathNitrogenDeadMassPlt + IntrNodeNitrogenDeadMassPlt+ RootNitrogenMassPlt;
			
		//Z plant based leaf N content (aged leaves release N)
		//  in photosynthesis gas_exchange model, unit mg N per m^2 leaf area
		if (GreenLfAreaPlt <= 0.1) {
			// default N mass (g) per leaf area m^2, with N mass fraction 3.5% and specific leaf weight 0.000045 g/mm^2
			LeafNitrogenContentPlt = 1.5750 / 3.5 * MAX_N_PCT;
		}
		else {
			// convert leafN mass from mg to g
			// convert green leaf area to m^2
			LeafNitrogenContentPlt = (LeafNitrogenMassPlt / 1000.0) / (GreenLfAreaPlt / 10000.0);
		}
		
		//Z maintanence respiration rate
		calcMaintRespiration(weather);

		//Z biomass allocation
		//  in this case, biomass pools are all from seeds, so need to decrease SeedMass step by step
		//  then set mass to the plant organs
		calsBiomassAllocation(weather);
		
		//Z nitrogen allocation
		//  in this case, nitrogen are all from seeds, with seed N fraction (3.4%)
		//  then set nitrogen mass to the plant organs
		//Z Nitrogen release from seed should be propotional to the BiomassSupply based on (3.4%)
		calsNitrogenAllocation();

		SeedMass = __max(0.0, SeedMass - BiomassSupply);
		SeedNitrogenMass = __max(0.0, SeedNitrogenMass - RootNitrogenAvailiableAllocation - ShootNitrogenAvailiableAllocation);
		
		calsSetMass();

		//Z sum the mass component at this time
		//  including biomass and nitrogen mass
		//  note seed mass is changing, but the nitrogen fraction in the seed mass is a constant
		RyeMass = SeedMass + LeafMassPlt + SheathMassPlt + IntrNodeMassPlt + RootMass;
		ShootMass = LeafMassPlt + SheathMassPlt + IntrNodeMassPlt;
		NitrogenMassPlt = SeedNitrogenMass + LeafNitrogenMassPlt + SheathNitrogenMassPlt + IntrNodeNitrogenMassPlt
			+ LeafNitrogenDeadMassPlt + SheathNitrogenDeadMassPlt + IntrNodeNitrogenDeadMassPlt+ RootNitrogenMassPlt;
	}
	//Z after emergence
	//Z need to obtain hourly root nitrogen uptake from 2DSOIL as nitrogen input
	if ((develop->is_emerge()) && (!develop->is_dead()))
	{
		//Z from this point, plant will grow in a "normal way" or say, on its "path"
		//  update the temperature and then compute the R:FR ratio, which will furtherly used as a limiter for tiller number
		temperature = develop->get_Tcur();
		calcRed_FRedRatio(weather);

		//Z update morph features, like leaf area, stem length, mass and N etc.
		//Z after update the growth, also compute the biomass requirement under a very ideal condition
		//		+ potential leaf area increase * 0.0045 g/cm^2 (specific weight)
		//      + potential internode increase * 0.015 g/cm (specific weight)
		//Z after update the growth, also compute the nitrogen requirement under a very ideal condition
		//		+ potential leaf area increase * 0.0045 g/cm^2 (specific weight) * 0.035 (N mass fraction) * 1000 (convert to mg)
		//      + potential internode increase * 0.015 g/cm (specific weight) * 0.035 (N mass fraction) * 1000 (convert to mg)
		RyeRecursiveUpdateSummary();

		//Z potential growth can be consider as the demand
		//  biomass g
		//  nitrogen mass mg
		if (develop->is_accel())
		{
			BiomassLeafGrowth_ptn = PtnLeafMassIncreasePlt * 2.0;
			BiomassIntrNodeGrowth_ptn = PtnIntrNodeMassIncreasePlt;
			BiomassShootGrowth_ptn = BiomassLeafGrowth_ptn + BiomassIntrNodeGrowth_ptn;
			NitrogenLeafGrowth_ptn = PtnLeafNitrogenMassIncreasePlt * 2.0;
			NitrogenIntrNodeGrowth_ptn = PtnIntrNodeNitrogenMassIncreasePlt;
			NitrogenShootGrowth_ptn = NitrogenLeafGrowth_ptn + NitrogenIntrNodeGrowth_ptn;
		}
		else
		{
			BiomassLeafGrowth_ptn = PtnLeafMassIncreasePlt * 1.5;
			BiomassIntrNodeGrowth_ptn = PtnIntrNodeMassIncreasePlt;
			BiomassShootGrowth_ptn = BiomassLeafGrowth_ptn + BiomassIntrNodeGrowth_ptn;
			NitrogenLeafGrowth_ptn = PtnLeafNitrogenMassIncreasePlt * 1.5;
			NitrogenIntrNodeGrowth_ptn = PtnIntrNodeNitrogenMassIncreasePlt;
			NitrogenShootGrowth_ptn = NitrogenLeafGrowth_ptn + NitrogenIntrNodeGrowth_ptn;
		}
		RootMass += weather.pcrs;
		NitrogenRootGrowth_ptn = __max((RootMass * MAX_N_PCT_ROOT * 10.0 - RootNitrogenMassPlt), 0.0);
		
		//Z plant mass (above ground part)
		//  including biomass and nitrogen mass
		RyeMass = SeedMass + LeafMassPlt + SheathMassPlt + IntrNodeMassPlt + RootMass;
		ShootMass = LeafMassPlt + SheathMassPlt + IntrNodeMassPlt;
		NitrogenMassPlt = SeedNitrogenMass + LeafNitrogenMassPlt + SheathNitrogenMassPlt + IntrNodeNitrogenMassPlt
			+ LeafNitrogenDeadMassPlt + SheathNitrogenDeadMassPlt + IntrNodeNitrogenDeadMassPlt+ RootNitrogenMassPlt;

		//Z plant based leaf N content (aged leaves release N)
		if (GreenLfAreaPlt <= 1.0)
		{
			// default N mass (g) per leaf area m^2, with N mass fraction 3.5% and specific leaf weight 0.000045 g/mm^2
			LeafNitrogenContentPlt = 1.5750 / 3.5 * MAX_N_PCT;
		}
		else
		{
			// convert leafN mass from mg to g
			// convert green leaf area cm^2 to m^2
			LeafNitrogenContentPlt = (LeafNitrogenMassPlt / 1000.0) / (GreenLfAreaPlt / 10000.0);
		}

		//Z judge if plant is nearly dead
		//  currently, we use if the elongation start as the temporary criterion
		if ((GreenLfAreaPlt <= 0.00005 * LeafAreaPlt) && develop->is_startEnlongation())
		{
			develop->set_maturity(true, weather.daytime);
			develop->set_death(true, weather.daytime);
		}

		//Z call the photosynthesis method
		//  the criterion for maize is "GreenLfAreaPlt>10 cm^2", but that could be too large for rye
		if (GreenLfAreaPlt > 1.0)
		{
			calcGasExchange(weather, gasExparam);
		}
		
		//Z call maintainance respiration
		calcMaintRespiration(weather);

		//Z biomass partitioning between roots and shoots
		//  then partitioning shoots part to leaf (with sheath) and internode
		calsBiomassAllocation(weather);

		//Z nitrogen allocation
		//  in this case, nitrogen are from plant nitrogen release and root N uptake
		//  then set nitrogen mass to the plant organs
		//Z at this time, biomass and nitrogen may not be allocated at a fixed ratio, 
		//  and hence, plant organs will have 
		calsNitrogenAllocation();

		//Z assign photosynthesis to the biomass pools
		//  this comes after the biomass allocation, that means, the photosynthesis biomass will be used in the next time step
		if (abs(weather.time) < 0.0001)
		{
			//Z by the end of a calender day, zero short term biomass pool
			//  add short term biomass pool to the long term biomass reservior
			BiomassReserve += __max(0.0, BiomassPool);
			BiomassPool = 0.0; //reset shorterm C_pool to zero at midnight, needs to be more mechanistic
		}
		else
		{
			//Z if not at the end of a calender day, take photosynthesis biomass output
			//  should be assigned to the fast/short term biomass pool
			BiomassPool += __max(photosynthesis_gross,0.0);
		}
		calsSetMass();

		//Z sum the mass component at this time
		//  including biomass and nitrogen mass
		//  note seed mass is changing, but the nitrogen fraction in the seed mass is a constant
		RyeMass = SeedMass + LeafMassPlt + SheathMassPlt + IntrNodeMassPlt + RootMass;
		ShootMass = LeafMassPlt + SheathMassPlt + IntrNodeMassPlt;
		NitrogenMassPlt = SeedNitrogenMass + LeafNitrogenMassPlt + SheathNitrogenMassPlt + IntrNodeNitrogenMassPlt
			+ LeafNitrogenDeadMassPlt + SheathNitrogenDeadMassPlt + IntrNodeNitrogenDeadMassPlt+ RootNitrogenMassPlt;

	}
}

//Z make recursive summary of plant morph features
//  cumulate leaf and internode features tiller by tiller
void RyePlant::RyeRecursiveUpdateSummary()
{
	//Z need to reset the plant parameters before sum
	//  recall that tiller is a box that hold leaves and internodes
	// 
	//Z reset leaf number, area, mass
	LeafNumPlt = 0.0;
	GreenLfNumPlt = 0.0;
	DropLfNumPt = 0.0;

	LeafAreaPlt = 0.0;
	GreenLfAreaPlt = 0.0;
	GreenLfLengthPlt = 0.0;
	GreenLfWidthPlt = 0.0;
	SenescentLfAreaPlt = 0.0;
	DropLfAreaPlt = 0.0;

	LeafMassPlt = 0.0;
	SheathMassPlt = 0.0;
	DropLeafMassPlt = 0.0;
	DropSheathMassPlt = 0.0;
	LeafNitrogenMassPlt = 0.0;
	SheathNitrogenMassPlt = 0.0;
	LeafNitrogenReleasePlt = 0.0;
	SheathNitrogenReleasePlt = 0.0;

	LeafNitrogenDeadMassPlt = 0.0;
	SheathNitrogenDeadMassPlt = 0.0;
	LeafNitrogenYellowMassPlt = 0.0;
	SheathNitrogenYellowMassPlt = 0.0;
	LeafNitrogenDroppedMassPlt = 0.0;
	SheathNitrogenDroppedMassPlt = 0.0;

	PtnLeafMassIncreasePlt = 0.0;
	PtnLeafNitrogenMassIncreasePlt = 0.0;
	
	//Z reset tiller and internode number, mass, length
	TillerNumSingle = 0;
	TillerNumPlt = 0.0;
	IntrNodeLengthPlt = 0.0;
	
	IntrNodeMassPlt = 0.0;
	IntrNodeNitrogenMassPlt = 0.0;
	IntrNodeNitrogenReleasePlt = 0.0;

	IntrNodeNitrogenDeadMassPlt = 0.0;

	PtnIntrNodeMassIncreasePlt = 0.0;
	PtnIntrNodeNitrogenMassIncreasePlt = 0.0;

	//Z remove the geometical based potential growth quantity
	//PtnLeafAreaIncreasePlt = 0.0;
	//PtnIntrNodeLengthIncreasePlt = 0.0;

	//Z recursive cumulate leaf and internode from the mainstem
	this->PltCallTillerRecursiveUpdate(this->mainstem);

	//Z now we take average of green leaf length and width
	//  need to note that green leaf length and width are first cumulative values and then averaged values
	//  in tiller, it is only cumulative values, and finally, take average here
	//  should be done after sum all the plant mass
	GreenLfLengthPlt = GreenLfLengthPlt / GreenLfNumPlt;
	GreenLfWidthPlt = GreenLfWidthPlt / GreenLfNumPlt;

	//Z set a tiller number to development, once tiller number > 4, tiller will shade each other and reduce existing tiller living fraction
	develop->set_plantTillerNum(TillerNumSingle);
}

//Z Recursive UPDATE and SUMMARY of leaf mass and internode mass
//     recursive is done for all the tillers
//  
void RyePlant::PltCallTillerRecursiveUpdate(RyeTiller* rt)
{
	//Z call the update eailier since the parent tiller will affect the child tiller
	//  "rT" in the argument means ryeTiller

	//Z morph updates of leaf and internode
	rt->RyeTillerSingleMorph();
	//Z existing tiller organs growth (leaf and internode)
	rt->RyeTillerSingleUpdate();
	//Z sum the tiller features
	rt->RyeTillerSingleSummary();
	//Z sum the tiller features to the plant for integrated evaluation of leaf area, or potential growth quantities
	//  thus, the data flow is tiller scale data -> plant scale data
	//  need to zero all the plant scale variables before this process

	//Z leaf number, area, mass
	LeafNumPlt += rt->get_tillerLeafNum();
	GreenLfNumPlt += rt->get_tillerGreenLeafNum();
	DropLfNumPt += rt->get_tillerDropLeafNum();

	LeafAreaPlt += rt->get_tillerLeafArea();
	GreenLfAreaPlt += rt->get_tillerGreenLeafArea();
	GreenLfLengthPlt += rt->get_tillerGreenLeafLength();
	GreenLfWidthPlt += rt->get_tillerGreenLeafWidth();
	SenescentLfAreaPlt += rt->get_tillerSeneLeafArea();
	DropLfAreaPlt += rt->get_tillerDropLeafArea();
	
	LeafMassPlt += rt->get_tillerLeafMass();
	SheathMassPlt += rt->get_tillerSheathMass();
	DropLeafMassPlt += rt->get_tillerDropLeafMass();
	DropSheathMassPlt += rt->get_tillerDropSheathMass();
	LeafNitrogenMassPlt += rt->get_tillerLeafNitrogenMass();
	SheathNitrogenMassPlt += rt->get_tillerSheathNitrogenMass();
	LeafNitrogenReleasePlt += rt->get_tillerLeafNitrogenRelease();
	SheathNitrogenReleasePlt += rt->get_tillerSheathNitrogenRelease();

	LeafNitrogenDeadMassPlt += rt->get_tillerLeafNitrogenDeadMass();
	SheathNitrogenDeadMassPlt += rt->get_tillerSheathNitrogenDeadMass();
	LeafNitrogenYellowMassPlt += rt->get_tillerLeafNitrogenYellowMass();
	SheathNitrogenYellowMassPlt += rt->get_tillerSheathNitrogenYellowMass();
	LeafNitrogenDroppedMassPlt += rt->get_tillerLeafNitrogenDroppedMass();
	SheathNitrogenDroppedMassPlt += rt->get_tillerSheathNitrogenDroppedMass();

	PtnLeafMassIncreasePlt += rt->get_tillerLeafMassIncrease_ptn();
	PtnLeafNitrogenMassIncreasePlt += rt->get_tillerLeafNitrogenMassIncrease_ptn();

	//Z tiller and internode number, mass, length
	//  use tiller "livingFrac" as the effective tiller number
	TillerNumSingle += rt->is_tillerLiving();
	TillerNumPlt += rt->get_tillerLivingFrac();
	IntrNodeLengthPlt += rt->get_IntrNodeLength();
	
	IntrNodeMassPlt += rt->get_IntrNodeMass();
	IntrNodeNitrogenMassPlt += rt->get_IntrNodeNitrogenMass();
	IntrNodeNitrogenReleasePlt += rt->get_IntrNodeNitrogenRelease();

	IntrNodeNitrogenDeadMassPlt += rt->get_IntrNodeNitrogenDeadMass();

	PtnIntrNodeMassIncreasePlt += rt->get_tillerIntrMassIncrease_ptn();
	PtnIntrNodeNitrogenMassIncreasePlt += rt->get_tillerIntrNitrogenMassIncrease_ptn();

	//Z the geometrical based potential growth should be replaced by mass based potential growth quantity
	//PtnLeafAreaIncreasePlt += rt->get_tillerLeafAreaIncrease_ptn();
	//PtnIntrNodeLengthIncreasePlt += rt->get_IntrNodeLengthIncrease_ptn();

	//Z recursive process, starting with the mainstem
	for (int ii = 0; ii < MAXLEAFNUM; ii++)
	{
		//Z first name this subtiller "subrt"
		RyeTiller* subrt=rt->get_subtiller(ii);
		if (subrt != NULL)
		{
			PltCallTillerRecursiveUpdate(subrt);
		}
	}
}

//Z Coupled Photosynthesis and Stomatal Conductance Model
//  Farquhar-von Caemmerer-Berry (FvCB), Ball-Woodrow-Berry (BWB), leaf-level energy balance model 
void RyePlant::calcGasExchange(const TWeather& weather, const TGasExSpeciesParam& photoparam)
{
	const double tau = 0.50;					//Z atmospheric transmittance, to be implemented as a variable => done
	const double LAF = 1.035;					//Z leaf angle factor for rye grass leaf, Campbell and Norman (1998), Table 15.1, for rye (NOT ryegrass)
	const double leafwidth = GreenLfWidthPlt;	//Z greenleaf width compute during plant growth every hour
	const double atmPressure = 101.3;			//Z kPa, to be predicted using altitude

	//Z GreenLfAreaPlt increases and decreases as plant growth and getting old
	//  LeafAreaPlt increases as plant growth but will not decrease, even leaf drops, leaf area counts
	double activeLeafRatio = this->GreenLfAreaPlt / this->LeafAreaPlt;

	//Z not all seeds germinate, not all germinated seed emerge, ......
	//  any accident can happen, so put "initInfo.plantDensity * PlantLivingFraction" as the real plant density
	//  10000.0 cm^2=1 m^2, which convert the plant leaf area unit from cm^2 to m^2
	double LAI = this->GreenLfAreaPlt * __max(1.0, initInfo.plantDensity * PlantLivingFraction) / 10000.0; 
	
	if ((weather.doy <= 30) || (weather.doy >= 200))
	{
		BiomassLeftover = 0.0;
	}
	else 
	{
		BiomassLeftover = (BiomassPool + BiomassReserve) / LAI;
	}

	CGasExchange* sunlit = new CGasExchange(photoparam);
	CGasExchange* shaded = new CGasExchange(photoparam);

	CSolar* sun = new CSolar();
	CRadTrans* light = new CRadTrans();

	//Z determine DOY and set it into the sun(CSolar) and light(CRadTrans) category
	Timer timer;
	int mm, dd, yy;
	timer.caldat(weather.jday, mm, dd, yy);
	int jday = timer.julday(1, 1, yy);
	sun->SetVal(weather.jday - jday + 1, weather.time, initInfo.latitude, initInfo.longitude, initInfo.altitude, weather.solRad);
	light->SetVal(*sun, LAI, LAF);

	sunlit_PFD = light->Qsl();		//Z CRadTrans object => mean photon flux density on sunlit leaves (umol m-2 s, unit converted from PAR w/m^2)
	shaded_PFD = light->Qsh();		//Z CRadTrans object => mean photon flux density on sunshade leaves 
	sunlit_LAI = light->LAIsl();	//Z CRadTrans object => LAI portion under sun shine
	shaded_LAI = light->LAIsh();	//Z CRadTrans object => LAI portion under sun shade

	double leaf_psi = weather.LeafWP;
    //??? cauwzj test 
	//leaf_psi = -0.05;
	//this->LeafNitrogenContentPlt = 1.5;
	//Calculating transpiration and photosynthesis with stomatal controlled by leaf water potential LeafWP Y
	sunlit->SetVal(sunlit_PFD, initInfo, weather.airT, weather.CO2,
		weather.RH, weather.wind, atmPressure, leafwidth,
		leaf_psi, this->LeafNitrogenContentPlt, BiomassLeftover);
	shaded->SetVal(shaded_PFD, initInfo, weather.airT, weather.CO2,
		weather.RH, weather.wind, atmPressure, leafwidth,
		leaf_psi, this->LeafNitrogenContentPlt, BiomassLeftover);

	//Z Get immediate values from photosynthesis
	sunlit_A_net = sunlit->A_net;
	shaded_A_net = shaded->A_net;
	sunlit_A_gross = sunlit->A_gross;
	shaded_A_gross = shaded->A_gross;
	sunlit_gs = sunlit->gs;
	shaded_gs = shaded->gs;

	//plantsPerMeterSquare units are umol CO2 m-2 ground s-1 ;
	photosynthesis_gross = sunlit_A_gross * sunlit_LAI + shaded_A_gross * shaded_LAI;
	photosynthesis_net = sunlit_A_net * sunlit_LAI + shaded_A_net * shaded_LAI;
	
	//when outputting the previous step transpiration is compared to the current step's water uptake
	transpirationOld = transpiration;  
	transpiration = 0;
	if (sunlit_LAI > 0) transpiration = sunlit->ET * sunlit_LAI;
	if (shaded_LAI > 0) transpiration += shaded->ET * shaded_LAI;
	// plantsPerMeterSquare units are grams per plant per hour;
	// Transpiration g (Water) plant^-1 hour^-1
	// Units of Transpiration from sunlit->ET/shaded->ET are mmol m-2 (leaf area) s-1
	// Calculation of transpiration from ET involves the conversion to gr per plant per hour 
	transpiration = transpiration * (60.0 * initInfo.timeStep) * 18.01 / 1000.0
		/ __max(1.0, initInfo.plantDensity * PlantLivingFraction);
	
	//psi_l = (sunlit->get_psi()*sunlitLAI + shaded->get_psi()*shadedLAI)/LAI;
	this->VPD = sunlit->VPD;

	// photosynthesis_gross is umol CO2 m-2 leaf s-1
	// in the following we convert to g C plant-1 per hour
	// assimilate: grams CO2 per plant per hour
	assimilate = (photosynthesis_gross * CO2_MW / 1.0e6) * (60.0 * initInfo.timeStep) / __max(1.0, initInfo.plantDensity * PlantLivingFraction);
	// photosynthesis_gross: grams biomass per plant per hour
	photosynthesis_gross = photosynthesis_gross * CH2O_MW / 1.0e6 * (60.0 * initInfo.timeStep) / __max(1.0, initInfo.plantDensity * PlantLivingFraction);
	// photosynthesis_net: grams biomass per plant per hour
	photosynthesis_net = photosynthesis_net * CH2O_MW / 1.0e6 * (60.0 * initInfo.timeStep) / __max(1.0, initInfo.plantDensity * PlantLivingFraction);

	if (LAI != 0)
	{
		//Z leaf temperature
		temperature = (sunlit->Tleaf * sunlit_LAI + shaded->Tleaf * shaded_LAI) / LAI;
		//YY average stomatal conductance
		this->conductance = __max(0.0, ((sunlit_gs * sunlit_LAI + shaded_gs * shaded_LAI) / LAI));
	}
	else
	{
		temperature = sunlit->Tleaf;
		this->conductance = 0;
	}

	//static double cumu_photo_gross;
	//cumu_photo_gross = cumu_photo_gross + photosynthesis_gross;
	//cout << " temperature " << photosynthesis_gross << '\n';

	delete sunlit;
	delete shaded;
	delete sun;
	delete light;
}

//Z maintainance respiration, output "g biomass plant^-1"
// based on McCree's paradigm, See McCree(1988), Amthor (2000), Goudriaan and van Laar (1994)
// units very important here, be explicit whether dealing with gC, gCH2O, or gCO2
void RyePlant::calcMaintRespiration(const TWeather& wthr)
{
	const double Q10 = 2.0; // typical Q10 value for respiration, Loomis and Amthor (1999) Crop Sci 39:1584-1596 - could try 1.8
	double dt = initInfo.timeStep * DAYPERMINUTES;
	//	const double maintCoeff = 0.015; // gCH2O g-1DM day-1 at 20C for young plants, Goudriaan and van Laar (1994) Wageningen textbook p 54, 60-61
	const double maintCoeff = 0.018;
	double agefn = (this->GreenLfAreaPlt + 1.0) / (this->LeafAreaPlt + 1.0); // as more leaves senesce maint cost should go down, added 1 to both denom and numer to avoid division by zero. 
	//no maint cost for dead materials but needs to be more mechanistic, SK
	//agefn=1.0;
	double q10fn = pow(Q10, (wthr.airT - 20.0) / 10.0); // should be soil temperature or leaf or combination of use as --> (-stemMass*stem_coef) to reduce
												// total mass. Implement later after testing
												
    //Z stem_coef is not used 
	//double stem_coef = __min(1.0, this->DropLeafMassPlt / this->LeafMassPlt);

	maintRespiration = q10fn * maintCoeff * agefn * (this->ShootMass - this->DropLeafMassPlt - this->DropSheathMassPlt) * dt;// gCH2O dt-1, agefn effect removed. 11/17/14. SK.
}

//Z biomass allocation
//  before doing any computation, C from photosynthesis must be assgined to pools
//  then two steps, 
//		first determine how much biomass is assigned to leaf (including sheath), internode (stem) and root
//		second determine the biomass rates for one unit leaf area, stem length
void RyePlant::calsBiomassAllocation(const TWeather& wthr)
{
	double BiomassReserver2PoolFrac = 0.2;	//Z during BiomassPool send mass out, BiomassReserve also wants to send 20% (refer to apsim wheat)
	double BiomassRootReleaseFrac = 0.8;	//Z root biomass pool want to add 20% storage to root partition if biomass is sufficient

	// ********* step 1 ************************
	//Z organize the budge between reserviors and pools
	//  we will mainly argue the biomass split among the "BiomassSupply", "BiomassDemand", "BiomassPool" and "BiomassReserve".

	double b1 = 2.325152587; // Normalized (0 to 1) temperature response fn parameters, Pasian and Lieth (1990)
						// Lieth and Pasian Scientifica Hortuculturae 46:109-128 1991
						// parameters were fit using rose data -
	double b2 = 0.185418876; // I'm using this because it can have broad optimal region unlike beta fn or Arrhenius eqn
	double b3 = 0.203535650;
	const double Td = 48.6; //High temperature compensation point

	double g1 = 1.0 + exp(b1 - b2 * wthr.airT);
	double g2 = 0.0;
	if (wthr.airT < Td) g2 = 1 - exp(-b3 * (Td - wthr.airT));

	double tmprEffect = g2 / g1; //Z must <1 since g2<1<g1
	double grofac = 1.0 / (5.0 * 60.0 / initInfo.timeStep); // translocation limitation and lag, assume it takes 1 hours to complete, 0.2=5 hrs

	BiomassSupply = 0.0;

	//Z start to rebalance crop biomass budget
	//  always remember these
	//  BiomassPool ---- Checking Account
	//  BiomassReserve ---- Saving Account
	//  BiomassSupply ---- Withdraw

	//Z the first "if statement" handles maintrespiration
	//  the current maintrespiration is computed in "calcMaintRespiration"
	//  which is placed one step higher than the "calsBiomassAllocation", therefore, we can always have the current maintrespiration
	
	//Z feed some biomass from reservior to rapid pool
	BiomassPool += BiomassReserver2PoolFrac * __max(BiomassReserve, 0.0);
	BiomassReserve -= BiomassReserver2PoolFrac * __max(BiomassReserve, 0.0);

	//Z first "pay for" the maint respiration
	if (maintRespiration > 0)
	{
		//Z short term pool can pay for the maintrespiration
		if (BiomassPool > maintRespiration)
		{
			BiomassSupply = maintRespiration;
			BiomassPool -= BiomassSupply;
		}
		//Z short term pool cannot pay for the maintrespiration
		//  but long term pool can pay for the maintrespiration
		else if (BiomassReserve > maintRespiration)
		{
			//Z in this block, we already assume that "BiomassPool <= maintRespiration"
			//  Thus, BiomassReserve is used to satisfy the maintainance respiration demand.
			//  Then BiomassPool is used to replenish BiomassReserve and then BiomassPool is set to 0. 
			//  In other words, BiomassPool is depleted first and then BiomassReserve is used to supplement whatever demand that's left
			BiomassSupply = maintRespiration;
			BiomassReserve -= BiomassSupply;
			BiomassReserve += BiomassPool;	// send remaining C (not enough to meet the demand) in shorterm pool to reserve
			BiomassPool = 0.0;				// empty the BiomassPool
		}
		//Z both short term and long term pools cannot pay for the maintrespiration
		//  then combine them and pay as much as possible
		//  the worst case is that "BiomassPool+BiomassReserve" are exhausted
		else
		{
			BiomassReserve += BiomassPool;
			BiomassPool = 0.0;
			BiomassSupply = __min(BiomassReserve, maintRespiration);
			BiomassReserve -= BiomassSupply;
		}
	}
	//Z no maintRespiration, then there is no BiomassSupply based on it
	else
	{
		BiomassSupply = 0.0;
	}

	//Z the second "if statement" handles the biomass supply to heading, grain filling or relatied issues will also be added to BiomassSupply
	//  Since BiomassDemand==0 based on the current RYESIM design (never enter reproduction stages), 
	//  just reserve this function here
	//Z BiomassDemand does not enter into equations until grain fill
	//  upto this point BiomassSupply only has the possible maintrespiration
	//  or say, maintRepsiration is taken out from the short or/and long term pools
	if (BiomassDemand > 0)
	{
		//Z pay for heading, grain filling or relatied issues using short term pool first and then use long term reserve
		//  but based on the payment, it is proportional to "BiomassPool" rather than "BiomassDemand"
		//  why use "__max(BiomassPool * tmprEffect * grofac, 0)"
		//    rather than "__max(BiomassDemand * tmprEffect * grofac, 0)"
		//  that is because the plant is RICH in biomass and generous
		if (BiomassPool > BiomassDemand)
		{
			double biomassSupplyGrain = __max(BiomassPool * tmprEffect * grofac, 0.0); //CADD from Grant
			BiomassSupply += biomassSupplyGrain;
			BiomassPool -= biomassSupplyGrain;
		}
		//Z if the short term pool BiomassPool not able to pay for the demand, 
		//  then we have to use the long term Reserve pool
		else if (BiomassReserve > BiomassDemand)
		{
			//Z first add BiomassPool to BiomassReserve
			//  BiomassPool always pay first as the short term pool, so it will be zeroed definitely
			//  so by the end, only BiomassReserve will have positive storage
			//  must be a positive storage? yes, because "BiomassReserve > BiomassDemand"
			BiomassReserve += BiomassPool; // BiomassPool negative here, add instead of subtract
			BiomassPool = 0.0;

			//Z this supply is based on "BiomassDemand"
			//  that means the plant can give sufficient biomass but not that generous since it needs to use the saving account
			double biomassSupplyGrain = __max(BiomassDemand * tmprEffect * grofac, 0.0);
			BiomassSupply += biomassSupplyGrain;
			BiomassReserve -= biomassSupplyGrain;
		}
		//Z now neither short term and long term pool cannot satisfy the "BiomassDemand"
		else
		{
			//Z first add BiomassPool to BiomassReserve
			//  BiomassPool always pay first as the short term pool, so it will be zeroed definitely
			//  so by the end, only BiomassReserve will have positive storage
			//  must be a positive storage? yes, because "BiomassReserve > BiomassDemand"
			BiomassReserve += BiomassPool; // BiomassPool negative here, add instead of subtract
			BiomassPool = 0.0;

			//Z this supply is based on "BiomassReserve"
			//  that means the plant may not be even able to supply based on "BiomassDemand"
			//  the plant should really careful about its biomass savings
			double biomassSupplyGrain = __max(__min(BiomassReserve, BiomassDemand) * tmprEffect * grofac, 0.0); //CADD from Grant
			BiomassSupply += biomassSupplyGrain;
			BiomassReserve -= biomassSupplyGrain;
		}
	}
	//Z no reproductive growth and BiomassDemand, then there is no BiomassSupply based on it
	else
	{
		BiomassSupply += 0.0;
	}

	//Z the third "if statement" handles the biomass supply to plant part, since "BiomassDemand" are only for grains
	//  we just need to determine if the plant is willing to give
	//  
	//Z may also need to put rootMassHere
	//  then total ptn mass will be root + shoot 
	double RootShootFractionController = __min(0.95, 0.7 + 0.3 * develop->get_TTd_Joint() / 1000.0);
	double Yg = 1.0; // synthesis efficiency
	BiomassPltGrowth_ptn = __max(BiomassShootGrowth_ptn, 0.0) / RootShootFractionController / Yg;
	{
		//Z ideally after supplying the BiomassPool, 
		//  biomasspool is sufficient to pay for the shoot growth
		if (BiomassPool > BiomassPltGrowth_ptn)
		{
			//Z use BiomassPool to pay for max biomass demand 
			BiomassPool -= BiomassPltGrowth_ptn;
			BiomassSupply += BiomassPltGrowth_ptn;
		}
		//  biomasspool is not sufficient, but the overall biomass is good to pay for the shoot growth
		else if ((BiomassPool+BiomassReserve) > BiomassPltGrowth_ptn)
		{
			//Z plant almost has no short term storage to spend
			BiomassReserve += BiomassPool;
			BiomassPool = 0.0;
			BiomassSupply += BiomassPltGrowth_ptn;
			BiomassReserve -= BiomassPltGrowth_ptn;
		}
		//  not sufficient to pay for the demand
		//  shoot and leaf can still grow due to water pressure
		else
		{
			//Z plant almost has no storage to spend later
			//  if BiomassReserve<= 0 (possibly true), there may be even no C supply to plant growth
			BiomassReserve += BiomassPool;
			BiomassPool = 0.0;
			double BiomassSupply_PltGrowth = __max(__min(BiomassPltGrowth_ptn, BiomassReserve), 0.0);
			BiomassSupply += BiomassSupply_PltGrowth;
			BiomassReserve -= BiomassSupply_PltGrowth;
		}
	}

	//static double cauwzj111;
	//cauwzj111 = BiomassReserve + BiomassPool;
	//cout << cauwzj111 << '\n';
	
	// ********* step 2 ************************
	//Z after get BiomassSupply, and shoot and root parts, 
	//  will need to distribute the values among plant organs
	//  we do not have a lot of growing stage, such as tassel..., so we just run until maturaty

	//Z first determine the supply, and assign 0;
	double shootPart_real = 0.0;
	double rootPart_real = 0.0;

	if (!develop->is_germinInit())
	{
		//Z no seed germinate yet (actually < 1% seeds germinate), pass and wait for the first seed germinating
		return;
	}
	else
	{
		//Z biomass (g) partitioned to shoot and root after the first seed germinates
		//  why not emergence? because when seed germinates, there is already shoot/root growth
		//  start with 0.5,0.5 and then 
		shootPart = __max(0.0, Yg * (RootShootFractionController * (BiomassSupply - maintRespiration)));
		rootPart = __max(0.0, Yg * ((1.0 - RootShootFractionController) * (BiomassSupply - maintRespiration)));
		shootPart_real = shootPart;
		rootPart_real = rootPart;

		if (!develop->is_mature())
		{
			//Z in time step t-1, the value of pcrs is higher than that of pcrl, i.e., use more than expectation, or say already use some from pcrq
			if (wthr.pcrs > rootPart_old)
			{
				//Z the root biomass seems to be overcharged, however, there may exists some root storage,
				//  we add a test to see if the root growth can be paid by the root biomass storage before using the current shoot part
				//
				double rootMassOvercharge = wthr.pcrs - rootPart_old;
				//Z root biomass assignment is sufficient (at least combine with the BiomassRoot Pool), so return leftover to the biomassRoot
				//  then the biomassRoot will release a certion fraction (may be 20%) of the biomass storage to root growth
				if (BiomassRoot >= rootMassOvercharge)
				{
					BiomassRoot -= rootMassOvercharge;
					BiomassRoot = __max(BiomassRoot, 0.0);
					shootPart_real = shootPart;
					rootPart_real = rootPart + BiomassRootReleaseFrac * BiomassRoot;
					BiomassRoot = (1.0 - BiomassRootReleaseFrac) * BiomassRoot;
					rootPart_old = rootPart_real;
				}
				// root biomass storage is not sufficient
				// take some shoot part to pay for the root overcharge
				//Z always make the "rootPart_real" as the root growth values assigned at this time step. 
				else
				{
					// first let the shoot part pay for the overcharge
					if (rootPart > (rootMassOvercharge - BiomassRoot))
					{
						rootPart_real = __max(rootPart - (rootMassOvercharge - BiomassRoot), 0.0);
						shootPart_real = shootPart;
						BiomassRoot = 0.0;
						rootPart_old = rootPart_real;
					}
					// second let the shoot part and root part pay for the overcharge
					else if (shootPart > (rootMassOvercharge - BiomassRoot - rootPart))
					{
						rootPart_real = 0.0;
						shootPart_real = __max(shootPart - (rootMassOvercharge - BiomassRoot - rootPart), 0.0);
						BiomassRoot = 0.0;
						rootPart_old = rootPart_real;
					}
					// eventually previous day overcharge is huge, stop growth
					else
					{
						//Z rootPart_old will leave a negative number
						//  need to be fed back in the following step
						//  during the fed back period, technically no root can grow
						//  the apparent root assignment will be (rootPart_real+BiomassRoot)=0
						rootPart_old = -rootMassOvercharge + BiomassRoot + shootPart + rootPart;
						shootPart_real = 0.0;
						rootPart_real = 0.0;
						BiomassRoot = 0.0;
					}
				}
			}
			//Z root biomass assignment is sufficient, so return leftover to the biomassRoot
			//  then the biomassRoot will release a certion fraction (may be 20%) of the biomass storage to root growth
			else
			{
				double rootMassLeftOver = rootPart_old - wthr.pcrs;
				BiomassRoot += rootMassLeftOver;
				rootPart_real = rootPart + BiomassRootReleaseFrac * BiomassRoot;
				BiomassRoot = (1.0 - BiomassRootReleaseFrac) * BiomassRoot;
				rootPart_old = rootPart_real;
			}

			//Z fill the shoot part again
			if (shootPart_real >= shootPart)
			{
				//Z shoot growth is satisfied, do nothing
			}
			else if ((shootPart_real + BiomassPool) >= shootPart)
			{
				BiomassPool -= (shootPart - shootPart_real);
				shootPart_real = shootPart;
			}
			else if ((shootPart_real + BiomassPool + BiomassReserve) >= shootPart)
			{
				BiomassReserve -= (shootPart - shootPart_real - BiomassPool);
				shootPart_real = shootPart;
				BiomassPool = 0.0;
			}
			else
			{
				shootPart_real += (BiomassReserve + BiomassPool);
				BiomassPool = 0.0;
				BiomassReserve = 0.0;
			}

			//Z partition carbon to plant organs
			leafPart = shootPart_real * __max(BiomassLeafGrowth_ptn / BiomassShootGrowth_ptn, 0.0);
			internodePart = __max(shootPart_real - leafPart, 0.0);
				
			//Z judge if plant give leaf or internode too much
			if (leafPart <= BiomassLeafGrowth_ptn) {
				// potential amount >= actual amount, then use actual amount
			}
			else {
				// potential amount < actual amount, so the plant gives leaf more than leaf can take
				BiomassPool += ((leafPart - BiomassLeafGrowth_ptn) / Yg);
				leafPart = BiomassLeafGrowth_ptn;
			}

			if (internodePart <= BiomassIntrNodeGrowth_ptn) {
				// potential amount >= actual amount, then use actual amount
			}
			else {
				// potential amount < actual amount, so the plant gives leaf more than leaf can take
				BiomassPool += ((internodePart - BiomassIntrNodeGrowth_ptn) / Yg);
				internodePart = BiomassIntrNodeGrowth_ptn;
			}
		}
	}
	//??? cauwzj test if root biomass is correctly assigned
	//rootPart_real = 0.05;

	ActuralRootBiomassAssignment_PCRL = rootPart_real;
	ActuralShootBiomassAssignment = shootPart_real;
}

//Z Nitrogen allocation
//		first get the availiable nitrogen either from seeds or from plant root uptake
//		second assign the nitrogen to plant organs
void RyePlant::calsNitrogenAllocation()
{
	//Z first determine the nitrogen quantity that can be allocated
	//  similar to the biomass allocation
	ShootNitrogenAvailiableAllocation = 0.0;
	RootNitrogenAvailiableAllocation = 0.0;
	leafPartNitrogen = 0.0;
	internodePartNitrogen = 0.0;

	if (!develop->is_germinInit())
	{
		return;
	}
	else 
	{
		if (!develop->is_emerge())
		{
			//Z nitrogen solely comes from seed at 3.4%
			//  measured in mg (1000.0 converts g to mg)
			//  NitrogenPool += leafPart * MAX_N_PCT * 10.0;
			NitrogenPool = SeedNitrogenMass;
			
			//Z first root
			if (NitrogenPool >= NitrogenRootGrowth_ptn)
			{
				RootNitrogenAvailiableAllocation = NitrogenRootGrowth_ptn;
				RootNitrogenMassPlt += RootNitrogenAvailiableAllocation;
				NitrogenPool -= RootNitrogenAvailiableAllocation;
			}
			else
			{
				RootNitrogenAvailiableAllocation = NitrogenPool;
				RootNitrogenMassPlt += RootNitrogenAvailiableAllocation;
				NitrogenPool = 0.0;
			}

			//Z second shoot
			if (NitrogenPool >= NitrogenShootGrowth_ptn)
			{
				ShootNitrogenAvailiableAllocation = NitrogenShootGrowth_ptn;
				NitrogenPool -= ShootNitrogenAvailiableAllocation;
			}
			else
			{
				ShootNitrogenAvailiableAllocation = NitrogenPool;
				NitrogenPool = 0.0;
			}
		}
		else
		{
			//Z nitrogen comes from root uptake during this hour and the N released from plant aging
			NitrogenPool += (HourlyNitrogenSoilUptake + LeafNitrogenReleasePlt + SheathNitrogenReleasePlt + IntrNodeNitrogenReleasePlt);
			
			//Z first root
			if (NitrogenPool >= NitrogenRootGrowth_ptn)
			{
				RootNitrogenAvailiableAllocation = NitrogenRootGrowth_ptn;
				RootNitrogenMassPlt += RootNitrogenAvailiableAllocation;
				NitrogenPool -= RootNitrogenAvailiableAllocation;
			}
			else
			{
				RootNitrogenAvailiableAllocation = NitrogenPool;
				RootNitrogenMassPlt += RootNitrogenAvailiableAllocation;
				NitrogenPool = 0.0;
			}

			//Z second shoot
			if (NitrogenPool >= NitrogenShootGrowth_ptn)
			{
				ShootNitrogenAvailiableAllocation = NitrogenShootGrowth_ptn;
				NitrogenPool -= ShootNitrogenAvailiableAllocation;
			}
			else
			{
				ShootNitrogenAvailiableAllocation = NitrogenPool;
				NitrogenPool = 0.0;
			}
		}
	}

	//Z second assign nitrogen to shoot organs
	if (!develop->is_germinInit())
	{
		return;
	}
	else
	{
		if (!develop->is_mature())
		{
			//Z partition carbon to plant organs

			//Z internode is growing, both leaves and internode receive nitrogen
			//  acceleration, leaf and sheath grow at the same speed
			//  pre-accel, sheath grow at 0.5 speed of the leaf grow
			leafPartNitrogen = ShootNitrogenAvailiableAllocation * __max(NitrogenLeafGrowth_ptn / NitrogenShootGrowth_ptn, 0.0);
			internodePartNitrogen = __max(ShootNitrogenAvailiableAllocation - leafPartNitrogen, 0.0);
		}
	}
}

//Z after allocating the biomass, set the allocated biomass "AND N"
//  back to plant organs for the growth
void RyePlant::calsSetMass()
{
	//Z first compute the assignment rate
	//  for example, we have the potential leaf area incremental and the biomass assigned for leaves,
	//  then, it will be easy to compute the leaf biomass per leaf area
	
	//Z note that we do not need to distinct the "acceleration" growth or not, 
	//  just use "leafPart/(leaf area incremental)" although some leafPart will go to sheath
	//  the leaf/sheath ratio will be adjusted in the RyeLeaf class
	LeafBiomassAssignmentRate = __max(leafPart / PtnLeafMassIncreasePlt, 0.0);
	LeafNitrogenAssignmentRate = __max(leafPartNitrogen / PtnLeafNitrogenMassIncreasePlt, 0.0);

	IntrNodeBiomassAssignmentRate = __max(internodePart / PtnIntrNodeMassIncreasePlt, 0.0);
	IntrNodeNitrogenAssignmentRate = __max(internodePartNitrogen / PtnIntrNodeNitrogenMassIncreasePlt, 0.0);

	//Z we have the "mass rate"
	//  then assign the mass rate to grow plant in a recursive way
	this->PltCallTillerRecursiveMassDistribution(this->mainstem);
}

//Z Recursive UPDATE and SUMMARY of leaf mass and internode mass
//     recursive is done for all the tillers
//  
void RyePlant::PltCallTillerRecursiveMassDistribution(RyeTiller* rt)
{
	rt->RyeTillerSingleMassDistribution(LeafBiomassAssignmentRate, LeafNitrogenAssignmentRate,
		IntrNodeBiomassAssignmentRate, IntrNodeNitrogenAssignmentRate);
	//Z recursive process, starting with the mainstem
	for (int ii = 0; ii < MAXLEAFNUM; ii++)
	{
		//Z first name this subtiller "subrt"
		RyeTiller* subrt = rt->get_subtiller(ii);
		if (subrt != NULL)
		{
			PltCallTillerRecursiveMassDistribution(subrt);
		}
	}
}

void RyePlant::calcRed_FRedRatio(const TWeather& weather)
// this function calculates an estimate of the Red to Far red light ratio from sunlit and shaded ET. This 
// ration is used to estimate the effects of plant density on leaf expansion and LAI. 
// A daily mean ratio is calculated. We use a 3 parameter sigmoid function to model the effect
{
	//double Xo=0.43, B=0.05, A=1.2;
	//double Xo=0.6, B=0.13, A=2.0; original
	//double Xo=0.9, B=0.43, A=2.0;
	double Xo = 0.85, B = 0.65, A = 2.0;
	double dt = initInfo.timeStep * DAYPERMINUTES;
	double C2_effectTemp;
	//First set counter to 0 if it is the beginning of the day. 
	if (abs(weather.time) < 0.0001)
	{// have to rename C2_effect to Light_effect
		//Zhu et al. Journal of Experimental Botany, Vol. 65, No. 2, pp. 641–653, 2014
		C2_effectTemp = exp(-(SunlitRatio - Xo) / B);
		C2_effect = __min(1.0, A / (1.0 + C2_effectTemp));
		develop->set_shadeEffect(C2_effect);
		SunlitRatio = 0.0;
	}
	else
	{
		// calculate from emergence
		if (sunlit_LAI / (sunlit_LAI + shaded_LAI) > 0.05) { SunlitRatio += sunlit_LAI / (sunlit_LAI + shaded_LAI) * dt; }
		else { SunlitRatio += 1.0 * dt; }
	}
}