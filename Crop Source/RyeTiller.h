#pragma once

#ifndef _RYE_TILLER_H_
#define _RYE_TILLER_H_
#include "weather.h"
#include "RyeThermalTime.h"
#include "RyeDevelopment.h"
#include "RyeInterNode.h"
#include "RyeLeaf.h"

#define MAXLEAFNUM 20

class RyeTiller
{
public:
public:
	// Original constructor
	RyeTiller(int rank, int order, int cumuRank, bool mainstem, RyeDevelopment* dv, double livingFrac);

	// ADD THIS: Overloaded constructor with seedMass
	RyeTiller(int rank, int order, int cumuRank, bool mainstem, RyeDevelopment* dv, double livingFrac, double seedMass);
    ~RyeTiller();
	//Z update the growth of the current tiller
	void RyeTillerSingleMorph(void);
	void RyeTillerSingleDeath(void);
	void RyeTillerSingleUpdate(void);
	void RyeTillerSingleSummary(void);
	void RyeTillerSingleMassDistribution(double, double, double, double);

	//Z global update based on a recursive method
	//  embed in the ryeplant class, because ryeplant class is a "upper class" and ryetiller class is a "lower class"
	//  plant class has the authority to control the tiller class
	//  while tiller class can control (get_information from) even lower classes, such as leaf and internode.


	//***** Ordinary IO function for tiller-leaf-internode *********
	double get_physAge() { return physAge; }
	double get_FlagLfPhyllchron() { return FlagLfPhyllchron; }
	RyeTiller* get_subtiller(int ii) { return this->SubTiller[ii]; }

	//***** mainstem development *********
	bool get_mainstemInitilization() { return mainstemInitialization; }
	bool is_living() { return living; }
	void MainstemInitialize();

	//***** Group of functions used in Plant class recursive process *********************
	//Z In the recursive process
	//  IO function that cumulate leaf and internode features tiller by tiller
	//  those functions will be called by the tillers using their parentPlant pointer
	//  some are of int type and some are of double type
	//  some are for leaves and some are for internode (stem)

//Z tiller living fraction
	//  can be used as "representative" tiller number if living, not necessarily integer because it report a statistical value
	double get_tillerLivingFrac() { return livingFrac; }
	//  should be a boolean function, but want to use it as an integer number of living tiller in a plant
	int is_tillerLiving() 
	{ 
		if (living) { return 1; } 
		else { return 0; }
	}
	void set_tillerLivingFrac(double x)
	{
		if (!living) { livingFrac = 0.0; }
		else { livingFrac = x; }
	}

//Z for all the leaves
	//  numbers
	double get_tillerLeafNum() { return this->LeafNum_Rep; }
	double get_tillerGreenLeafNum() { return this->GreenLfNum_Rep; }
	double get_tillerDropLeafNum() { return this->DropLfNum_Rep; }

	//  double variables
	double get_tillerLeafArea() { return this->TillerLfArea; }
	double get_tillerGreenLeafArea() { return this->TillerGreenLfArea; }
	double get_tillerGreenLeafLength() { return this->TillerGreenLfLength; }
	double get_tillerGreenLeafWidth() { return this->TillerGreenLfWidth; }
	//  always remember that DropLeaf is included in Seneleaf, so DropLeaf <= Seneleaf and = at leaf drop stage
	double get_tillerSeneLeafArea() { return this->TillerSeneLfArea; }
	double get_tillerDropLeafArea() { return this->TillerDropLfArea; }
	
	double get_tillerLeafMass() { return this->TillerLfMass; }
	double get_tillerSheathMass() { return this->TillerSheathMass; }
	double get_tillerDropLeafMass() { return this->TillerDropLfMass; }
	double get_tillerDropSheathMass() { return this->TillerDropSheathMass; }
	double get_tillerLeafNitrogenMass() { return this->TillerLfNitrogenMass; }
	double get_tillerSheathNitrogenMass() { return this->TillerSheathNitrogenMass; }

	//Z nitrogen mass redistribution when organ is dead, 
	//  partition for leaf is more dedicated, since leaf has senescence function
	//  need to include that for internode?
	double get_tillerLeafNitrogenDeadMass() { return this->TillerLfNitrogenDeadMass; }
	double get_tillerSheathNitrogenDeadMass() { return this->TillerSheathNitrogenDeadMass; }
	double get_tillerLeafNitrogenRelease() { return this->TillerLfNitrogenRelease; }
	double get_tillerSheathNitrogenRelease() { return this->TillerSheathNitrogenRelease; }
	double get_tillerLeafNitrogenYellowMass() { return this->TillerLfNitrogenYellowMass; }
	double get_tillerSheathNitrogenYellowMass() { return this->TillerSheathNitrogenYellowMass; }
	double get_tillerLeafNitrogenDroppedMass() { return this->TillerLfNitrogenDroppedMass; }
	double get_tillerSheathNitrogenDroppedMass() { return this->TillerSheathNitrogenDroppedMass; }
	
	
	//Z potential growth should always based on mass
	//  rather than area/length
	double get_tillerLeafMassIncrease_ptn() { return this->TillerLfMassIncrease_ptn; }
	double get_tillerLeafNitrogenMassIncrease_ptn() { return this->TillerLfNitrogenMassIncrease_ptn; }
	//double get_tillerLeafAreaIncrease_ptn() { return this->TillerLfAreaIncrease_ptn; }

//Z for all the internode
	double get_IntrNodeLength() { return this->TillerIntrNodeLength; }
	double get_IntrNodeMass() { return this->TillerIntrNodeMass; }
	double get_IntrNodeNitrogenMass() { return this->TillerIntrNodeNitrogenMass; }
	double get_IntrNodeNitrogenRelease() { return this->TillerIntrNodeNitrogenRelease; }
	double get_IntrNodeNitrogenDeadMass() { return this->TillerIntrNitrogenDeadMass; }

	//Z potential growth should always based on mass
	//  rather than area/length
	double get_tillerIntrMassIncrease_ptn() { return this->TillerIntrMassIncrease_ptn; }
	double get_tillerIntrNitrogenMassIncrease_ptn() { return this->TillerIntrNitrogenMassIncrease_ptn; }
	//double get_IntrNodeLengthIncrease_ptn() { return this->TillerIntrNodeLengthIncrease_ptn; }

private:

	//Z trace the main plant development status
	RyeDevelopment* develop;

	//Z LINKED-LIST for tiller hierarchy
	int SubTillerNum;
	RyeTiller* SubTiller[MAXLEAFNUM];
	//Z Internode and leaf are derived from the tiller
	int InterNodeNum;
	RyeInterNode* SubInterNode[MAXLEAFNUM];
	int LeafNum;
	int GreenLfNum;
	int DropLfNum;
	double seedMass;
	RyeLeaf* SubLeaf[MAXLEAFNUM];
    public:
    // Add this method to RyeTiller to access a leaf by index
    RyeLeaf* get_subLeaf(int ii) { return this->SubLeaf[ii]; }
	int rank;						//Z number of the tiller at its parent tiller (maybe mainstem) counted from bottom
	int order;						//Z the order of tiller, 0 for mainstem, 1 for the tiller branching from the mainstem, 2 for the tiller branching from the order 1 tiller...
	int cumuRank;					//Z rank of parent tiller + rank of this tiller, maybe better measure the distance of the tiller from root
	                                //  recall that "parent tiller"'s cumurank is the grand parent tiller + rank of this parent tiller
	                                //  so cumuRank is a recursive sum

	bool singleRidge;				// the singleRidge stage is reached for the whole plant
	bool singleRidge_already;		// combine with "singleRidge", s.t. tiller death operation on single ridge date only called/determined once
	bool singleRidge_kill;          // this tiller should be killed at single ridge stage, a boolean marker for this tiller
	bool born_after_singleRidge;
	
	bool elongationStart;			// receive start elongation from the plant and development class
	bool elongationStart_already;	// if elongation already occurs, combine with "elongationStart" to trace the first entry of elongation
	
	bool living;					// reserve for counting tiller death, living=1, death=0
	double livingFrac;				// livingFrac of the tiller, 
	double livingFrac_ini;			// livingfrac when this tiller is initialized
	bool force_to_death_current_step;// marker for tiller death due to external/ambient reasons, force to die, can only be true on one step when tiller death function is called
	
	bool mainstem;					// if this tiller is the mainstem or not
	bool mainstemInitialization;	// a special case for mainstem, immediately grow two leaves

	bool death_2_finalize;          // a computational marker, when tiller is killed or dead, use this to update its leaf and internode for the last time

	int elongation_first;			// the first internode for elongation
	int elongation_last;			// the last internode for elongation
	int rank_2_elongation;			// mark the rank of internode, which will elong at the current time step
	
	double TillerLfPhyllchron[MAXLEAFNUM];

	//***** Age and temperature *********
	double TillerAge; // chronological age based on time, days
	double physAge;   // gdd age (in reference to endGrowth and lifeSpan, days)
	//Z Thermal time for this step (1 hour)
	double Cur_TTd;
	double FlagLfPhyllchron;
	//Z Thermal time start from some critical time points
	double TTd_since_singleRidge;

	//***** Environmental Factors *********
	double TmprTiller;		// organ temperature, C
	double StressTiller;	//Z Tiller water and nitrogen stress, 1 means no stress and 0 means full stress
	double StressTillerMin; //Z the minimal stress value ever recorded

	double N_effect;
	double water_effect;
	double shade_effect;
	double tmpr_effect;
	double tmpr_effect_terminal;
	double Cold_Time;


	//Z this block will determine (parent leaf stress, which determine the "one-step" tiller emerge fraction)
	//  the stress from previous steps (>5 days or 120 hours) will be forgot
	//  0 for N stress and 1 for water stress (when tiller emerge, leaf is young, so only use "water_effect_expand" from leaf)
	// 
	//Z because we only account ONE tiller branching for a given tiller, 
	//  therefore, it is sufficient to ONLY compute the values within "that phyllochron" and then reset for the next tiller
	double LeafStressData_N;
	double LeafStressData_W;
	double LeafStressDataCount;
	double LeafStressDataCountMax;
	bool NewTillerIndicator;		//Z this will indicate the emerge of new tiller, and then the stress factor will be reset and move to next tiller position
	
	//-----------------------------------------------------------------
	/*Z below this point, all the variables should be understood as "Representative" Values
	    SubTillerNum, InterNodeNum, LeafNum, GreenLfNum and DropLfNum have specific functions such as index
		so we need to define "double" version of those parameters for the "Representative" values
	*/

	double LeafNum_Rep;
	double GreenLfNum_Rep;
	double DropLfNum_Rep;

	//***** Leaf Group (g biomass, cm^2 area) ONLY for this tiller *********
	double TillerLfArea;
	double TillerGreenLfArea;
	double TillerGreenLfLength;
	double TillerGreenLfWidth;
	double TillerSeneLfArea;
	double TillerDropLfArea;
	double TillerLfMass;
	double TillerSheathMass;
	double TillerGreenLfMass;
	double TillerGreenSheathMass;
	double TillerDropLfMass;
	double TillerDropSheathMass;
	double TillerLfNitrogenMass;
	double TillerSheathNitrogenMass;
	double TillerLfNitrogenRelease;
	double TillerSheathNitrogenRelease;
	double TillerLfMassIncrease_ptn;           //Z should be enlarged 1.5* or 2* to include leaf & sheath 
	double TillerLfNitrogenMassIncrease_ptn;   //Z should be enlarged 1.5* or 2* to include leaf & sheath 

	double TillerLfNitrogenDeadMass;
	double TillerSheathNitrogenDeadMass;
	double TillerLfNitrogenYellowMass; 
	double TillerSheathNitrogenYellowMass; 
	double TillerLfNitrogenDroppedMass; 
	double TillerSheathNitrogenDroppedMass; 

	//***** InterNode Group (g biomass) ONLY for this tiller *********
	double TillerIntrNodeLength;
	double TillerIntrNodeMass;
	double TillerIntrNodeLivingMass;
	double TillerIntrNodeNitrogenMass;
	double TillerIntrNodeNitrogenRelease;
	double TillerIntrMassIncrease_ptn;
	double TillerIntrNitrogenMassIncrease_ptn;

	double TillerIntrNitrogenDeadMass;

	//Z potential growth should based on mass
	//  rather than area/length
	//double TillerLfAreaIncrease_ptn;
	//double TillerIntrNodeLengthIncrease_ptn;
};

#endif