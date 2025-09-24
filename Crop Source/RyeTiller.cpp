#include "stdafx.h"
#include "RyeTiller.h"
#include "weather.h"
#include "initinfo.h"
#include <cmath>
#include <algorithm>
#define MINUTESPERDAY (24*60)
#define DAYPERMINUTES 0.00069444444
#define DAYPERHOUR 0.041666666667
#define MAXLEAFNUM 20
#define PHYLLOCHRON 106.0
#define CO2_MW 44.0098
#define C_MW 12.011
#define CH2O_MW 30.03
#define MAX_N_PCT 3.5
#define MIN_N_PCT 0.35

using namespace std;

//Z constructor, need to identify if it is the mainstem
// ADD THIS: New overloaded constructor
RyeTiller::RyeTiller(int n, int o, int cr, bool m, RyeDevelopment* dv, double tillerLivingFrac)
{
	//Z initialize the growing stage
	rank = n;
	order = o;
	cumuRank = cr;
	mainstem = m;
	develop = dv;
	mainstemInitialization = false;	//Z a special case for mainstem, immediately grow two leaves
	death_2_finalize = false;

	//Z Growing Stage
	//  for the whole plant, such that each tiller follow that stage
	singleRidge = develop->is_singleRidge();//Z the singleRidge stage is reached for the whole plant, if tiller branching after single ridge, keep it,
	singleRidge_already = singleRidge;		//Z combine with "singleRidge", s.t. tiller death operation on single ridge date only called once
	singleRidge_kill = false;
	born_after_singleRidge = false;
	if (singleRidge && singleRidge_already)
	{
		born_after_singleRidge = true;
	}

	elongationStart = false;			    //Z the enlongation is initiated
	elongationStart_already = false;		//Z if elongation already occurs, combine with "elongationStart" to trace the first entry of elongation
	living = true;						    //Z the tiller is living, reserve for counting tiller death, living=1, death=0
	livingFrac = tillerLivingFrac;		    //Z DIFFERENT FROM LEAF,INTERNODE, leaf and internode livingfrac follows the tiller, while tiller and subtiller livingfrac management independently
										    //  but the initial subtiller livingfrac follows the tiller livingfrac, after initialization, subtiller will be on its own, EVEN the parent tiller is dead.
	                                        //  This setting is meant to mimic "perenial cases" where in the next year, all tillers are new subtillers 

	livingFrac_ini = tillerLivingFrac;
	force_to_death_current_step = false;    //Z marker for tiller death due to external/ambient reasons, force to die
	TillerAge = 0.0;
	physAge = 0.0;
	TTd_since_singleRidge = 0.0;

	//Z initialize tiller temperature
	TmprTiller = 25.0;
	Cur_TTd = 0.0;
	StressTiller = 1.0;
	StressTillerMin = 1.0;

	N_effect = 1.0;
	water_effect = 1.0;
	shade_effect = 1.0;
	tmpr_effect = 1.0;
	tmpr_effect_terminal = 1.0;
	Cold_Time = 0.0;
	
	//Z 0 for N stress and 1 for water stress
	LeafStressData_N = 1.0;
	LeafStressData_W = 1.0;
	LeafStressDataCount = 0.0;
	LeafStressDataCountMax = 72.0;
	NewTillerIndicator = false;

	for (int ii = 0; ii < MAXLEAFNUM; ii++)
	{
		TillerLfPhyllchron[ii] = (double)ii * PHYLLOCHRON + PHYLLOCHRON;
		SubTiller[ii] = NULL;
		SubInterNode[ii] = NULL;
		SubLeaf[ii] = NULL;
	}
	FlagLfPhyllchron = TillerLfPhyllchron[MAXLEAFNUM - 1];

	//Z integer counts for the tillers
	//  tiller is just a "box", leaf and internode are the real growing parts
	LeafNum = 0;					//Z Leaf number = green leaf number + dropped leaf number
	GreenLfNum = 0;
	DropLfNum = 0;
	InterNodeNum = 0;
	SubTillerNum = 0;

	//Z Cumulating Leaf features for this tiller
	TillerLfArea = 0.0;							//Z tiller leaf area cm^2 = green leaf area + senescent leaf area
	TillerSeneLfArea = 0.0;						//Z tiller senecent leaf area, including drop leaf area
	TillerGreenLfArea = 0.0;					//Z tiller green leaf area
	TillerGreenLfLength = 0.0;				    //Z tiller green leaf length, a little bit special here because it is the sum (not average) of green leaf length on this tiller
	TillerGreenLfWidth = 0.0;					//Z tiller green leaf width, a little bit special here because it is the sum (not average) of green leaf width on this tiller
	TillerDropLfArea = 0.0;						//Z tiller drop leaf area must be senecent leaf, but senecent (part) of the leaf may not dropped
	TillerLfMass = 0.0;							//Z even dropped, leaf mass still hold the number
	TillerGreenLfMass = 0.0;					//Z tiller sum of no scenescent (or dropped) leaf mass
	TillerGreenSheathMass = 0.0;				//Z tiller sum of no scenescent (or dropped) sheath mass
	TillerSheathMass = 0.0;						//Z even dropped, sheath mass still hold the number
	TillerDropLfMass = 0.0;						//Z after dropped, = leaf mass, otherwise 0
	TillerDropSheathMass = 0.0;					//Z after dropped, = sheath mass, otherwise 0
	TillerLfNitrogenMass = 0.0;					//Z leaf N mass in mg
	TillerSheathNitrogenMass = 0.0;				//Z sheath N mass in mg
	TillerLfNitrogenRelease = 0.0;				//Z during senecent, N released from leaf, mg
	TillerSheathNitrogenRelease = 0.0;			//Z during senecent, N released from sheath, mg
	TillerLfMassIncrease_ptn = 0.0;				//Z potential leaf biomass increase (not include sheath, need compute depending on elongationStart later) for this tiller, g
	TillerLfNitrogenMassIncrease_ptn = 0.0;		//Z potential leaf nitrogen mass increase (not include sheath, need compute depending on elongationStart later) for this tiller, mg

	//Z partition of the non-mobile N in plant senescent and dead organs:
	//  Dead: total amount of N that cannot be recycled for new organs
	//  Yellow: N in senescent but not dropped leaf
	//  Dropped: N in dropped leaf
	//    Dead=Yellow+Dropped

	//Z Different from biomass computation, 
	//  For Biomass, Senescent biomass includes dropped biomass
	//  For N, Yellow and Dropped are exclusive

	//Z Different from internode N
	//  For internode, there is no senescent function, so internode just has "dead nitrogen mass"
	TillerLfNitrogenDeadMass = 0.0;				//Z Nitrogen put in dead tissues (leaf), not removable but carried to residue
	TillerSheathNitrogenDeadMass = 0.0;			//Z Nitrogen put in dead tissues (sheath), not removable but carried to residue
	TillerLfNitrogenYellowMass = 0.0;
	TillerSheathNitrogenYellowMass = 0.0;
	TillerLfNitrogenDroppedMass = 0.0;
	TillerSheathNitrogenDroppedMass = 0.0;

	//Z Cumulating Internode features for this tiller
	//Z Assume internode will not die unless the whole tiller dies
	TillerIntrNodeLength = 0.0;					//Z sum of internode length in cm
	TillerIntrNodeMass = 0.0;					//Z internode mass in g
	TillerIntrNodeLivingMass = 0.0;				//Z internode living portion of the mass, g
	TillerIntrNodeNitrogenMass = 0.0;			//Z internode N mass in mg
	TillerIntrNodeNitrogenRelease = 0.0;		//Z internode N release in mg
	TillerIntrMassIncrease_ptn = 0.0;			//Z potential internode biomass increase, g
	TillerIntrNitrogenMassIncrease_ptn = 0.0;	//Z potential internode nitrogen mass increase, mg

	TillerIntrNitrogenDeadMass = 0.0;           //Z Nitrogen put in dead tissues (internode), not removable but carried to residue

	/*Z the new code is mass based
	    remove all the variable baesd on area or lenght
	TillerLfAreaIncrease_ptn = 0.0;				//Z potential leaf area increase, used to distribute leaf (sheath) biomass assignemnt from assimilation
	TillerIntrNodeLengthIncrease_ptn = 0.0;		//Z internode increase at this step, use to distribute internode biomass assignment from assimilation
    */

	//Z marker for the elongation internode ranks
	elongation_first = -1;
	elongation_last = -1;
	rank_2_elongation = -1;

	//----------------------------------------------------------
	/*Z initialize the representative value
	    this changes "one single tiller" to "representative tiller"
		NOTE that other (double/float) variables such as leaf area, mass, internode length, mass are all default as representative values
		or say, the difference between representative/single values are in leaf and internode classes
	*/
	LeafNum_Rep = 0.0;
	GreenLfNum_Rep = 0.0;
	DropLfNum_Rep = 0.0;
  #if _DEBUG
	cout << "Tiller Emerging: Order = " << order << ", Rank = " << rank << ", CumuRank = " << cumuRank << "\n";
  #endif
}
RyeTiller::RyeTiller(int n, int o, int cr, bool m, RyeDevelopment* dv, double tillerLivingFrac, double seedMass)
	: RyeTiller(n, o, cr, m, dv, tillerLivingFrac)  // Delegate to original constructor
{
	this->seedMass = seedMass;  // Override seedMass for mainstem
	rank = n;
	order = o;
	//TillerLfMass = seedMass*0.50;							//Z even dropped, leaf mass still hold the number
	//illerGreenLfMass=TillerLfMass;					//Z tiller sum of no scenescent (or dropped) leaf mass
}
//Z destructor, need to be recursive
//  in all the models, death does not mean destruction,
//  dead organ still holds informations
//Z Therefore, this destructor will not involve mass transfer issues
//  Recall for any recursive process, only call it from the root element, aka the mainstem
RyeTiller::~RyeTiller()
{
	//Z first kill all the leaf and internode
	for (int ii = 0; ii < MAXLEAFNUM; ii++)
	{
		if (SubLeaf[ii] != NULL)
		{
			delete SubLeaf[ii];
			SubLeaf[ii] = NULL;
		}
		if (SubInterNode[ii] != NULL)
		{
			delete SubInterNode[ii];
			SubInterNode[ii] = NULL;
		}
	}
	//Z need to kill all the subtiller in a recursive way
	//  therefore, do not destruct tiller during model running, even that tiller is dead
	//  because its subtiller or sub-subtiller may still functioning
	for (int ii = 0; ii < MAXLEAFNUM; ii++)
	{
		if (SubTiller[ii] != NULL)
		{
			delete SubTiller[ii];
			SubTiller[ii] = NULL;
		}
	}
}

/*Z special case of mainstem initializaiton between germinationand emergence
    apsim initialized 2 leaves, therefore in the function
    a. initial 2 leaves and internodes
    b. revise the TillerLfPhyllchron order
*/
/*Z why this is an isolated function separate from the constructor?
    because constructor for mainstem is in the constructor of plant
    and this function is called between germination and emergence
*/
/*Z because this is the mainstem, "livingFrac=1.0" for all its leaves and internodes initially
*/
void RyeTiller::MainstemInitialize()
{
	//Z this function only called once for mainstem
	//  but we still put if-statment for protection
	if (living && (!singleRidge_kill))
	{
		if (mainstem && !mainstemInitialization)
		{
			mainstemInitialization = true;
			// first 1 leaves
			// recall parameter, index of leaf on the tiller; bool if tiller is mainstem; development obj; initial livingFrac
			SubLeaf[0] = new RyeLeaf(1, this->order, mainstem, develop, 1.0, seedMass);
			LeafNum = 1;
			GreenLfNum = 1;
			DropLfNum = 0;
			// the associate internode
			// recall parameter, index of internode on the tiller; bool if tiller is mainstem; development obj; initial livingFrac
			SubInterNode[0] = new RyeInterNode(1, this->order, mainstem, develop, 1.0);
			InterNodeNum = 1;

			TillerLfPhyllchron[0] = 0.0;
			for (int ii = 1; ii < MAXLEAFNUM; ii++)
			{
				TillerLfPhyllchron[ii] = (double)(ii) * PHYLLOCHRON;
			}
			FlagLfPhyllchron = TillerLfPhyllchron[MAXLEAFNUM - 1];
		}
	}
}

// *************** GROWTH AND DEVELOP *******************************************
// As OOP, one tiller should only care about itself, leaves and internodes
// multiple tiller recursive operation should be done in a upper class -> RyePlant

//Z NEW leaf and internode only on this tiller
//     include new leaf, internode, and subtillers
//     because new organs will emerge, the tiller has to be living
// 
//Z this is the first function to call when tiller is operatored
void RyeTiller::RyeTillerSingleMorph(void)
{
	if(living)
	{
		//Z get the TTd for this time step and adding to the age of this tiller
		//     TTd is the thermal time for this 1 hour period
		TmprTiller = develop->get_Tcur();
		Cur_TTd = develop->get_TTd();
		elongationStart = develop->is_startEnlongation();
		physAge += Cur_TTd;
		TillerAge += develop->get_TimeStep() * DAYPERMINUTES;

		//Z new leaf, tiller, and stem
		//  initiate leaf, internode and tillers based on Phyllchron
		if (physAge >= TillerLfPhyllchron[LeafNum])
		{
			NewTillerIndicator = true;

			// first emerge leaf
			SubLeaf[LeafNum] = new RyeLeaf(LeafNum+1, this->order, mainstem, develop, this->livingFrac);
			LeafNum += 1;

			// second emerge sub-tiller
			// if leaf appears, then tiller at pre-specified position will also appear
			// first argument means "mainstem=false", that is correct because subtiller cannot be the mainstem
			//		1. if elongation starts, no new tillers
			//		2. if tiller killed during singleRidge, no subtiller from it but existing subtillers will continue
			if ((!elongationStart) && (!singleRidge_kill))
			{
				if (mainstem && LeafNum >= 3)
				{
					if (SubTiller[LeafNum - 3] == NULL)
					{
						/*Z recall the input parameters, bool if sub - tiller is the mainstem(of course not)
						                                 develop obj
														 initial sub-tiller livingfrac follow the tiller livingfrac
						*/
						double tillerStress = this->livingFrac;
						float xx = __min(tillerStress, sqrt(LeafStressData_W * water_effect));
						tillerStress = __max(__min(xx, 1.0), 0.4);
						xx = __min(tillerStress, sqrt(LeafStressData_N * N_effect));
						tillerStress = __max(__min(xx, 1.0), 0.4);
#if _DEBUG
						cout << LeafStressData_N << "," << LeafStressData_W << "," << tillerStress << ",";
#endif
						SubTiller[LeafNum - 3] = new RyeTiller(LeafNum - 2, this->order + 1, this->rank + LeafNum - 2, false, develop, tillerStress);
						SubTillerNum += 1;
					}
				}
				else if ((!mainstem && LeafNum >= 2) && (!singleRidge_kill))
				{
					if (SubTiller[LeafNum - 2] == NULL)
					{
						/*Z recall the input parameters, bool if sub - tiller is the mainstem(of course not)
														 develop obj
														 initial sub-tiller livingfrac follow the tiller livingfrac
						*/
						double tillerStress = this->livingFrac;
						float xx = __min(tillerStress, sqrt(LeafStressData_W * water_effect));
						tillerStress = __max((__min(xx, 1.0)), 0.4);
						xx = __min(tillerStress, sqrt(LeafStressData_N * N_effect));
						tillerStress = __max(__min(xx, 1.0), 0.4);
#if _DEBUG
						cout << LeafStressData_N << "," << LeafStressData_W << "," << tillerStress << ",";
#endif
						SubTiller[LeafNum - 2] = new RyeTiller(LeafNum - 1, this->order + 1, this->rank + LeafNum - 1, false, develop, tillerStress);
						SubTillerNum += 1;
					}
				}
			}

			// third emerge internode
			// recall the inputs, internode index on this tiller; bool if the tiller is mainstem; develop obj; initial internode livingfrac following its tiller
			SubInterNode[InterNodeNum] = new RyeInterNode(InterNodeNum + 1, this->order, mainstem, develop, this->livingFrac);
			InterNodeNum += 1;
		}
	}
	//Z TODO growth stage, HAUN, FEEKERS, ZODAK
}

//Z DEATH of leaf and internode only for this tiller
//     not a destruction of the tiller, internode or leaf at this step
//     because those classes should exists to hold drop leaf and dead internode mass
void RyeTiller::RyeTillerSingleDeath(void)
{
	//Z prevent kill the same tiller twice
	if (living)
	{
		/*Z tiller (and its leaf and internode) dies within the "tiller update" function
		    therefore, although we set living = false and livingFrac = 0.0 here
		    it will still call leaf and internode update to finish the death
		*/
		/*Z even the tiller is dead, 
		    single tiller summary function will be always called to account living AND dead tissues
		*/
		living = false;
		livingFrac = 0.0;
		force_to_death_current_step = true;
		//Z call leaf and internode to die
		for (int ii = 0; ii < MAXLEAFNUM; ii++)
		{
			//Z leaf death, which will also set leaf's "force_to_death_current_step=true", 
			//  s.t., later leaf area update and leaf senescence can be called normally.
			if (SubLeaf[ii] != NULL)
			{
				SubLeaf[ii]->RyeLeafDead();
			}
			//Z internode death, which will also set internode's "force_to_death_current_step=true", 
			//  s.t., later internode update (and possibly future senescence) can be called normally.
			if (SubInterNode[ii] != NULL)
			{
				//Z tiller controls which internode (or rank) should elongate
				SubInterNode[ii]->RyeInterNodeDeath();
			}
		}
	}
}

//Z GROWTH of leaf and internode only for this tiller
//     include leaf, internode, but DO NOT sub-tiller, 
//     sub-tiller manages its own growth 
void RyeTiller::RyeTillerSingleUpdate(void)
{
	TmprTiller = develop->get_Tcur();
	double TmprTiller_movingAve = develop->get_Tave();
	//Z this function call the update function for each leaf and internode in this tiller
	//  shape growth only
	//  no mass growth, mass growth will be based on the assignment later
	if (living)
	{
		/*Z adjust tiller living fraction based on external stress
		    external stress are based on nitrogen and water
		    in SHOOTGRO, this process is done in "Evalst.for"
		*/
		singleRidge = develop->is_singleRidge();

		/*Z this is just compariable but it is simplified, Z does not know how to correctly use it, and hence use the livingfrac functions
		    1. first compute nitrogen, water stress, and temperature stresses for the tiller
		    2. then determine if a new minimal stress value (new highest stress occurs)
		    3. use the new stress value to reduce the livingfraction
		    4. if livingfraction < 0.01 and it is not the mainstem, then call tiller death function
        */
		double tillerNitrogenMass = TillerLfNitrogenMass + TillerSheathNitrogenMass + TillerIntrNodeNitrogenMass;
		double tillerMass = TillerGreenLfMass + TillerGreenSheathMass + TillerIntrNodeLivingMass;
		double tillerNitrogenMassContent = 0.10 * tillerNitrogenMass / tillerMass;	//Z note that N is in mg and biomass is in g, and the content is based on %
		double tillerCriticalNitrogen = __max(MIN_N_PCT, tillerNitrogenMassContent);
		if (tillerMass < 0.0001)
		{
			tillerCriticalNitrogen = MAX_N_PCT;
		}
		
		//Z to give more Discrimination at high N end
		N_effect = __max(0.0, (2.0 / (1.0 + exp(-1.0 * (tillerCriticalNitrogen - MIN_N_PCT))) - 0.6));
		N_effect = __min(1.0, N_effect);
		N_effect = __max(0.1, N_effect);

		water_effect = develop->get_tillerWaterEffectAve();


		shade_effect = 1.0;
		//Z shaded effects on tiller, marked based on R:FR ratio, only affect plant with more than 4 tillers.
		if (develop->get_plantTillerNum() >= 4) { shade_effect = develop->get_shadeEffect(); }


		Cold_Time = develop->get_ColdTime();
		if (singleRidge)
		{
			TTd_since_singleRidge = develop->get_TTd_sinceSingleRidge();
			tmpr_effect = __min(__max(1.0 - (1.0 - tmpr_effect_terminal) / PHYLLOCHRON * (TTd_since_singleRidge - PHYLLOCHRON), tmpr_effect_terminal), 1.0);
		}
		tmpr_effect = __min(tmpr_effect, 1.0);
		tmpr_effect = __max(tmpr_effect, 0.01);
	
		StressTiller = 1.0;
		StressTiller = tmpr_effect * water_effect * N_effect * shade_effect;

		if (StressTiller < StressTillerMin)
		{
			livingFrac = livingFrac * (StressTiller / StressTillerMin);
			StressTillerMin = StressTiller;
			if (mainstem)
			{
				livingFrac = __max(livingFrac, 0.90);
			}
			else
			{
				if (livingFrac < 0.01) { RyeTillerSingleDeath(); }
			}
		}

		//Z call tiller death on single-ridge stage
		//  when single-ridge occurs, tiller with leaf number < 4 is killed
		//  This will not prevent new tiller branching based on SHOOTGRO

		//Z tiller emerge after singleRidge should be kept
		//  so we use a marker "singleRidge_kill" to name which tiller should be killed
		
		if (singleRidge && (!singleRidge_already))
		{
			int singleRidgeKill_LeafNum = 3;
			double gdddiff= develop->get_TTd_Plant()- develop->get_TTd_Joint();
			double cold_time_frac = develop->get_ColdTimeRatioJoint();
			if (Cold_Time <= 30.0) 
			{ 
				singleRidgeKill_LeafNum = 0; 
				tmpr_effect_terminal = 1.0;
			}
			else 
			{ 
				singleRidgeKill_LeafNum = 1; 
				tmpr_effect_terminal = 2.0 / (1.0 + exp(4.0 * __max(cold_time_frac - 0.1, 0.0)));
			}

			singleRidge_already = singleRidge;
			TTd_since_singleRidge = 0.0;
			if (LeafNum < singleRidgeKill_LeafNum && (!mainstem)) { singleRidge_kill = true; }
		}
		if (singleRidge && singleRidge_already && born_after_singleRidge)
		{
			double cold_time_frac = develop->get_ColdTimeRatioJoint();
			if (Cold_Time <= 30.0)
			{
				tmpr_effect_terminal = 1.0;
			}
			else
			{
				tmpr_effect_terminal = 2.0 / (1.0 + exp(5.0 * __max(cold_time_frac - 0.1, 0.0)));
			}
			born_after_singleRidge = false;
		}

		//Z technically, when crop is at single ridge, some tiller should be killed.
		//  this statement means, if the tiller should be killed, let it not happen immediately, but in a period of, say, 1 PHYLLOCHRON
		if (singleRidge_kill && livingFrac >= 0.01)
		{
			livingFrac = livingFrac * (1.0 * PHYLLOCHRON - TTd_since_singleRidge - Cur_TTd) / (1.0 * PHYLLOCHRON - TTd_since_singleRidge);
			if (livingFrac < 0.01) { RyeTillerSingleDeath(); }
		}

		//Z internode elongation or not
		//  only run this part once to determine the fist and the last internode to elongate for this tiller
		//  internode status is for the whole plant, but the internode that can grow is based on each tiller's situation
		if (elongationStart && (!elongationStart_already))
		{
			if (Cold_Time >= 50.0)
			{
				tmpr_effect_terminal = 0.1;
			}

			//Z first get gdd parameters from development
			double mingdf = develop->get_TTd_FlagLf_min();
			double gddj = develop->get_TTd_Joint();
			double gddp = develop->get_TTd_Plant();
			double gdde = develop->get_TTd_Elong();

			//Z how many leaves can grow on this tiller in future
			int lnum = __max((int)((mingdf - TillerLfPhyllchron[LeafNum - 1]) / PHYLLOCHRON), 0.0);
			//Z flag leaf may not be the MAXLEAF, it depends on ambient and plant structure
			//  see "gddf" in the shootgro
			FlagLfPhyllchron = TillerLfPhyllchron[LeafNum - 1] + (double)(lnum + 1) * PHYLLOCHRON;

			for (int ii = LeafNum + lnum + 1; ii < MAXLEAFNUM; ii++)
			{
				//Z no more leaf can grow beyond current leaf num + future leaf number "lnum", set Phyllchron for leaf emergence to arbitrary large
				TillerLfPhyllchron[ii] = 9999.0;
			}

			//Z first (lowest) internode that can elongate
			//  the idea is if the current leaf is complete, then use "LeafNum - 1"
			//  if the current leaf is still expanding, the use "LeafNum - 2"
			//  finally, the first two internode should never elongate
			if ((TillerLfPhyllchron[LeafNum - 1] + PHYLLOCHRON) <= (gddp - gddj + gdde))
			{
				elongation_first = LeafNum - 1;
			}
			else
			{
				elongation_first = LeafNum - 2;
			}
			elongation_first = __max(elongation_first, 3);

			//Z last (highest) internode that can elongate
			elongation_last = elongation_first + lnum;
			elongation_last = __min(elongation_last, MAXLEAFNUM);

			//Z make sure this block only run once.
			elongationStart_already = elongationStart;
		}

		//Z under elongation status, control which internode should grow
		if (elongationStart)
		{
			//Z first collect plant gdd number need to use
			double gddp = develop->get_TTd_Plant();

			//case 1: If the last internode on the culm (= peduncle) has completed elongation, after initiation of the last leaf (= flag leaf).
			//		  It is 3*PCHRON because internodes grow for 1 phyllochron after a lag of 2 phyllochrons from when the associated leaf appears.
			if (elongation_last > 0)
			{
				//Z the current leaf is already the last one, so the even the last internode loses its chance 
				if ((TillerLfPhyllchron[LeafNum - 1] + PHYLLOCHRON) >= FlagLfPhyllchron)
				{
					if (gddp >= (TillerLfPhyllchron[elongation_last - 1] + 3.0 * PHYLLOCHRON))
					{
						//Z all the internodes that can elongate should already finish their elongation growth
						rank_2_elongation = -1;
					}
				}
			}

			//case 2: Grow the peduncle (or last internode to elongate on the culm) if appropriate.
			if (elongation_last > 0)
			{
				//Z the current leaf is already the last one, the last one can still grow
				if ((TillerLfPhyllchron[LeafNum - 1] + PHYLLOCHRON) >= FlagLfPhyllchron)
				{
					if ((gddp >= (TillerLfPhyllchron[elongation_last - 1] + 2.0 * PHYLLOCHRON)) && (LeafNum > 1))
					{
						//Z at this time, the flagleaf is growing and the last internode should grow
						rank_2_elongation = LeafNum;
					}
				}
			}

			//case 3: Grow the penultimate internode on the culm if appropriate.
			if (elongation_last - 1 > 0)
			{
				if (LeafNum > 1)
				{
					//Z gdd is sufficient for plant growth and the corresponding leaf exists
					if ((gddp >= (TillerLfPhyllchron[elongation_last - 2] + 2.0 * PHYLLOCHRON)) && SubLeaf[elongation_last - 2] != NULL)
					{
						//Z at this time, the last-but-one leaf could grow and the last-but-two internode can grow 
						rank_2_elongation = LeafNum - 1;
					}
				}
			}

			//case 4: If only two leaves have appeared on the culm at the beginning of today, no internode will grow today
			//        That tiller is not ready to elongate
			if (LeafNum <= 2)
			{
				if (gddp >= (TillerLfPhyllchron[LeafNum - 1] + PHYLLOCHRON))
				{
					rank_2_elongation = -1;
				}
			}

			//case 5: The currently growing internode is the one corresponding to the currently growing leaf - 2 and is not the peduncle or penultimate internode.
			if (LeafNum >= 3)
			{
				rank_2_elongation = LeafNum - 2;
			}
		}
		
		//Z call leaf and internode to grow
		for (int ii = 0; ii < MAXLEAFNUM; ii++)
		{
			//Z leaf area update for existing leaves
			if (SubLeaf[ii] != NULL)
			{
				SubLeaf[ii]->set_LfLivingFrac(this->livingFrac);
				SubLeaf[ii]->RyeLeafAreaUpdate();
			}
			//Z internode update for existing internodes
			if (SubInterNode[ii] != NULL)
			{
				//Z tiller controls which internode (or rank) should elongate
				SubInterNode[ii]->set_ElongationRank(rank_2_elongation);
				SubInterNode[ii]->set_IntrLivingFrac(this->livingFrac);
				SubInterNode[ii]->RyeInterNodeLengthUpdate();
			}
		}
	}
	else
	{
		//Z not living anymore
		//  still need to call once to finalize suborgan, leaf and internode, calculations for once
		if (death_2_finalize == false)
		{
			for (int ii = 0; ii < MAXLEAFNUM; ii++)
			{
				//Z leaf area update for existing leaves
				if (SubLeaf[ii] != NULL)
				{
					SubLeaf[ii]->set_LfLivingFrac(0.0);
					SubLeaf[ii]->RyeLeafAreaUpdate();
				}
				//Z internode update for existing internodes
				if (SubInterNode[ii] != NULL)
				{
					//Z tiller controls which internode (or rank) should elongate
					SubInterNode[ii]->set_IntrLivingFrac(0.0);
					SubInterNode[ii]->RyeInterNodeLengthUpdate();
				}
			}
			death_2_finalize = true;
		}
	}
}


/*Z SUM of leaf and internode only for this tiller
       include leaf, internode, but DO NOT sub-tiller, 
       sub-tiller manages its own growth
  NOTE all the quantities here mush be "representative" values that will be reported to plant later. 
*/
void RyeTiller::RyeTillerSingleSummary(void)
{
	//Z LEAF
	//  reset everything that needs to be summed
	//  cumulating Numbers for this tiller
	TillerLfArea = 0.0;
	TillerSeneLfArea = 0.0;
	TillerGreenLfArea = 0.0;
	TillerGreenLfLength = 0.0;
	TillerGreenLfWidth = 0.0;
	TillerDropLfArea = 0.0;
	
	TillerLfMass = 0.0;
	TillerSheathMass = 0.0;
	TillerGreenLfMass = 0.0;
	TillerGreenSheathMass = 0.0;
	TillerDropLfMass = 0.0;
	TillerDropSheathMass = 0.0;
	TillerLfNitrogenMass = 0.0;
	TillerSheathNitrogenMass = 0.0;
	TillerLfNitrogenRelease = 0.0;
	TillerSheathNitrogenRelease = 0.0;

	TillerLfMassIncrease_ptn = 0.0;
	TillerLfNitrogenMassIncrease_ptn = 0.0;

	TillerLfNitrogenDeadMass = 0.0;
	TillerSheathNitrogenDeadMass = 0.0;
	TillerLfNitrogenYellowMass = 0.0;
	TillerSheathNitrogenYellowMass = 0.0;
	TillerLfNitrogenDroppedMass = 0.0;
	TillerSheathNitrogenDroppedMass = 0.0;

	//Z leaf number
	//  update integer as index, update double as representative values
	LeafNum = 0;
	GreenLfNum = 0;
	DropLfNum = 0;

	LeafNum_Rep = 0.0;
	GreenLfNum_Rep = 0.0;
	DropLfNum_Rep = 0.0;

	//Z INTERNODE
	//  reset everything that needs to be summed
	//  cumulating Numbers for this tiller
	TillerIntrNodeLength = 0.0;
	TillerIntrNodeMass = 0.0;
	TillerIntrNodeLivingMass = 0.0;
	TillerIntrNodeNitrogenMass = 0.0;
	TillerIntrNodeNitrogenRelease = 0.0;

	TillerIntrMassIncrease_ptn = 0.0;
	TillerIntrNitrogenMassIncrease_ptn = 0.0;

	TillerIntrNitrogenDeadMass = 0.0;

	for (int ii = 0; ii < MAXLEAFNUM; ii++)
	{
		if (SubLeaf[ii] != NULL) 
		{
			//Z leaf number
			//  when leaf is initialized, leaf number should be kept, 
			//  but living fraction is dynamic, so living fraction can be used for green leaf number
			LeafNum_Rep += SubLeaf[ii]->get_LfLivingFracIni();
			GreenLfNum_Rep += SubLeaf[ii]->get_LfLivingFrac();

			//Z integer leaf number as index
			LeafNum += 1;
			GreenLfNum += 1;
			if (SubLeaf[ii]->is_Dropped()) { GreenLfNum -= 1; }

			//Z leaf area 
			TillerLfArea += SubLeaf[ii]->get_leafArea();
			TillerGreenLfArea += SubLeaf[ii]->get_greenArea();
			TillerGreenLfLength += SubLeaf[ii]->get_greenLfLength();
			TillerGreenLfWidth += SubLeaf[ii]->get_greenLfWidth();
			TillerSeneLfArea += SubLeaf[ii]->get_senescentArea();
			//Z even the leaf is not dropped, from the statistical or representative perspective,
			//  there may be some dropped leaf due to tiller death
			TillerDropLfArea += SubLeaf[ii]->get_dropLfArea();

			//Z leaf/sheath mass
			TillerLfMass += SubLeaf[ii]->get_leafMass();
			TillerSheathMass += SubLeaf[ii]->get_sheathMass();
			TillerGreenLfMass += SubLeaf[ii]->get_greenLfMass();
			TillerGreenSheathMass += SubLeaf[ii]->get_greenSheathMass();
			TillerDropLfMass += SubLeaf[ii]->get_dropLfMass();
			TillerDropSheathMass += SubLeaf[ii]->get_dropSheathMass();
			TillerLfNitrogenMass += SubLeaf[ii]->get_leafNitrogenMass();
			TillerSheathNitrogenMass += SubLeaf[ii]->get_sheathNitrogenMass();
			TillerLfNitrogenRelease += SubLeaf[ii]->get_leafNitrogenRelease();
			TillerSheathNitrogenRelease += SubLeaf[ii]->get_sheathNitrogenRelease();

			TillerLfNitrogenDeadMass += SubLeaf[ii]->get_leafNitrogenDeadMass();
			TillerSheathNitrogenDeadMass += SubLeaf[ii]->get_sheathNitrogenDeadMass();
			TillerLfNitrogenYellowMass += SubLeaf[ii]->get_leafNitrogenYellowMass();
			TillerSheathNitrogenYellowMass += SubLeaf[ii]->get_sheathNitrogenYellowMass();
			TillerLfNitrogenDroppedMass += SubLeaf[ii]->get_leafNitrogenDroppedMass();
			TillerSheathNitrogenDroppedMass += SubLeaf[ii]->get_sheathNitrogenDroppedMass();

			TillerLfMassIncrease_ptn += SubLeaf[ii]->get_potentialLfMassIncrease();
			TillerLfNitrogenMassIncrease_ptn += SubLeaf[ii]->get_potentialLfNitrogenMassIncrease();
						
		}
		if (SubInterNode[ii] != NULL)
		{
			TillerIntrNodeLength += SubInterNode[ii]->get_InterNodeLength();
			TillerIntrNodeMass += SubInterNode[ii]->get_InterNodeMass();
			TillerIntrNodeLivingMass += SubInterNode[ii]->get_InterNodeLivingMass();
			TillerIntrNodeNitrogenMass += SubInterNode[ii]->get_InterNodeNitrogenMass();
			TillerIntrNodeNitrogenRelease += SubInterNode[ii]->get_InterNodeNitrogenRelease();

			TillerIntrNitrogenDeadMass += SubInterNode[ii]->get_InterNodeNitrogenDeadMass();

			TillerIntrMassIncrease_ptn += SubInterNode[ii]->get_potentialInterNodeMassIncrease();
			TillerIntrNitrogenMassIncrease_ptn += SubInterNode[ii]->get_potentialInterNodeNitrogenMassIncrease();
		}
	}

	//Z only count total leaf num and green leaf num
	//  drop leaf num = total leaf num - green leaf num
	DropLfNum = LeafNum - GreenLfNum;
	DropLfNum_Rep = LeafNum_Rep - GreenLfNum_Rep;

	//Z Leaf Stresses Moving Average
	if (NewTillerIndicator)
	{
		//Z recall 0 for N stress and 1 for water stress (when tiller emerge, leaf is young, so only use "water_effect_expand" from leaf)
		LeafStressData_N = 1.0;
		LeafStressData_W = 1.0;
		LeafStressDataCount = 0.0;
		NewTillerIndicator = false;
	}
	LeafStressDataCount += 1.0;
	if (mainstem && LeafNum >= 2)
	{
		//Z the tiller from leaf number - 3 is already emerged, now we need to compute the next one
		if (LeafStressDataCount <= LeafStressDataCountMax)
		{
			LeafStressData_N = (LeafStressData_N * LeafStressDataCount + SubLeaf[LeafNum - 2]->get_leafNitrogenEffect()) / (LeafStressDataCount + 1.0);
			LeafStressData_W = (LeafStressData_W * LeafStressDataCount + SubLeaf[LeafNum - 2]->get_leafWaterEffectExpand()) / (LeafStressDataCount + 1.0);
		}
		else
		{
			LeafStressData_N = (LeafStressData_N * (LeafStressDataCountMax - 1.0) + SubLeaf[LeafNum - 2]->get_leafNitrogenEffect()) / LeafStressDataCountMax;
			LeafStressData_W = (LeafStressData_W * (LeafStressDataCountMax - 1.0) + SubLeaf[LeafNum - 2]->get_leafWaterEffectExpand()) / LeafStressDataCountMax;
		}
	}
	else if (!mainstem && LeafNum >= 1)
	{
		//Z the tiller from leaf number - 1 is already emerged, now we need to compute the next one, i.e., tiller is "producing" at the current leaf
		if (LeafStressDataCount <= LeafStressDataCountMax)
		{
			LeafStressData_N = (LeafStressData_N * LeafStressDataCount + SubLeaf[LeafNum - 1]->get_leafNitrogenEffect()) / (LeafStressDataCount + 1.0);
			LeafStressData_W = (LeafStressData_W * LeafStressDataCount + SubLeaf[LeafNum - 1]->get_leafWaterEffectExpand()) / (LeafStressDataCount + 1.0);
		}
		else
		{
			LeafStressData_N = (LeafStressData_N * (LeafStressDataCountMax - 1.0) + SubLeaf[LeafNum - 1]->get_leafNitrogenEffect()) / LeafStressDataCountMax;
			LeafStressData_W = (LeafStressData_W * (LeafStressDataCountMax - 1.0) + SubLeaf[LeafNum - 1]->get_leafWaterEffectExpand()) / LeafStressDataCountMax;
		}
	}
	LeafStressData_N = __max(__min(LeafStressData_N, 1.0), 0.01);
	LeafStressData_W = __max(__min(LeafStressData_W, 1.0), 0.01);
}


/*Z SINGLE MASS (NITROGEN) DISTRIBUTION
       assign biomass and N to each of the leaves and internodes on this tiller
       DO NOT include subtillers, they are on their own
       NOTE that biomass is in "g", nitrogen is in "mg"
       to be specific, lfBiomassRate (include sheath)   g per unit leafarea cm^2
                       lfNitrogenRate (include sheath)  g per unit leafarea cm^2
                       intrnodeBiomassRate              g per unit intrnodeLength cm
                       intrnodeNitrogenRate             g per unit intrnodeLength cm
*/
//Z the biomass/nitrogen rates will be adjusted based on living fraction in leaf and internode classes
void RyeTiller::RyeTillerSingleMassDistribution(double lfBiomassRate, double lfNitrogenRate, double intrnodeBiomassRate, double intrnodeNitrogenRate)
{
	for (int ii = 0; ii < MAXLEAFNUM; ii++)
	{
		//Z if leaf exists and grow, add biomass and N to it
		if (SubLeaf[ii] != NULL)
		{
			SubLeaf[ii]->RyeLeafMassUpdate(lfBiomassRate, lfNitrogenRate);
		}
		//Z if internode starts to grow with a positive potential length increase, add biomass and N to it
		if (SubInterNode[ii] != NULL)
		{
			SubInterNode[ii]->RyeInterNodeMassUpdate(intrnodeBiomassRate, intrnodeNitrogenRate);
		}
	}
}

