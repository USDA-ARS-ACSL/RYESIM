#pragma once
#ifndef _RYE_LEAF_
#define _RYE_LEAF_
#include "weather.h"
#include "RyeDevelopment.h"

class RyeLeaf
{
public:
	RyeLeaf(int rank, int order, bool mainstem, RyeDevelopment* dv, double livingFrac); // take current leaf rank and total leaf number to calculate potentialArea
	// New overloaded constructor with seedmass
	RyeLeaf(int rank, int order, bool mainstem, RyeDevelopment* dv, double LfLivingFrac, double seedMass);

	~RyeLeaf();
	//Z IO bool functions for leaf stages
	//  unused function are listed at the end and comment-out
	// 	DO NOT REMOVE, for those stages may be used in future
	//  leaf death and dropping are off by one time step, say, one-hour
	bool is_Dropped() { return dropped; }
	//bool is_Initiated() { return initiated; }
	//bool is_Appeared() { return appeared; }
	//bool is_Growing() { return growing; }
	//bool is_Mature() { return mature; }
	//bool is_Aging() { return aging; }
	//bool is_Dead() { return dead; }

//Z IO living fraction, initialize, setting, and getting
/*Z "get_LfLivingFrac()" is the current fraction, can be considered as the active leaf or greenleaf
	"get_LfLivingFracIni()" is the livingfraction when leaf is initialized, so it can be considered as the representative value of the leafnum, since leafnum include everything even die earlier.
*/
	double get_LfLivingFrac() { return livingFrac; }
	double get_LfLivingFracIni() { return livingFrac_ini; }
	void set_LfLivingFrac(double x)
	{
		if (dropped)
		{
			livingFrac = 0.0;
		}
		else
		{
			livingFrac = x;
		}
	}

	//Z IO Functions
	// 	ONLY REPRESENTATIVE values are output and seen by upperlevel structures, e.g., tiller and plant
	// 	unused function are listed at the end and comment-out temporarily
		//leaf scale, area or gdd
	double get_leafArea() { return LeafArea_Rep; }					//leaf area includes every thing, i.e., green + senescent
	double get_greenArea() { return GreenLfArea_Rep; }
	double get_greenLfLength() { return GreenLfLength_Rep; }
	double get_greenLfWidth() { return GreenLfWidth_Rep; }
	double get_senescentArea() { return SenescentLfArea_Rep; }		//recall that senescentArea include dropArea, 
	double get_dropLfArea() { return DropLfArea_Rep; }				//dropped leaf area is a subset of senescent leaf area

	//for mass
	//total mass includes dropped mass
	//dropped tissue has 0 N mass, so N released back to pool
	double get_leafMass() { return LeafMass_Rep; }
	double get_sheathMass() { return SheathMass_Rep; }
	double get_greenLfMass() { return GreenLfMass_Rep; }
	double get_greenSheathMass() { return GreenSheathMass_Rep; }
	double get_dropLfMass() { return DropLfMass_Rep; }
	double get_dropSheathMass() { return DropSheathMass_Rep; }
	double get_leafNitrogenMass() { return LeafNitrogenMass_Rep; }
	double get_sheathNitrogenMass() { return SheathNitrogenMass_Rep; }
	double get_leafNitrogenRelease() { return LeafNitrogenRelease_Rep; }
	double get_sheathNitrogenRelease() { return SheathNitrogenRelease_Rep; }
	double get_leafNitrogenDeadMass() { return DeadLeafNitrogenMass_Rep; }
	double get_sheathNitrogenDeadMass() { return DeadSheathNitrogenMass_Rep; }
	double get_leafNitrogenYellowMass() { return YellowLeafNitrogenMass_Rep; }
	double get_sheathNitrogenYellowMass() { return YellowSheathNitrogenMass_Rep; }
	double get_leafNitrogenDroppedMass() { return DroppedLeafNitrogeMass_Rep; }
	double get_sheathNitrogenDroppedMass() { return DroppedSheathNitrogeMass_Rep; }
	//double get_leafNitrogenContent() { return LeafNitrogenContent; } //Z remove this function, for all calculation/sum should be mass based

	/*Z  leaf potential growth for biomass/nitrogen allocation
	     mass increase only for leaf
	     sheath mass increase = 0.5* or 1.0* leaf mass increase depending on accleration
	     therefore, some computation will shown in "plant.cpp"
	*/

	double get_potentialLfMassIncrease() { return ptnLfMassIncrease_Rep; }
	double get_potentialLfNitrogenMassIncrease() { return ptnLfNitrogenMassIncrease_Rep; }
	//Z gradually remove potential area increases, and then rely on mass incremental
	//double get_potentialLfAreaIncrease() { return ptnLfAreaIncrease_Rep; }

	//leaf stress factor
	double get_leafNitrogenEffect() { return N_effect; }
	double get_leafWaterEffectExpand() { return water_effect_expand; }

/*Z  leaf growth and senescent operations
     first compute geometrical growth
                   potential mass and nitrogen requests
     then "RyeLeafMassUpdate" is the actural assigned biomass
*/
	void RyeLeafAreaUpdate();
	void LeafAreaExpand();
	void LeafAreaSenescence();
	void LeafLivingFractorAdjustment();
	void RyeLeafMassUpdate(double biomassIncomeRate, double nitrogenIncomeRate);
	void RyeLeafDead();
	
//Z  leaf environmental functions
//   tmpr
//   leaf/canopy water potential (now in development.cpp)
	double TmprEffect(double Tmpr);
//	double LWPeffect(double predawn_psil, double threshold);

private:
	
	RyeDevelopment* develop;

	//Z Growing Stage
	int rank;			//Z number of the leaf at that tiller (mainstem) counted from root to top
	int order;          //Z the order of tiller that this leaf grows on, 0 means mainstem, 1 means the tillerS from mainstem
	bool initiated;     //Z following are growing status
	bool appeared;
	bool growing; 
	bool mature;
	bool aging;
	bool dead;
	bool dropped;
	bool mainstem;						//Z mainstem leaf yes/no, equivalent to "order == 0"
	bool acceleration;					//Z yes: sheath growth=leaf growth; no: sheath growth=0.5 leaf growth
	bool force_to_death_current_step;	//Z when "RyeLeafDead()" is called, used to bypass the leaf senescence function
	bool is_singleRidge;

	/*Z living fraction
	    trace the living fraction changes
	    for one leaf, it is living or dead; 
	    for a representative plant, need to consider statistical effects -- therefore, a leaf can be partially living.
	    when living fraction changes, caused by tiller living fraction, those portion should be considered as dead and dropped
	    therefore, computation of droping should be done and N should be recycled
	*/
	double livingFrac;         
	double livingFrac_old;
	double livingFrac_ini;
	
	//Z Growing time/gdd measure
	double LeafAge;				//Z chronological ordinary time, day
	double physAge;				//Z gdd time (in reference to endGrowth and lifeSpan, days)
	double GDD2mature;			//Z gdd at mature, after fully expansion
	double activeAge;			//Z gdd since mature, the total functioning time period
	double seneAge;				//Z gdd since the end of activeAge, period for senescence
	double Cur_TTd;				//Z gdd current incremental (after photoperiod and vernalisation adjustment)
	double seneFrac;			//Z current senefraction, based on the matured (max) leaf area
	double seneFrac_old;		//Z old senefraction for the previous iteration

	//Z Leaf scale and area recorder (cm or cm^2)
	double LeafLength;			//Z leaf length in cm
	double LeafWidth;			//Z leaf width in cm
	double GreenLfLength;		//Z leaf length when it is not dropped, if dropped, will be 0
	double GreenLfWidth;		//Z leaf width when it is not dropped, if dropped, will be 0
	double maxLfLength;			//Z max leaf length based on leaf rank and if it is at the mainstem
	double maxLfWidth;			//Z max leaf width based on leaf rank and if it is at the mainstem
	double ptnLfLength;			//Z potential leaf length, the values after one step of potential growth
	double ptnLfWidth;			//Z potential leaf width, the values after one step of potential growth
	double ptnLfLengthIncrease;	//Z potential leaf length increase at the current time step
	double ptnLfWidthIncrease;	//Z potential leaf width increase at the current time step
	
	double A_LW;				//Z shape factor "area = shape factor * length * width", dimensionless
	double LeafArea;			//Z total leaf area changing over time, cm^2
	double MaturedLfArea;		//Z leaf area at the end of expansion, i.e., the actual max leaf area at one time (different from LeafArea)
	double SenescentLfArea;		//Z aged leaf area changing over time, cm^2
	double DropLfArea;			//Z drop leaf area, the leaf area at drop (dead) time point
	double GreenLfArea;			//Z green (active) leaf area, i.e., leaf area before aging, and part of the leaf area after aging
	double maxLeafArea;			//Z max leaf area based on leaf rank and if it is at the mainstem
	double ptnLeafArea;			//Z potential leaf area, the values after one step of potential growth
	double ptnLfAreaIncrease;	//Z potential leaf area increase, the values after one step of potential growth
	double ptnLfAreaDecrease;	//Z potential leaf area decrease, the values after one step of aging

	/*Z leaf biomass + sheath biomass (g)
	    dropped mass is still mass, so even if the leaf is dropped, it should be counted for leafmass
	    therefore, if needed, LeafMass-DropLfMass = attached leafmass 
	    (for single leaf, either attached leafmass 0 or = leafmass, since leaf is either attached or dropped)
	*/
	double slw = 0.0045;		    //Z specific leaf weight (g/cm^2)
	double slw_max;					//Z leaf can uptake more carbon if carbon is sufficient
	double slw_min = 0.0011;  //  min slw value when plant is small, refer to the slw range in (DOI: 10.1080/01904167.2015.1017051)
	double slw_run;                 //  running value of slw;
	double LeafMass;				//Z leaf biomass no matter it is green or not, dropped or not
	double SheathMass;				//Z sheath biomass no matter it is green or not, dropped or not
	double GreenLfMass;				//Z green (living) leaf mass, exclude senescent (and dropped) mass g
	double GreenSheathMass;			//Z green (living) sheath mass, exclude senescent (and dropped) mass g
	double DropLfMass;				//Z dropped leaf biomass, 0 or = leafmass depending on leaf drop or not g
	double DropSheathMass;			//Z dropped sheath biomass, 0 or = sheathmass depending on leaf drop or not mg zN
	double LeafBiomassIncrease;		//Z biomass increase after plant assign biomass based on photosynthesis and the biomass request
	double SheathBiomassIncrease;   //Z biomass increase after plant assign biomass based on photosynthesis and the biomass request
	double ptnLfMassIncrease;		//Z potential leaf mass increase due to leaf area increase at this time step, and biomass is not fully assigned in previous steps
									//		maximum leaf mass = "potential leaf area"*slw
									//		actual leaf mass = LeafMass
									//      potential leaf mass increase = maximum leaf mass - current leaf mass
									//  potential sheath mass increase = 0.5 or 1.0 leaf mass increase depending on accleration
	double seedMass;  // Seedmass to pass when intiating the first leaves

	//Z Leaf Sheath Nitrogen (mg): income part
	double LeafNitrogenContent;		//Z leaf N% fraction in % (e.g., 3.5% as the max value)
	double SheathNitrogenContent;	//Z sheath N% fraction in % (e.g., 3.5% as the max value)
	double LeafNitrogenMass;		//Z leaf nitrogen mass mg N
	double SheathNitrogenMass;		//Z sheath nitrogen mass mg N
	double LeafNitrogenIncrease;	//Z leaf nitrogen mass increase after nitrogen assignment mg N
	double SheathNitrogenIncrease;	//Z sheath nitrogen mass increase after nitrogen assignment mg N
	double ptnLfNitrogenMassIncrease;//Z potential leaf N mass increase due to leaf area increase at this time step, and N is not fully assigned in previous steps
									//		maximum leaf N mass = "potential leaf area"*slw*3.5% (need unit adjustment)
									//		actual leaf N mass = LeafMass
									//      potential leaf N mass increase = maximum leaf N mass - current leaf N mass
									//  potential sheath N mass increase = 0.5 or 1.0 leaf N mass increase depending on accleration

	//Z Leaf Sheath Nitrogen (mg): output after organ senecent or dropped
	double LfNitrogenReleaseThreshold;	//Z the min leaf N mass content, below that, no N will be released and all N goes to dropped tissue, assume to be 0.5%
	double LfNitrogenMaxReleasePtge;	//Z the leaf N release percentage, default to be 80%, that means 80% of leaf N will be recycled/remobilized
	                                    //  but need to make conpromise with the LfNitrogenReleaseThreshold, that means temproery percentge values should be defined
	
	/*Z Reduce_Single : single leaf N output due to senecense and dropping, need to divide
	*   Decline : at that step, reduced N goes to residue and immobilized
	*   Release : at that step, remobilized and can reused for new plant organs
	*   DeadNitrogen: cumulated N mass in the senecense/dead plant organs
	*/
	double LeafNitrogenReduce_Single;
	double SheathNitrogenReduce_Single;
	double LfNitrogenDecline;			//Z stepwise dead portion of leaf nitrogen, will be in the dead plant tissue and not removable
	double SheathNitrogenDecline;		//Z stepwise dead portion of sheath nitrogen, will be in the dead plant tissue and not removable
	double LeafNitrogenRelease;			//Z release nitrogen when a portion of leaf releases nitrogen (all the nitrogen is mobile and very valuable)
	double SheathNitrogenRelease;		//Z release nitrogen when a portion of sheath releases nitrogen
	double DeadLeafNitrogenMass;		//Z the dead portion of leaf nitrogen, will be in the dead plant tissue and not removable
	double DeadSheathNitrogenMass;		//Z the dead portion of sheath nitrogen, will be in the dead plant tissue and not removable

	//Z environemtal data
	double TmprLeaf;
	double PreDawnLWP;
	double N_effect;
	double water_effect_expand;
	double water_effect_senesc;
	double shade_effect;
	double Tmpr_effect;
	double slw_effect;

	//Z Tmpr dependent Leaf Expansion Parameters from K.Paff
	double LeafExpan_TmprMax;
	double LeafExpan_TmprMin;
	double LeafExpan_TmprOpt;

	//-----------------------------------------------------------------
	/*Z Leaf geometry and mass properties after considering "livingFrac",
	    this will make this plant and leaf become "representative plant" or "representative leaf" from the statistical prespective
	    In the maincode, leaf growth and shape is computated based on "one single leaf"
	    while we use "livingFrac" to adjust the "one single leaf" to be a "representative leaf"
	    we use "_Rep" to indicate "Representative".
	*/

	//Z why we adjust "livingfrac" based on leaf, a very deep/elementry level class?
	//  that is because our model design that leaf handles its own shape and growth STRICTLY within the leaf class

	double GreenLfArea_Rep;
	double GreenLfLength_Rep;
	double GreenLfWidth_Rep;
	double SenescentLfArea_Rep;
	double DropLfArea_Rep;
	double LeafArea_Rep;

	double DropLfMass_Rep;
	double DropSheathMass_Rep;
	double LeafMass_Rep;
	double SheathMass_Rep;
	double GreenLfMass_Rep;
	double GreenSheathMass_Rep;
	double LeafNitrogenMass_Rep;
	double SheathNitrogenMass_Rep;
	
	double LeafNitrogenRelease_Rep;
	double SheathNitrogenRelease_Rep;
	double DeadLeafNitrogenMass_Rep;
	double DeadSheathNitrogenMass_Rep;
	//Z for dead leaf/sheath nitrogen mass rep, we further partition that into
	//  "Yellow": N in the senecent portion but not dropped from the stem
	//  "Dropped" N in the dropped leaf, so already in the residue
	//  Note that for nitrogen: 
	//     Dead = Yellow + Dropped
	double YellowLeafNitrogenMass_Rep;
	double YellowSheathNitrogenMass_Rep;
	double DroppedLeafNitrogeMass_Rep;
	double DroppedSheathNitrogeMass_Rep;

	double ptnLfMassIncrease_Rep;
	double ptnLfNitrogenMassIncrease_Rep;

	double relative_growth;

};
#endif