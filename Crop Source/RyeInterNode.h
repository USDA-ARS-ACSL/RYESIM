#pragma once
#ifndef _RYE_INTERNODE_H_
#define _RYE_INTERNODE_H_
#include "weather.h"
#include "initinfo.h"
#include "RyeDevelopment.h"

class RyeInterNode
{
public:
//Z take current leaf rank and total leaf number to calculate potentialArea
	RyeInterNode(int rank, int order, bool mainstem, RyeDevelopment* dv, double livingFrac); 
	~RyeInterNode();
	void RyeInterNodeLengthUpdate(void);
	void RyeInterNodeMassUpdate(double biomassIncomeRate, double nitrogenIncomeRate);
	void IntrLengthEnlongation(void);
	void RyeInterNodeDeath();
	void InterNodeLivingFractorAdjustment();

	//Z stress functions
	//  tmpr
	//  canopy water potential (now in development.cpp)
	double TmprEffect(double Temperature);
//	double LWPeffect(double predawn_psi_bars, double threshold);

//Z IO living fraction, initialize, setting, and getting
	double get_IntrLivingFrac() { return livingFrac; }
	void set_IntrLivingFrac(double x)
	{
		if (!living) { livingFrac = 0.0; }
		else { livingFrac = x; }
	}

//Z Output Functions
	// This should be representative values
	double get_InterNodeLength() { return IntrLength_Rep; }
	double get_InterNodeMass() { return IntrMass_Rep; }
	double get_InterNodeLivingMass() { return IntrLivingMass_Rep; }
	double get_InterNodeNitrogenMass() { return IntrNitrogenMass_Rep; }
	double get_InterNodeNitrogenRelease() {return IntrNitrogenRelease_Rep;}
	double get_InterNodeNitrogenDeadMass() { return DeadIntrNitrogenMass_Rep; }

	//Z remove this output function
	// we are making the model based on mass incremental
	// double get_potentialInterNodeLengthIncrease() { return ptnIntrLengthIncrease_Rep; }
	double get_potentialInterNodeMassIncrease() { return ptnIntrMassIncrease_Rep; }
	double get_potentialInterNodeNitrogenMassIncrease() { return ptnIntrNitrogenMassIncrease_Rep; }

/*Z Input Function
    set three variables from parent tiller
    set this only one before the whole plant/tiller is update, but not after tiller update because there maybe new leaves appear
*/   
	void set_ElongationRank(int x) { rank_2_elongation = x; }

private:

	RyeDevelopment* develop;

	//Z Growing Stage
	int  rank;								//Z number of the internode at that tiller (mainstem) counted from bottom
	int  order;								//Z the order of tiller that this internode grows on
	bool mainstem;
	bool elongationStart;					//Z if elongation is started
	bool living;
	bool force_to_death_current_step;
	bool is_singleRidge;

	/*Z living fraction
		trace the living fraction changes
		for one internode, it is living or dead; for a representative plant, need to consider statistical effects
		therefore, an internode can be partially living.
		when living fraction changes, caused by tiller living fraction, those portion should be considered as dead
		therefore, computation of droping should be done and N should be recycled
	*/
	double livingFrac;
	double livingFrac_old;

	//Z environmental conditions
	double TmprIntrNode;
	double PreDawnLWP;
	double Cur_TTd;

	//Z Growing time/gdd measure
	double IntrAge;					//Z chronological ordinary time, day
	double physAge;					//Z gdd time (in reference to endGrowth and lifeSpan, days)

	//Z Internode Length (cm)
	double IntrLength;				//Z internode length cm
	double ptnIntrLength;			//Z internode length after one time step of growth cm
	double maxIntrLength;			//Z max internode length based on 
	double ptnIntrLengthIncrease;	//Z potnetial internode enlongation, cm
	double ptnIntrLengthDecrease;   //Z potential internode decrease due to death, cm
	double DeadIntrLength;			//Z dead internode length, cm, if the whole tiller is dead

	//Z Internode biomass (g)
	double stl = 0.015;		//Z internode specific weight g cm^-1
	double IntrMass;				//Z internode biomass in g
	double IntrMassIncrease;		//Z biomass assigned to this internode in g
	double ptnIntrMassIncrease;
	double DeadIntrMass;			//Z dead internode mass in g


	//Z Internode Nitrogen (mg)
	double IntrNitrogenContent;		//Z mass fraction in % (e.g., 3.0%)
	double IntrNitrogenMass;		//Z internode nitrogen mass mg
	double IntrNitrogenIncrease;	//Z nitrogen mass assigned to this internode in mg
	double ptnIntrNitrogenMassIncrease;
	double IntrNitrogenRelease;     //z nitrogen release due to internode decreasing or dying in mg

	double IntrNitrogenReleaseThreshold; //Z the min N mass content, below that, no N will be released even that organ is dead, assume to be 0.5%
	double IntrNitrogenMaxReleasePtge;	 //Z during plant organ dying, max percent of N releasing
	double IntrNitrogenDecline;			 //Z stepwise dead portion of nitrogen, in the dead plant tissue and not removable
	double DeadIntrNitrogenMass;		 //Z cumulated dead portion of nitrogen, in the dead plant tissue and not removable

	//Z Parent tiller information, which internode should elongate
	int rank_2_elongation;

	//-----------------------------------------------------------------
	/*Z Internode geometry (length) and mass properties after considering "livingFrac",
		this will make this plant and internode become "representative plant" or "representative leaf" from the statistical prespective
		In the maincode, internode growth and length is computated based on "one single leaf"
		while we use "livingFrac" to adjust the "one single leaf" to be a "representative leaf"
		we use "_Rep" to indicate "Representative".
	*/

	//Z why we adjust "livingfrac" based on internode, this very deep/elementry level class?
	//  that is because our model design that internode should handle its own shape and growth STRICTLY within the its own class

	//Z we do not count internode death without "livingfrac", 
	//  but if we use livingfrac, we have to define a new variable because tiller/internode could die any time, so we need to count those correctly
	double IntrDeadLength_Rep;
	double IntrLivingLength_Rep;
	double IntrLength_Rep;

	double IntrDeadMass_Rep;
	double IntrLivingMass_Rep;
	double IntrMass_Rep;

	double IntrNitrogenMass_Rep;
	double IntrNitrogenRelease_Rep;

	double DeadIntrNitrogenMass_Rep;

	double ptnIntrMassIncrease_Rep;
	double ptnIntrNitrogenMassIncrease_Rep;

	double relative_growth;
};

#endif _RYE_INTERNODE_H_
