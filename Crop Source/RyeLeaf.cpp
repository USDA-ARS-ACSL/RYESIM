#include "stdafx.h"
#include "RyeLeaf.h"
#include "weather.h"
#include "initinfo.h"
#include <cmath>
#include <algorithm>
#define DAYPERMINUTES 0.00069444444
#define MAXLEAFNUM 20
#define PHYLLOCHRON 106.0
#define MAX_N_PCT 3.5
#define MIN_N_PCT 0.35

using namespace std;

RyeLeaf::RyeLeaf(int n, int o, bool m, RyeDevelopment* dv, double LfLivingFrac)
{
	//Z ****** set default information or initialize variables *******
	// DT the first leaves on the mainstem that come out of the seed need to be intialized
	// differently as they have already started growing. We will assume their leaf lenght
	// is 30% of the maximum
	rank = n;
	order = o;
	mainstem = m;
	develop = dv;

	//Z temperature parameter for Tmpr Effects, from K.Paff
	LeafExpan_TmprMax = 16.0;
	LeafExpan_TmprMin = 0.0;
	LeafExpan_TmprOpt = 8.0;
	TmprLeaf = 25.0;
	PreDawnLWP = -0.05;
	N_effect = 1.0;
	water_effect_expand = 1.0;
	water_effect_senesc = 1.0;
	shade_effect = 1.0;
	Tmpr_effect = 1.0;
	slw_effect = 1.0;
	relative_growth = 1.0;

	//Z Growing Stage
	initiated = false;
	appeared = false;
	growing = false;
	mature = false;
	aging = false;
	dead = false;
	dropped = false;
	acceleration = false;
	force_to_death_current_step = false;
	is_singleRidge = false;
	livingFrac = LfLivingFrac;
	livingFrac_old = LfLivingFrac;
	livingFrac_ini = LfLivingFrac;

	//Z Leaf area recorder (cm or cm^2)
	// actual values; max (potential) values based on rank; potential (incremental) value for current step
	LeafLength = 0.0001;		//Z leaf length in cm
	LeafWidth = 0.0001;			//Z leaf width in cm
	GreenLfLength = 0.0001;		//Z leaf length when it is not dropped, if dropped, will be 0
	GreenLfWidth = 0.0001;		//Z leaf width when it is not dropped, if dropped, will be 0
	maxLfLength = 0.0001;		//Z max leaf length based on leaf rank and if it is at the mainstem
	maxLfWidth = 0.0001;		//Z max leaf width based on leaf rank and if it is at the mainstem
	ptnLfLength = 0.0;			//Z potential leaf length, the values after one step of potential growth
	ptnLfWidth = 0.0;			//Z potential leaf width, the values after one step of potential growth
	ptnLfLengthIncrease = 0.0;	//Z potential leaf length increase at the current time step
	ptnLfWidthIncrease = 0.0;	//Z potential leaf width increase at the current time step

	//Z shape factor "area = shape factor * length * width", dimensionless
	if (rank == 1 && mainstem == true) {
		A_LW = 0.83;		// Z: leaf area / leaf L*W for the first leaf in the mainstem
	}
	else {
		A_LW = 0.74;		// Z: leaf area / leaf L*W for other leaves
	}

	//Z ****** initialize the leaf after established from the tiller *******
		//Z initialize the leaf length, width based on the leaf rank
		//  max leaf lenght and width (cm) 
		//  for leaf rank >10, leaf size becomes constant

		//  assume leaf can be up to 30 cm long and 1.2 cm wide (based on some ag extension papers)
	if (rank <= 10) {
		double leaf_a_coef = exp(0.15 * (double)rank);
		maxLfLength = 11 * leaf_a_coef;
		maxLfWidth = 0.6 * leaf_a_coef;
	}
	else {
		double leaf_a_coef = exp(0.15 * 10.0);
		maxLfLength = 11 * leaf_a_coef;
		maxLfWidth = 0.6 * leaf_a_coef;
	}
	//Z max(reference) leaf area at that tiller and rank(cm^2)
	maxLeafArea = A_LW * maxLfLength * maxLfWidth;

	LeafArea = A_LW * LeafLength * LeafWidth;	//Z total leaf area changing over time, cm^2
	MaturedLfArea = 0.0;						//Z leaf area at the end of expansion, i.e., the actual max leaf area at one time (different from LeafArea)
	SenescentLfArea = 0.0;						//Z aged leaf area changing over time, cm^2
	DropLfArea = 0.0;							//Z drop leaf area, the leaf area at drop (dead) time point
	GreenLfArea = 0.0;							//Z green (active) leaf area, i.e., leaf area before aging, and part of the leaf area after aging
	ptnLeafArea = 0.0;							//Z potential leaf area, the values after one step of potential growth
	ptnLfAreaIncrease = 0.0;					//Z potential leaf area increase, the values after one step of potential growth
	ptnLfAreaDecrease = 0.0;					//Z potential leaf area increase, the values after one step of aging

	//Z Growing time/gdd measure
	LeafAge = 0.0;				// chronological ordinary time, day
	physAge = 0.0;				// gdd time (in reference to endGrowth and lifeSpan, days)
	GDD2mature = 0.0;			// gdd at mature, after fully expansion
	activeAge = 0.0;			// gdd since mature, the total functioning time period
	seneAge = 0.0;				// gdd since the end of activeAge, period for senescence
	Cur_TTd = 0.0;				// gdd current incremental (after photoperiod and vernalisation adjustment)
	seneFrac = 0.0;				// current seneFrac during aging period
	seneFrac_old = 0.0;			// old seneFrac during aging period

	/*Z leaf biomass + sheath biomass (g)
		dropped mass is still mass, so even if the leaf is dropped, it should be counted for leafmass
		therefore, if needed, LeafMass-DropLfMass = "attached leafmass"
		(for single leaf, either attached leafmass 0 or = leafmass, since leaf is either attached or dropped)
		(attached leafmass = green leaf mass + aged but not dropped leaf mass)
	*/
	seedMass = 0.0; // needed when initializing the two initial leaves
	slw_max = slw * 1.0;
	slw_run = slw;
	LeafMass = LeafArea * slw;
	SheathMass = 0.5 * LeafMass;
	GreenLfMass = LeafMass;
	GreenSheathMass = SheathMass;
	LeafBiomassIncrease = 0.0;
	SheathBiomassIncrease = 0.0;
	DropLfMass = 0.0;
	DropSheathMass = 0.0;

	//Z Leaf Sheath Nitrogen (mg): income
	//  assume 4.0% will be the max accessible value (at no stress)
	LeafNitrogenContent = MAX_N_PCT;								// leaf N % 
	SheathNitrogenContent = MAX_N_PCT;								// sheath N % 
	LeafNitrogenMass = LeafMass * LeafNitrogenContent * 10.0;		// leaf nitorgen mass (mg N), 10.0=/100.0*1000.0
	SheathNitrogenMass = SheathMass * SheathNitrogenContent * 10.0;	// sheath nitrogen mass
	LeafNitrogenIncrease = 0.0;			// leaf nitrogen mass increase after nitrogen assignment 
	SheathNitrogenIncrease = 0.0;		// sheath nitrogen mass increase after nitrogen assignment 
	
	//Z Leaf Sheath Nitrogen (mg): output portion, and N redistribution
	/*Z Reduce_Single : single leaf N output due to senecense and dropping, need to divide
	*   Decline : at that step, reduced N goes to residue and immobilized
	*   Release : at that step, remobilized and can reused for new plant organs
	*   DeadNitrogen: cumulated N mass in the senecense/dead plant organs
	*/
	LfNitrogenReleaseThreshold = 0.5;
	LfNitrogenMaxReleasePtge = 80.0;

	LeafNitrogenReduce_Single = 0.0;
	SheathNitrogenReduce_Single = 0.0;
	LfNitrogenDecline = 0.0;
	SheathNitrogenDecline = 0.0;
	LeafNitrogenRelease = 0.0;
	SheathNitrogenRelease = 0.0;
	DeadLeafNitrogenMass = 0.0;
	DeadSheathNitrogenMass = 0.0;

	//Z when tiller class create the new leaf object, its initiated and appeared
	initiated = true;
	appeared = true;

	//Z potential pool for carbon (g) and nitrogen (mg) that can be assigned to
	ptnLfMassIncrease = 0.0;
	ptnLfNitrogenMassIncrease = 0.0;

	//----------------------------------------------------------
	//Z initialize the representative value
	//  this changes "one single leaf" to "representative leaf"
	GreenLfArea_Rep = 0.0;
	GreenLfLength_Rep = 0.0;
	GreenLfWidth_Rep = 0.0;
	SenescentLfArea_Rep = 0.0;
	DropLfArea_Rep = 0.0;
	LeafArea_Rep = 0.0;

	DropLfMass_Rep = 0.0;
	DropSheathMass_Rep = 0.0;
	GreenLfMass_Rep = 0.0;
	GreenSheathMass_Rep = 0.0;
	LeafMass_Rep = 0.0;
	SheathMass_Rep = 0.0;
	LeafNitrogenMass_Rep = 0.0;
	SheathNitrogenMass_Rep = 0.0;

	LeafNitrogenRelease_Rep = 0.0;
	SheathNitrogenRelease_Rep = 0.0;
	DeadLeafNitrogenMass_Rep = 0.0;
	DeadSheathNitrogenMass_Rep = 0.0;
	YellowLeafNitrogenMass_Rep = 0.0;
	YellowSheathNitrogenMass_Rep = 0.0;
	DroppedLeafNitrogeMass_Rep = 0.0;
	DroppedSheathNitrogeMass_Rep = 0.0;

	ptnLfMassIncrease_Rep = 0.0;
	ptnLfNitrogenMassIncrease_Rep = 0.0;
}
// New constructor implementation
RyeLeaf::RyeLeaf(int n, int o, bool m, RyeDevelopment* dv, double LfLivingFrac, double seedMass)
	: RyeLeaf(n, o, m, dv, LfLivingFrac)

{
	// Additional initialization specific to this constructor
	this->seedMass = seedMass;
	double shootmass = 0.6 * this->seedMass;
	if (rank == 1 && mainstem == true) {
		LeafLength = maxLfLength * 0.3;
		LeafWidth = maxLfWidth * 0.3;
		LeafArea = A_LW * LeafLength * LeafWidth;
		// Recalculate dependent values

		LeafMass = 0.8 * shootmass;
		SheathMass = 0.5 * (shootmass-LeafMass);
		GreenSheathMass = SheathMass;
		GreenLfArea = LeafArea;
		GreenLfMass = LeafMass;
	}
	// Call existing constructor logic or duplicate initialization
	// ... existing initialization code ...
}

//Z leaf death does not means the objective is deleted
//  drop leaf mass variable still hold values
//  no memory allocation need to be released at this step
RyeLeaf::~RyeLeaf() {}

//Z update for the area
//  but under leaf senescence, N release based on the senescence fraction is 
void RyeLeaf::RyeLeafAreaUpdate()
{
	TmprLeaf = develop->get_Tcur();						//Z ambient temperature
	Cur_TTd = develop->get_TTd();						//Z ggd incremental at this (hourly) time step, unit in "degree day" 
	PreDawnLWP = develop->get_PredawnLWP();				//Z predawn leaf water potential MPa
	LeafAge += develop->get_TimeStep() * DAYPERMINUTES;	//Z this is actually 1/24 (one hour), but TimeStep is in minute, therefore, convert to the number of day
	physAge += Cur_TTd;

	//Z computation of stress factors
	//Z nitrogen argument (%mass)
	//  give a mini growth rate effects 10%
	double CriticalNitrogen;
	CriticalNitrogen = __max(MIN_N_PCT, LeafNitrogenContent);
	N_effect = __max(0.0, (2.0 / (1.0 + exp(-2.9 * (CriticalNitrogen - MIN_N_PCT))) - 1.0));
	N_effect = __min(1.0, N_effect);
	N_effect = __max(0.1, N_effect);
	//N_effect = 1.0;
	//cout << "RyeLeaf: " << N_effect << endl;

//	if (!develop->is_singleRidge()) {
//		N_effect = __max(0.8, sqrt(N_effect));
//	}
//	if (develop->is_singleRidge()) {
//		N_effect = __max(0.6, sqrt(N_effect));
//	}
	if (!develop->is_startEnlongation()) {
		N_effect = __max(0.3, sqrt(N_effect));
		N_effect = 0.3 + 0.7 * N_effect;
	}
	else
	{
//		N_effect = __max(0.6, sqrt(N_effect));
	}
	
	if (develop->is_singleRidge() && (!is_singleRidge)) {
	    if (develop->get_ColdTime() < 16.0) {   //30.0
			relative_growth = 1.0;  //1.8    EJ
		}
		else if (develop->get_ColdTime() < 25.0) {   //30.0
			relative_growth = 2.3;  //EJ
		}
		else if (develop->get_ColdTime() < 50.0) {   //EJ never reach to this point for AL sites
			relative_growth = 2.3; // 1.8;
		}
		else {
			relative_growth = 1.0;
		}
		is_singleRidge = true;
	}
	if (!develop->is_singleRidge()) {
		relative_growth = 1.4;
	}
	

//	if (develop->is_accel()) {
//		N_effect = __max(0.7, N_effect);
//	}
//	if (develop->is_startEnlongation()) {
//		N_effect = __max(0.5, N_effect);
//	}

	//Z water potential argument (%mass)
	//  1. water stress for expand
	water_effect_expand = develop->get_leafWaterEffectExpand();
	//  2. water stress for senesc
	water_effect_senesc = develop->get_leafWaterEffectSenesc();

	//Z shaded effect (for the whole plant)
	shade_effect = develop->get_shadeEffect();

	//Z tmpr adjustment
	Tmpr_effect = TmprEffect(TmprLeaf);

	//Z slw limiter
	slw_run = LeafMass / LeafArea;
	if (slw_run < slw_min) { slw_effect = 0.001; }			// slw is very very low, nearly no grow
	else if (slw_run < 1.5 * slw_min) { slw_effect = 0.5; } // warning, slw should not be that low. This is my guess while papers saying ~2 should be the warning level
	else { slw_effect = 1.0; }

	// ----------------------
	//Z the next two functions are for "one single leaf (not fractional)" morphology computation,
	//  but finally ajust the results based on the living fraction

	//Z leaf expansion is only for leaf area
	LeafAreaExpand();
	//Z leaf senescence is for leaf area and nitrogen release
	//  leaf biomass carbon is structure carbon and will not leave
	LeafAreaSenescence();

	//Z update leaf nitrogen content 
	//Z during senescence, leafMass is always the mass at maturity and not changes, so use (1.0 - seneFrac) to adjust the actual leaf mass after part of the leaf is yellow.
	if (seneFrac < 1.0)	{
		LeafNitrogenContent = __max(__min(0.10 * LeafNitrogenMass / LeafMass / (1.0 - seneFrac), MAX_N_PCT), 0.0);
	}
	else {
		LeafNitrogenContent = MAX_N_PCT;
	}

	// ----------------------
	//Z the next function adjust the living fractions
	//  if livingFrac<livingFrac_old based on stress, then put "livingFrac_old-livingfrac" to die
	//  in principle, livingfrac<=livingfrac_old
	LeafLivingFractorAdjustment();

}

//Z leaf expansion rate, only for leaf area
void RyeLeaf::LeafAreaExpand()
{
	/*Z this will be for single leaf geometry,
	so be careful that we do not use "living fraction" for adjustments
	we save living fraction adjustments in a separate function
	*/

	//Z reset for each time step
	ptnLfLengthIncrease = 0.0;
	ptnLfWidthIncrease = 0.0;
	ptnLfLength = 0.0;
	ptnLfWidth = 0.0;
	ptnLeafArea = 0.0;
	ptnLfAreaIncrease = 0.0;
	ptnLfMassIncrease = 0.0;
	ptnLfNitrogenMassIncrease = 0.0;

	//Z leaf and sheath expansion
	//     until the max leaf lenght is reached, and then set the leaf status to mature
	if (appeared && !mature)
	{
		growing = true;
		//N_effect = 1.0;
		double limit_factor = __min(__min(__min(N_effect, water_effect_expand), __min(shade_effect, Tmpr_effect)), slw_effect);
		limit_factor = __max(sqrt(limit_factor), 0.1);
		double grll = maxLfLength / PHYLLOCHRON * limit_factor * relative_growth;
		double grlw = maxLfWidth / PHYLLOCHRON * limit_factor * relative_growth;

		ptnLfLengthIncrease = grll * Cur_TTd;
		ptnLfWidthIncrease = grlw * Cur_TTd;
		ptnLfLength = LeafLength + ptnLfLengthIncrease;
		ptnLfWidth = LeafWidth + ptnLfWidthIncrease;

		//Z leaf area govened maturity
		if (ptnLfLength >= maxLfLength && ptnLfWidth >= maxLfWidth)
		{
			// Z compute the leaf area 
			ptnLfLengthIncrease = maxLfLength - LeafLength;
			ptnLfWidthIncrease = maxLfWidth - LeafWidth;
			ptnLeafArea = maxLeafArea;
			ptnLfAreaIncrease = maxLeafArea - LeafArea;

			LeafLength = LeafLength + ptnLfLengthIncrease;
			LeafWidth = LeafWidth + ptnLfWidthIncrease;
			LeafArea = ptnLeafArea;
			MaturedLfArea = LeafArea;

			mature = true;
			GDD2mature = physAge;
			growing = false;

			GreenLfLength = LeafLength;
			GreenLfWidth = LeafWidth;
		}
		//Z leaf physical age governed maturity
		//  need a time-based criterion here, otherwise leaf will grow unlessly
		else if (physAge >= 1.0 * PHYLLOCHRON)
		{
			ptnLeafArea = A_LW * ptnLfLength * ptnLfWidth;
			ptnLfAreaIncrease = ptnLeafArea - LeafArea;
			LeafLength = ptnLfLength;
			LeafWidth = ptnLfWidth;
			LeafArea = ptnLeafArea;
			MaturedLfArea = LeafArea;

			mature = true;
			GDD2mature = physAge;
			growing = false;

			GreenLfLength = LeafLength;
			GreenLfWidth = LeafWidth;
		}
		//Z leaf is not matured, grow normally
		else
		{
			ptnLeafArea = A_LW * ptnLfLength * ptnLfWidth;
			ptnLfAreaIncrease = ptnLeafArea - LeafArea;
			LeafLength = ptnLfLength;
			LeafWidth = ptnLfWidth;
			LeafArea = ptnLeafArea;

			growing = true;

			GreenLfLength = LeafLength;
			GreenLfWidth = LeafWidth;
		}
	}

	if (!aging)
	{
		//Z potential mass increase leaf (g)
		//  leaf not necessarily grow at the stage
		//  but can still fit biomass until the slw number is reached
		ptnLfMassIncrease = __max(LeafArea * slw_max - LeafMass, 0.0);
		//Z potential nitrogen increase leaf (mg)
		ptnLfNitrogenMassIncrease = __max(LeafArea * slw_max * MAX_N_PCT * 10.0 - LeafNitrogenMass, 0.0);
	}
	else 
	{
		ptnLfMassIncrease = 0.0;
		ptnLfNitrogenMassIncrease = 0.0;
	}

	GreenLfArea = __max(0.0, LeafArea - SenescentLfArea);

	return;
}

//Z leaf senescence, this is the very sloppy way from ShootGro
//  this function include the N release due to senescence
//  technically, leaf expanion and senescence will not occur at the same time, separated by a "long activating photosync" period
void RyeLeaf::LeafAreaSenescence()
{
	/*Z this will be for single leaf geometry,
		so be careful that we do not use "living fraction" for adjustments
		we save living fraction adjustments in a separated function
	*/

	/*Z if the leaf is force to death, DO NOT run "LeafAreaSenescence"
		if wrongly run this function, N release will be zeroed again WRONGLY
		"Senescence" and "Force to die" are "Mutually Exclusive"
		that is the reason the followling if statement is in the "LeafAreaSenescence()" function
		but NOT in the "LeafAreaExpand()" function
	*/
	if (force_to_death_current_step)
	{
		//Z leaf is force to death at this step.
		//  bypass the normal "LeafAreaSenescence()" function
		force_to_death_current_step = false;
		return;
	}

	//Z reset for each time step
	LeafNitrogenRelease = 0.0;
	SheathNitrogenRelease = 0.0;
	ptnLfAreaDecrease = 0.0;
	LfNitrogenDecline = 0.0;
	SheathNitrogenDecline = 0.0;

	double LfActiveLivingPhyll = 6.5;
	LfActiveLivingPhyll = __max(LfActiveLivingPhyll * __min(N_effect, water_effect_senesc), 2.5);

	//Z leaf senescence main code
	if (!mature && !aging && !dead)
	{
		//leaf not ready to senesce
	}
	else if (mature && !aging && !dead)
	{
		activeAge = physAge - GDD2mature;
		if (activeAge >= LfActiveLivingPhyll * PHYLLOCHRON)
		{
			aging = true;
		}
	}
	/*Z senecense will be finished within one PHYLLOCHRON
		senecense will NOT change leaf area, leaf area include green leaf, senecense and drop leaf
		senecent leaf is a cumulative value, it increases with respect to time and gradually reach leaf area
		when "senecent leaf = leaf area", leaf drop,
		then for one single leaf, leaf area = senecent leaf area = dropped leaf area
	*/
	else if (aging && !dead)
	{
		seneAge = physAge - activeAge - GDD2mature;
		seneFrac = seneAge / PHYLLOCHRON;
		if (seneFrac >= 1.0)
		{
			//Z leaf area senscene fraction
			ptnLfAreaDecrease = MaturedLfArea - SenescentLfArea;
			SenescentLfArea = MaturedLfArea;
			seneFrac = 1.0;

			//Z leaf totally dead, all N need to be reduced and redistributed
			LeafNitrogenReduce_Single = LeafNitrogenMass;
			SheathNitrogenReduce_Single = SheathNitrogenMass;
			LeafNitrogenMass = 0.0;
			SheathNitrogenMass = 0.0;

			//Z leaf and sheath mass assignment
			GreenLfMass = 0.0;
			GreenSheathMass = 0.0;

			//Z determine N release in mg
			//  all N will be released or decline (to dead portion), this leaf is dead
			if (LeafNitrogenContent <= LfNitrogenReleaseThreshold)
			{
				LeafNitrogenRelease = 0.0;
				SheathNitrogenRelease = 0.0;
				LfNitrogenDecline = LeafNitrogenReduce_Single;
				SheathNitrogenDecline = SheathNitrogenReduce_Single;

				DeadLeafNitrogenMass += LfNitrogenDecline;
				DeadSheathNitrogenMass += SheathNitrogenDecline;

			}
			else
			{
				//Z use "LfNitrogenMaxReleasePtge" to define a temporary fraction for 
				//  N decline -> stay in the dead/yellow organs or residue
				//  N release -> reused by other organs
				double NitrogenReleaseFrac = (LeafNitrogenContent - LfNitrogenReleaseThreshold) / LeafNitrogenContent;
				NitrogenReleaseFrac = __max(__min(NitrogenReleaseFrac, LfNitrogenMaxReleasePtge / 100.0), 0.0);
				double NitrogenDeclineFrac = 1.0 - NitrogenReleaseFrac;

				LeafNitrogenRelease = LeafNitrogenReduce_Single * NitrogenReleaseFrac;
				SheathNitrogenRelease = SheathNitrogenReduce_Single * NitrogenReleaseFrac;
				LfNitrogenDecline = LeafNitrogenReduce_Single * NitrogenDeclineFrac;
				SheathNitrogenDecline = SheathNitrogenReduce_Single * NitrogenDeclineFrac;

				DeadLeafNitrogenMass += LfNitrogenDecline;
				DeadSheathNitrogenMass += SheathNitrogenDecline;
				
			}

			seneFrac_old = 1.0;

			dead = true;

			GreenLfLength = 0.0;
			GreenLfWidth = 0.0;
		}
		else
		{
			//Z leaf area senscene fraction
			//  compute the senescent fraction comparing to the green leaf area last time
			double relative_sene_fraction = __min(1.0, (seneFrac - seneFrac_old) / (1.0 - seneFrac_old));
			ptnLfAreaDecrease = __max(0.0, relative_sene_fraction * (MaturedLfArea - SenescentLfArea));
			SenescentLfArea = SenescentLfArea + ptnLfAreaDecrease;

			//Z determine N release in mg
			//Z leaf is still living
			//  first compute how much N will leaf the living part of leaf,
			LeafNitrogenReduce_Single = relative_sene_fraction * LeafNitrogenMass;
			SheathNitrogenReduce_Single = relative_sene_fraction * SheathNitrogenMass;
			LeafNitrogenMass = LeafNitrogenMass - LeafNitrogenReduce_Single;
			SheathNitrogenMass = SheathNitrogenMass - SheathNitrogenReduce_Single;

			//  second redistribute those N
			if (LeafNitrogenContent <= LfNitrogenReleaseThreshold)
			{
				LeafNitrogenRelease = 0.0;
				SheathNitrogenRelease = 0.0;
				LfNitrogenDecline = LeafNitrogenReduce_Single;
				SheathNitrogenDecline = SheathNitrogenReduce_Single;

				DeadLeafNitrogenMass += LfNitrogenDecline;
				DeadSheathNitrogenMass += SheathNitrogenDecline;
			}
			else
			{
				//Z use "LfNitrogenMaxReleasePtge" to define a temporary fraction for 
				//  N decline -> stay in the dead/yellow organs or residue
				//  N release -> reused by other organs
				double NitrogenReleaseFrac = (LeafNitrogenContent - LfNitrogenReleaseThreshold) / LeafNitrogenContent;
				NitrogenReleaseFrac = __max(__min(NitrogenReleaseFrac, LfNitrogenMaxReleasePtge / 100.0), 0.0);
				double NitrogenDeclineFrac = 1.0 - NitrogenReleaseFrac;

				LeafNitrogenRelease = LeafNitrogenReduce_Single * NitrogenReleaseFrac;
				SheathNitrogenRelease = SheathNitrogenReduce_Single * NitrogenReleaseFrac;
				LfNitrogenDecline = LeafNitrogenReduce_Single * NitrogenDeclineFrac;
				SheathNitrogenDecline = SheathNitrogenReduce_Single * NitrogenDeclineFrac;

				DeadLeafNitrogenMass += LfNitrogenDecline;
				DeadSheathNitrogenMass += SheathNitrogenDecline;
			}

			dead = false;

			GreenLfLength = (1.0 - relative_sene_fraction) * GreenLfLength;
			GreenLfWidth = (1.0 - relative_sene_fraction) * GreenLfWidth;

			//Z leaf and sheath mass assignment
			GreenLfMass = (1.0 - relative_sene_fraction) * GreenLfMass;
			GreenSheathMass = (1.0 - relative_sene_fraction) * GreenSheathMass;

			seneFrac_old = seneFrac;
		}
	}
	else if ((dead && (!dropped)) && physAge >= GDD2mature)
	{
		//Z after dead, allow one more time step to be dropped
		dropped = true;
		DropLfMass = LeafMass;
		DropSheathMass = SheathMass;
		DropLfArea = SenescentLfArea;

		GreenLfLength = 0.0;
		GreenLfWidth = 0.0;
		GreenLfMass = 0.0;
		GreenSheathMass = 0.0;
	}

	GreenLfArea = __max(0.0, LeafArea - SenescentLfArea);

	return;
}

/*Z living factor adjustment
	after the "single leaf" morphology computation (geometry)
	use the living fraction to determine the "statistical/expectation" of the leaf morphology
*/
//Z we use "Rep" as "representative" or "statistically representative"
//  used for field scale mass or geometry computations
void RyeLeaf::LeafLivingFractorAdjustment()
{
	/*Z living fraction or living fraction old may NOT be necessary started at 1
		that is because if the tiller's living fraction is from, e.g. 0.5,
		then even the leaf is fully expanded, its living fraction is at most 0.5
		*/

	ptnLfMassIncrease_Rep = 0.0;
	ptnLfNitrogenMassIncrease_Rep = 0.0;
	LeafNitrogenRelease_Rep = 0.0;
	SheathNitrogenRelease_Rep = 0.0;

	if (livingFrac < livingFrac_old)
	{
		//Z livingfrac must be smaller than livingfrac old
		double livingDiff = livingFrac_old - livingFrac;

		//-------------------------
		//Z adjust living leaf geometry, mass and potential mass required for leaf growth 
		//  green part
		GreenLfArea_Rep = GreenLfArea * livingFrac;
		GreenLfLength_Rep = GreenLfLength * livingFrac;
		GreenLfWidth_Rep = GreenLfWidth * livingFrac;
		GreenLfMass_Rep = GreenLfMass * livingFrac;
		GreenSheathMass_Rep = GreenSheathMass * livingFrac;
		//  senecent part
		/*  senecent leaf portion has to parts
			1. senecent of living leaf
			2. dead/dropped whole leaf
			This is a cumulative value so must be initialized to 0
		*/
		SenescentLfArea_Rep += (ptnLfAreaDecrease * livingFrac + __max(MaturedLfArea - (SenescentLfArea - ptnLfAreaDecrease), 0.0) * livingDiff);
		//  drop part, for the equation below
		/*  1. if tiller dead, leaf will be dropped, so "DropLfArea_Rep" is a cumulative value
			2. current living leaf can drop, shown in the first term
			3. leafArea is everything, i.e., green, senecent and dead/dropped (exclusive), so no matter how those leaves shown, once tiller dead, it goes to the second term
			4. "exclusive" in 3 means if DropLfArea!=0, it must = LeafArea = SenescentLfArea, and livingFrac = livingDiff = 0
		*/
		DropLfArea_Rep += (DropLfArea * livingFrac + LeafArea * livingDiff);

		//  total leaf area, SenescentLfArea includes but not limit to drop leaf area
		LeafArea_Rep = GreenLfArea_Rep + SenescentLfArea_Rep;

		//-------------------------
		//Z adjust leaf/sheath biomass
		//  drop part, similar to the "DropLfArea", this should be a cumulative value
		//  note that usually (until the last step), DropLfMass should be always 0
		DropLfMass_Rep += (DropLfMass * livingFrac + LeafMass * livingDiff);
		DropSheathMass_Rep += (DropSheathMass * livingFrac + SheathMass * livingDiff);

		//  total leaf/sheath mass
		/*	1. if leaf is dropped, then rely on the second term
			2. if leaf is living, "DropLfMass" without caring the "livingfrac" will be 0, then both term shows
			3. similar to the sheath
		*/
		LeafMass_Rep = (LeafMass - DropLfMass) * livingFrac + DropLfMass_Rep;
		SheathMass_Rep = (SheathMass - DropSheathMass) * livingFrac + DropSheathMass_Rep;

		//-------------------------
		/*Z adjust leaf/sheath nitrogen mass and their contribution to nitrogen releasing
			1. mass computation is simple, since we assume only dead portion totally recycle nitrogen
		*/
		LeafNitrogenMass_Rep = LeafNitrogenMass * livingFrac;
		SheathNitrogenMass_Rep = SheathNitrogenMass * livingFrac;

		//-------------------------
		//Z adjust leaf/sheath nitrogen mass
		//  1. note that for the N decline to dead tissues due to livingfraction changes, still need to redistribute the N
		//  2. for releasing, the living portion contribute instantaneous releasing, will the dropped one contribute the total nitrogen
		//  3. NitrogenDecline is just a stepwise variable, and should eventually put into the dead leaf or dead sheath N portion

		DeadLeafNitrogenMass_Rep += LfNitrogenDecline * livingFrac;
		DeadSheathNitrogenMass_Rep += SheathNitrogenDecline * livingFrac;
		LeafNitrogenRelease_Rep = LeafNitrogenRelease * livingFrac;
		SheathNitrogenRelease_Rep = SheathNitrogenRelease * livingFrac;

		//  need to redistribute those N due to livingfraction changes
		//Z should be (LeafNitrogenMass+LfNitrogenDecline+LeafNitrogenRelease), 
		//  and those additional terms "LfNitrogenDecline+LeafNitrogenRelease=LeafNitrogenReduce_Single" can be higher order infinitesimal
		double LeafNitrogenFree_temp = (LeafNitrogenMass + LeafNitrogenReduce_Single) * livingDiff;
		double SheathNitrogenFree_temp = (SheathNitrogenMass + SheathNitrogenReduce_Single) * livingDiff;
		if (LeafNitrogenContent <= LfNitrogenReleaseThreshold)
		{
			// This part of decline N will be the dropped,
			// because for the change of "livingfraction", the portion belong to the differences means "dead" tiller 
			// so leaf on dead tiller means "dropped leaf"
			DeadLeafNitrogenMass_Rep += LeafNitrogenFree_temp;
			DeadSheathNitrogenMass_Rep += SheathNitrogenFree_temp;
			DroppedLeafNitrogeMass_Rep += LeafNitrogenFree_temp;
			DroppedSheathNitrogeMass_Rep += SheathNitrogenFree_temp;
		}
		else
		{
			//Z use "LfNitrogenMaxReleasePtge" to define a temporary fraction for 
				//  N decline -> stay in the dead/yellow organs or residue
				//  N release -> reused by other organs
			double NitrogenReleaseFrac = (LeafNitrogenContent - LfNitrogenReleaseThreshold) / LeafNitrogenContent;
			NitrogenReleaseFrac = __max(__min(NitrogenReleaseFrac, LfNitrogenMaxReleasePtge / 100.0), 0.0);
			double NitrogenDeclineFrac = 1.0 - NitrogenReleaseFrac;

			LeafNitrogenRelease_Rep += LeafNitrogenFree_temp * NitrogenReleaseFrac;
			SheathNitrogenRelease_Rep += SheathNitrogenFree_temp * NitrogenReleaseFrac;

			// This part of decline N will be the dropped,
			// because for the change of "livingfraction", the portion belong to the differences means "dead" tiller 
			// so leaf on dead tiller means "dropped leaf"
			double aaaa = LeafNitrogenFree_temp * NitrogenDeclineFrac;
			double bbbb = SheathNitrogenFree_temp * NitrogenDeclineFrac;
			DeadLeafNitrogenMass_Rep += aaaa;
			DeadSheathNitrogenMass_Rep += bbbb;
			DroppedLeafNitrogeMass_Rep += aaaa;
			DroppedSheathNitrogeMass_Rep += bbbb;

		}

		//Z Nitorgen in Yellow Leaf, dead (senescent and non-mobile anymore) but not dropped, should be the difference
		//  Yellow = Dead - Dropped
		YellowLeafNitrogenMass_Rep = DeadLeafNitrogenMass_Rep - DroppedLeafNitrogeMass_Rep;
		YellowSheathNitrogenMass_Rep = DeadSheathNitrogenMass_Rep - DroppedSheathNitrogeMass_Rep;

		//-------------------------
		//Z adjust potential leaf growth
		//  only living portion can grow, so this part is very simple
		//ptnLfAreaIncrease_Rep = ptnLfAreaIncrease * livingFrac;
		ptnLfMassIncrease_Rep = ptnLfMassIncrease * livingFrac;
		ptnLfNitrogenMassIncrease_Rep = ptnLfNitrogenMassIncrease * livingFrac;

		if (dropped)
		{
			livingFrac_old = 0.0;
			livingFrac = 0.0;

			GreenLfArea = 0.0;
			GreenLfLength = 0.0;
			GreenLfWidth = 0.0;
			GreenLfMass = 0.0;
			GreenSheathMass = 0.0;
			LeafNitrogenMass = 0.0;
			SheathNitrogenMass = 0.0;

			ptnLfMassIncrease_Rep = 0.0;
			ptnLfNitrogenMassIncrease_Rep = 0.0;

			//Z this representative leaf is totally dead, 
			YellowLeafNitrogenMass_Rep = 0.0;
			YellowSheathNitrogenMass_Rep = 0.0;
			DroppedLeafNitrogeMass_Rep = DeadLeafNitrogenMass_Rep;
			DroppedSheathNitrogeMass_Rep = DeadSheathNitrogenMass_Rep;

		}
		//Z finally update the living fraction numbers for this leaf
		livingFrac_old = livingFrac;
	}
	else
	{
		//Z basically do nothing here, but still need to count for the "Rep" values
		//  1. if "livingFrac=livingFrac_old = 0.0", protection path, i.e., ensure all the data is kept or zeroed
		//  2. if "livingFrac=livingFrac_old != 0.0", count for the "Rep" values

		//Z when the "tiller" is killed
		//  the last run of the leaf update function will enter this if statement
		//  and the "mass exchange variables" such as potential N demand or N release, will be set to 0
		if ((livingFrac == 0.0) && (livingFrac_old == 0.0))
		{
			//Z maybe (but not likely) the tiller is still living, but this leaf is dead and dropped anyway
			GreenLfArea_Rep = 0.0;
			GreenLfLength_Rep = 0.0;
			GreenLfWidth_Rep = 0.0;
			GreenLfMass_Rep = 0.0;
			GreenSheathMass_Rep = 0.0;

			//  senecent part
			//  senecent leaf portion is a cumulative quantity, so just stop the update
			//SenescentLfArea_Rep += (ptnLfAreaDecrease * livingFrac + __max(LeafArea + ptnLfAreaDecrease - SenescentLfArea, 0.0) * livingDiff);;

			//  drop part
			//  drop part is a cumulative quantity, so just stop the update
			//DropLfArea_Rep += (DropLfArea * livingFrac + LeafArea * livingDiff);

			//  total leaf area = SenescentLfArea
			LeafArea_Rep = SenescentLfArea_Rep;

			//-------------------------
			//Z adjust leaf/sheath biomass
			//  drop part, similar to the "DropLfArea", this should be a cumulative value
			//  therefore, just need to stop updating
			//DropLfMass_Rep += DropLfMass * livingFrac + LeafMass * livingDiff;
			//DropSheathMass_Rep += DropSheathMass * livingFrac + SheathMass * livingDiff;

			//  total leaf/sheath mass
			//	leaf is dead, so leaf/sheath mass = dropped leaf/sheath mass
			LeafMass_Rep = DropLfMass_Rep;
			SheathMass_Rep = DropSheathMass_Rep;

			//-------------------------
			//Z adjust leaf/sheath nitrogen mass and their contribution to nitrogen releasing
			//  leaf dead, all nitrogen is recycled
			LeafNitrogenMass_Rep = 0.0;
			SheathNitrogenMass_Rep = 0.0;
			LeafNitrogenRelease_Rep = 0.0;
			SheathNitrogenRelease_Rep = 0.0;

			//-------------------------
			//Z adjust potential leaf growth
			//  no leaf growth if leaf is dead
			//ptnLfAreaIncrease_Rep = 0.0;
			ptnLfMassIncrease_Rep = 0.0;
			ptnLfNitrogenMassIncrease_Rep = 0.0;

			//Z this representative leaf is totally dead, 
			YellowLeafNitrogenMass_Rep = 0.0;
			YellowSheathNitrogenMass_Rep = 0.0;
			DroppedLeafNitrogeMass_Rep = DeadLeafNitrogenMass_Rep;
			DroppedSheathNitrogeMass_Rep = DeadSheathNitrogenMass_Rep;

			//make sure the green or N part should be always 0 at this stage
			GreenLfArea = 0.0;
			GreenLfLength = 0.0;
			GreenLfWidth = 0.0;
			GreenLfMass = 0.0;
			GreenSheathMass = 0.0;
			LeafNitrogenMass = 0.0;
			SheathNitrogenMass = 0.0;

		}
		else
		{
			//Z Tiller induce no leaf death, but leaf can still die and drop itself
			//-------------------------
			//Z adjust living leaf geometry, mass and potential mass required for leaf growth 
			//  green part
			GreenLfArea_Rep = GreenLfArea * livingFrac;
			GreenLfLength_Rep = GreenLfLength * livingFrac;
			GreenLfWidth_Rep = GreenLfWidth * livingFrac;
			GreenLfMass_Rep = GreenLfMass * livingFrac;
			GreenSheathMass_Rep = GreenSheathMass * livingFrac;

			//  senecent part
			//  senecent leaf portion has to parts
			//  1. senecent of living leaf
			//  2. NO dead/dropped whole leaf DUE TO TILLER DEATH
			SenescentLfArea_Rep += (ptnLfAreaDecrease * livingFrac);

			//  drop part, for the equation below
			//  1. if tiller dead, leaf will be dropped, so "DropLfArea_Rep" is a cumulative value
			//  2. current living leaf can drop, shown in the first term
			//  3. NO dead/dropped whole leaf DUE TO TILLER DEATH
			//  4. note that both "DropLfArea" and "livingFrac" can be nonzero for ONLY ONE STEP
			DropLfArea_Rep += (DropLfArea * livingFrac);

			//  total leaf area, SenescentLfArea includes but not limit to drop leaf area
			LeafArea_Rep = GreenLfArea_Rep + SenescentLfArea_Rep;

			//-------------------------
			//Z adjust leaf/sheath biomass
			//  drop part, similar to the "DropLfArea", this should be a cumulative value
			//  note that both "DropLf/SheathMass" and "livingFrac" can be nonzero for ONLY ONE STEP
			DropLfMass_Rep += (DropLfMass * livingFrac);
			DropSheathMass_Rep += (DropSheathMass * livingFrac);

			//  total leaf/sheath mass
			//	1. if leaf is dropped, then rely on the second term
			//  2. if leaf is living, "DropLfMass" without caring the "livingfrac" will be 0, then both term shows
			//  3. similar to the sheath
			LeafMass_Rep = (LeafMass - DropLfMass) * livingFrac + DropLfMass_Rep;
			SheathMass_Rep = (SheathMass - DropSheathMass) * livingFrac + DropSheathMass_Rep;

			//-------------------------
			//Z adjust leaf/sheath nitrogen mass and their contribution to nitrogen releasing
			//  1. mass computation is simple, since we assume only dead portion totally recycle nitrogen
			//  2. for releasing, the living portion contribute instantaneous releasing, will the dropped one contribute the total nitrogen
			//  3. for decline/N in dead tissue
			//  4. note that "___NitrogenRelease" and "___NitrogenDecline" are stepwise variables
			LeafNitrogenMass_Rep = LeafNitrogenMass * livingFrac;
			SheathNitrogenMass_Rep = SheathNitrogenMass * livingFrac;
			LeafNitrogenRelease_Rep = LeafNitrogenRelease * livingFrac;
			SheathNitrogenRelease_Rep = SheathNitrogenRelease * livingFrac;
			DeadLeafNitrogenMass_Rep += LfNitrogenDecline * livingFrac;
			DeadSheathNitrogenMass_Rep += SheathNitrogenDecline * livingFrac;

			//Z no leaf drop in this stage, from representative prespective
			//  Thus, there is no change in the "DroppedLeafNitrogeMass_Rep" and "DroppedSheathNitrogeMass_Rep"
			//Z Nitorgen in Yellow Leaf, dead (senescent and non-mobile anymore) but not dropped, should be the difference
		    //  Yellow = Dead - Dropped
			YellowLeafNitrogenMass_Rep = DeadLeafNitrogenMass_Rep - DroppedLeafNitrogeMass_Rep;
			YellowSheathNitrogenMass_Rep = DeadSheathNitrogenMass_Rep - DroppedSheathNitrogeMass_Rep;

			//-------------------------
			//Z adjust potential leaf growth
			//  only living portion can grow, so this part is very simple
			//ptnLfAreaIncrease_Rep = ptnLfAreaIncrease * livingFrac;
			ptnLfMassIncrease_Rep = ptnLfMassIncrease * livingFrac;
			ptnLfNitrogenMassIncrease_Rep = ptnLfNitrogenMassIncrease * livingFrac;

			if (dropped)
			{
				livingFrac_old = 0.0;
				livingFrac = 0.0;

				GreenLfArea = 0.0;
				GreenLfLength = 0.0;
				GreenLfWidth = 0.0;
				GreenLfMass = 0.0;
				GreenSheathMass = 0.0;
				LeafNitrogenMass = 0.0;
				SheathNitrogenMass = 0.0;

				ptnLfMassIncrease_Rep = 0.0;
				ptnLfNitrogenMassIncrease_Rep = 0.0;

				//Z this representative leaf is totally dead, 
				YellowLeafNitrogenMass_Rep = 0.0;
				YellowSheathNitrogenMass_Rep = 0.0;
				DroppedLeafNitrogeMass_Rep = DeadLeafNitrogenMass_Rep;
				DroppedSheathNitrogeMass_Rep = DeadSheathNitrogenMass_Rep;
			}
			//Z finally update the living fraction numbers for this leaf
			livingFrac_old = livingFrac;
		}
	}
}

//Z leaf mass partition
//  should include leaf and sheath
//  the input parameters are based on the "representative" values of the leaf,
//  but leaf growth and mass should be baesd on "one single leaf"
//  Thus, we need to make a conversion between "Representative" and "one single leaf" using livingfrac
void RyeLeaf::RyeLeafMassUpdate(double biomassIncomeRate, double nitrogenIncomeRate)
{
	//Z reset for each time step
	//  store biomass (g) and N (mg) allocated from plant
	LeafBiomassIncrease = 0.0;
	SheathBiomassIncrease = 0.0;
	LeafNitrogenIncrease = 0.0;
	SheathNitrogenIncrease = 0.0;

	acceleration = develop->is_accel();
	if ((!dead) && (livingFrac > 0.0))
	{
		if (ptnLfMassIncrease_Rep > 0.0)
		{
			double biomassIncomeRate_adj = biomassIncomeRate / livingFrac;
			if (acceleration)
			{
				//Z under acceleration, leaf and sheath has the same growing speed
				LeafBiomassIncrease = ptnLfMassIncrease_Rep * biomassIncomeRate_adj * 0.5;		// mass in g
				SheathBiomassIncrease = ptnLfMassIncrease_Rep * biomassIncomeRate_adj * 0.5;	// mass in g
			}
			else
			{
				//Z before acceleration, leaf growing speed = 2 * sheath growing speed
				LeafBiomassIncrease = ptnLfMassIncrease_Rep * biomassIncomeRate_adj * 0.667;		// mass in g
				SheathBiomassIncrease = ptnLfMassIncrease_Rep * biomassIncomeRate_adj * 0.333;		// mass in g
			}
			//Z update leaf and sheath bio mass, N mass
			LeafMass = LeafMass + LeafBiomassIncrease;
			SheathMass = SheathMass + SheathBiomassIncrease;
			GreenLfMass = GreenLfMass + LeafBiomassIncrease;
			GreenSheathMass = GreenSheathMass + SheathBiomassIncrease;
		}
		if (ptnLfNitrogenMassIncrease_Rep > 0.0)
		{
			double nitrogenIncomeRate_adj = nitrogenIncomeRate / livingFrac;
			if (acceleration)
			{
				//Z under acceleration, leaf and sheath has the same growing speed
				LeafNitrogenIncrease = ptnLfNitrogenMassIncrease_Rep * nitrogenIncomeRate_adj * 0.5;	// mass in mg
				SheathNitrogenIncrease = ptnLfNitrogenMassIncrease_Rep * nitrogenIncomeRate_adj * 0.5;	// mass in mg
			}
			else
			{
				//Z before acceleration, leaf growing speed = 2 * sheath growing speed
				LeafNitrogenIncrease = ptnLfNitrogenMassIncrease_Rep * nitrogenIncomeRate_adj * 0.667;		// mass in mg
				SheathNitrogenIncrease = ptnLfNitrogenMassIncrease_Rep * nitrogenIncomeRate_adj * 0.333;	// mass in mg
			}
			LeafNitrogenMass = LeafNitrogenMass + LeafNitrogenIncrease;
			SheathNitrogenMass = SheathNitrogenMass + SheathNitrogenIncrease;
		}
		//Z update leaf and sheath N content
		//Z OK to over the total leafmass since for a single plant leaf, at this time, leaf is still growing
		if ((ptnLfMassIncrease_Rep > 0.0) || (ptnLfNitrogenMassIncrease_Rep > 0.0)) 
		{
			LeafNitrogenContent = 0.10 * LeafNitrogenMass / LeafMass; // 0.1=100/1000, 100 converts fraction to % and 1000 convert mg to g
			SheathNitrogenContent = 0.10 * SheathNitrogenMass / SheathMass;
			LeafNitrogenContent = __max(__min(LeafNitrogenContent, MAX_N_PCT), 0.0);
			SheathNitrogenContent = __max(__min(SheathNitrogenContent, MAX_N_PCT), 0.0);
		}
	}
	return;
}

//Z leaf death function, this is a manually call leaf death function
//  a combination of leaf senecence (to the end) and leaf drop
//  call when tiller dies externally, then will affect the leaf update function
void RyeLeaf::RyeLeafDead()
{
	LeafNitrogenRelease = 0.0;
	SheathNitrogenRelease = 0.0;
	ptnLfAreaDecrease = 0.0;
	LfNitrogenDecline = 0.0;
	SheathNitrogenDecline = 0.0;

	//Z prevent kill the leaf twice
	if (!dead)
	{
		//Z make the leaf senecence to the end, 100%
		if (mature) {
			ptnLfAreaDecrease = __max(0.0, (MaturedLfArea - SenescentLfArea));
		}
		else {
			ptnLfAreaDecrease = __max(0.0, LeafArea);
		}

		SenescentLfArea = SenescentLfArea + ptnLfAreaDecrease;

//Z DO NOT change "livingFrac" here
//  DO NOT set "GreenLf/Width/Length, LeafNitrogenMass, SheathNitrogenMass" to 0
//  That is because we need those values to calculate N return and convert greenlf to droplf in the "LivingFractionAdjustment" function
//  The correctness can be viewed by computing N mass fraction from the output.
//  The N mass% ranges 0.25% to 4.0%

		//Z determine N release/decline(to dead plant tissue) in mg
		if (LeafNitrogenContent <= LfNitrogenReleaseThreshold)
		{
			LeafNitrogenRelease = 0.0;
			SheathNitrogenRelease = 0.0;
			LfNitrogenDecline = LeafNitrogenMass;
			SheathNitrogenDecline = SheathNitrogenMass;

			DeadLeafNitrogenMass += LfNitrogenDecline;
			DeadSheathNitrogenMass += SheathNitrogenDecline;
		}
		else
		{
			//Z use "LfNitrogenMaxReleasePtge" to define a temporary fraction for 
				//  N decline -> stay in the dead/yellow organs or residue
				//  N release -> reused by other organs
			double NitrogenReleaseFrac = (LeafNitrogenContent - LfNitrogenReleaseThreshold) / LeafNitrogenContent;
			NitrogenReleaseFrac = __max(__min(NitrogenReleaseFrac, LfNitrogenMaxReleasePtge / 100.0), 0.0);
			double NitrogenDeclineFrac = 1.0 - NitrogenReleaseFrac;

			LeafNitrogenRelease = LeafNitrogenMass * NitrogenReleaseFrac;
			SheathNitrogenRelease = SheathNitrogenMass * NitrogenReleaseFrac;
			LfNitrogenDecline = LeafNitrogenMass * NitrogenDeclineFrac;
			SheathNitrogenDecline = SheathNitrogenMass * NitrogenDeclineFrac;

			DeadLeafNitrogenMass += LfNitrogenDecline;
			DeadSheathNitrogenMass += SheathNitrogenDecline;
		}


		//Z make the leaf drop
		DropLfMass = LeafMass;
		DropSheathMass = SheathMass;
		DropLfArea = SenescentLfArea;

		//Z set every growing stage to true, because we do not need this leaf anymore
		//  the last one mark the leaf death is forced by tiller death
		mature = true;
		aging = true;
		dead = true;
		dropped = true;
		force_to_death_current_step = true; //Z this only used when senecent function exists, so not for internode yet
	}

//Z DO NOT change "livingFrac" here
//  DO NOT set "GreenLf/Width/Length, LeafNitrogenMass, SheathNitrogenMass" to 0
//  That is because we need those values to calculate N return and convert greenlf to droplf in the "LivingFractionAdjustment" function
//  The correctness can be viewed by computing N mass fraction from the output.
//  The N mass% ranges 0.25% to 4.0%

	//Z Setting necessary values back to 0 can be done in "LivingFractionAdjustment" function
	//  under the condition that "drop==ture".
}

//Z tmpr effects on growth from K.Paff
double RyeLeaf::TmprEffect(double Tmpr)
{
	//These parameters were estimated for various spring wheat from the International Heat Stress Genotype Experiment 
	double AlphaGrowth = log(2.0) / (log((LeafExpan_TmprMax - LeafExpan_TmprMin) / (LeafExpan_TmprOpt - LeafExpan_TmprMin)));
	double Growth_Temperature_Effect = 0.0;
	//If the average daily temperature is outside of the acceptable range.
	if (Tmpr < LeafExpan_TmprMin || Tmpr > LeafExpan_TmprMax) {
		Growth_Temperature_Effect = 0.0;
	}
	else {
		double temp1 = pow((Tmpr - LeafExpan_TmprMin), AlphaGrowth);
		double temp2 = pow((LeafExpan_TmprOpt - LeafExpan_TmprMin), AlphaGrowth);
		Growth_Temperature_Effect = (2.0 * temp1 * temp2 - pow(temp1, 2.0)) / pow(temp2, 2.0);
	}
	Growth_Temperature_Effect = __max(Growth_Temperature_Effect, 0.0);
	Growth_Temperature_Effect = __min(Growth_Temperature_Effect, 1.0);
	return Growth_Temperature_Effect;
}
