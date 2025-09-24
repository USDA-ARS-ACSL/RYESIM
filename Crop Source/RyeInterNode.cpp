#include "stdafx.h"
#include "RyeInterNode.h"
#include "weather.h"
#include "initinfo.h"
#include <cmath>
#include <algorithm>
#define DAYPERMINUTES 0.00069444444
#define MAXLEAFNUM 20
#define PHYLLOCHRON 106.0
#define PLANTHEIGHT 1.20	   // 	
#define DEFAULTINTRLENGTH 0.1  // 0.1 cm for the max internode length without elongation
#define MAX_N_PCT 3.5
#define MIN_N_PCT 0.35

RyeInterNode::RyeInterNode(int n, int o, bool m, RyeDevelopment* dv, double IntrLivingFrac)
{
	rank = n;
	order = o;
	mainstem = m;
	develop = dv;

	//Z environment and thermal time
	Cur_TTd = 0.0;
	PreDawnLWP = -0.05;

	//Z Growing Stage, start enlongation
	elongationStart = false;
	living = true;
	livingFrac = IntrLivingFrac;
	livingFrac_old = IntrLivingFrac;
	force_to_death_current_step = false;
	is_singleRidge = false;
	TmprIntrNode = 25.0;
	relative_growth = 1.0;

	//Z Internode Length (cm)
	IntrLength = 0.1;				//Z internode length cm
	ptnIntrLength = 0.1;			//Z internode length after one time step of growth cm
	maxIntrLength = 0.1;			//Z theoretical max internode length based on internode position
	ptnIntrLengthIncrease = 0.0;	//Z potential internode length increase at this time step
	ptnIntrLengthDecrease = 0.0;	//Z potential internode length decrease, will be the internode length when force to death or reserve for internode dying
	DeadIntrLength = 0.0;			//Z dead internode length, internode dies when the tiller is forced to die

	//Z Growing time/gdd measure
	IntrAge = 0.0;					//Z chronological ordinary time, day
	physAge = 0.0;					//Z gdd time (in reference to endGrowth and lifeSpan, days)
	
	//Z internode biomass (g)
	// everything below is based on 0.1 cm default internode length
	// using STL=0.015 specific leaf area (g cm^-1)
	IntrMass = 0.000015;			//Z internode mass, 
	IntrMassIncrease = 0.0;			//Z internode mass increase based on biomass assignment
	DeadIntrMass = 0.0;				//Z dead internode length, cm

	//Z internode Nitrogen
	IntrNitrogenContent = MAX_N_PCT;		//Z internode N % at no stress (%)
	IntrNitrogenMass = 0.000525;	//Z mass mg, the default value will be 0.0015*(3/100)*1000 (change g to mg)
	IntrNitrogenIncrease = 0.0;		//Z internode nitrogen mass increase based on N uptake and assignment
	IntrNitrogenRelease = 0.0;		//Z release all the internode nitrogen when force to death

	IntrNitrogenReleaseThreshold = 0.5;
	IntrNitrogenMaxReleasePtge = 80.0;
	IntrNitrogenDecline = 0.0;
	DeadIntrNitrogenMass = 0.0;

	//Z Parent tiller information, which internode should elongate
	rank_2_elongation = -1;

	//Z ****** initialize the internode length/mass after established from the tiller *******
	//Z initialize the internode length length  based on the leaf rank

	//Z potential internode length
	if (rank >= 3) {
		maxIntrLength = 1.09 * pow((double)rank - 1.0, 1.73) * PLANTHEIGHT;
		maxIntrLength = __max(maxIntrLength, DEFAULTINTRLENGTH);
	}
	else {
		// assume the internode without elongation will maintian as the default level
		// the first two internode should never > 0.1 cm
		maxIntrLength = DEFAULTINTRLENGTH;
	}

	//Z biomass and nitrogen pool that can be assigned to
	//  max carbon and nitrogen request from plant
	ptnIntrMassIncrease = 0.0;
	ptnIntrNitrogenMassIncrease = 0.0;

	//----------------------------------------------------------
	//Z initialize the representative value
	//  this changes "one single internode" to "representative internode"
	IntrDeadLength_Rep = 0.0;
	IntrLivingLength_Rep = 0.0;
	IntrLength_Rep = 0.0;

	IntrDeadMass_Rep = 0.0;
	IntrLivingMass_Rep = 0.0;
	IntrMass_Rep = 0.0;

	IntrNitrogenMass_Rep = 0.0;
	IntrNitrogenRelease_Rep = 0.0;

	DeadIntrNitrogenMass_Rep = 0.0;

	//Z all the potential growth is based on mass incremental
	ptnIntrMassIncrease_Rep = 0.0;
	ptnIntrNitrogenMassIncrease_Rep = 0.0;
}
// overload to transfer seed mass for the first tiller
RyeInterNode::RyeInterNode(int n, int o, bool m, RyeDevelopment* dv, double IntrLivingFrac, double seedMass)
      : RyeInterNode(n, o, m, dv, IntrLivingFrac)
{
	this->seedMass = seedMass;
	IntrMass = 0.05 * seedMass;

}
//Z internode will likely not die
//  even dead, from the programming prespective, death does not means the objective is deleted
//  no memory allocation need to be released at this step
RyeInterNode::~RyeInterNode() 
{
}

void RyeInterNode::RyeInterNodeLengthUpdate(void)
{
	TmprIntrNode = develop->get_Tcur();					//Z ambient temperature
	Cur_TTd = develop->get_TTd();						//Z ggd incremental at this time step (hour)
	elongationStart = develop->is_startEnlongation();	//Z if enlongation occurs, only make computation at that 
	PreDawnLWP = develop->get_PredawnLWP();				//Z predawn leaf water potential MPa
	IntrAge += develop->get_TimeStep() * DAYPERMINUTES;	//Z this is actually 1/24 (one hour), but TimeStep=60 is in minute, therefore, convert to the number of day
	physAge += Cur_TTd;

	// ----------------------
	//Z the next function is for "one single internode (not fractional)" morphology computation,
	//  but finally ajust the results based on the living fraction
	IntrLengthEnlongation();

	//Z update internode nitrogen content 
	IntrNitrogenContent = __max(__min(0.10 * IntrNitrogenMass / IntrMass, MAX_N_PCT), 0.0);
	if (IntrNitrogenMass < 0.001) {
		IntrNitrogenContent = MAX_N_PCT;
	}

	if (develop->is_singleRidge() && (!is_singleRidge)) {
		if (develop->get_ColdTime() < 16.0) {   //30.0
			relative_growth = 1.0;  //EJ
		}
		else if (develop->get_ColdTime() < 25.0) {   //30.0
			relative_growth = 2.3;  //EJ
		}
		else if (develop->get_ColdTime() < 50.0) {
			relative_growth = 2.3; 
		}
		else {
			relative_growth = 1.0;
		}
		is_singleRidge = true;
	}
	if (!develop->is_singleRidge()) {
		relative_growth = 1.4;
	}
	// ----------------------
	//Z the next function adjust the living fractions
	//  if livingFrac<livingFrac_old based on stress, then put "livingFrac_old-livingfrac" to die
	//  in principle, livingfrac<=livingfrac_old
	InterNodeLivingFractorAdjustment();
}

//Z need to note that only one internode grow, 
//     this is only for the internode length, not for mass
// 
//     elongationStart = true
//     1. leaf number is less than 2 so no internode grow
//     2. leaf number >=3, and then only leafnum-2 (the first two will be <= 0.1 cm)
//     3. second last internode (penultimate internode),  
//     4. grow the last tiller
//Z all those scenarios are analyzed in tiller class as a "governing box".
void RyeInterNode::IntrLengthEnlongation()
{
	//Z reset for each time step for potential values
	ptnIntrLengthIncrease = 0.0;
	ptnIntrLength = 0.0;
	ptnIntrMassIncrease = 0.0;
	ptnIntrNitrogenMassIncrease=0.0;

	if (living) 
	{
		//Z only call this when elongation start
		//  note that when elongation start, there will not be new tillers
		if (elongationStart)
		{

			//Z only one internode grow on one tiller per time step based on leaf and internode rank
			//  if a new leaf starts, then the enlongation for the current internode stops
			//  and the enlongation for the next internode starts
			//Z Therefore, elongation has its own starting and terminating mechanism
			//  we use tiller information to govern which internode should elongate
			//  recall tiller is a "governing box" and internode should only focus its own growth with min info exchanges with the tiller box
			if (rank == rank_2_elongation)
			{
				//Z nitrogen argument (%mass)
				//  give a mini growing rate as 10%
				double CriticalNitrogen;
				CriticalNitrogen = __max(MIN_N_PCT, IntrNitrogenContent);
				double N_effect = __max(0.0, (2.0 / (1.0 + exp(-2.9 * (CriticalNitrogen - MIN_N_PCT))) - 1.0));
				//N_effect = 1.0;
				N_effect = __min(1.0, N_effect);
				N_effect = __max(0.1, N_effect);
				N_effect = __max(0.2, sqrt(N_effect));
				//cout << "Internode: " <<N_effect << endl;

				//Z water potential argument (%mass)
				double water_effect = develop->get_intrnodeWaterEffectGrowth();

				//Z tmpr adjustment
				double Tmpr_effect = TmprEffect(TmprIntrNode);

				//Z assume internode will not be affected by shaded
				double limit_factor = __min(N_effect, __min(water_effect, Tmpr_effect));
				double gril = maxIntrLength / PHYLLOCHRON * limit_factor * relative_growth;

				ptnIntrLengthIncrease = gril * Cur_TTd;
				ptnIntrLength = IntrLength + ptnIntrLengthIncrease;

				if (ptnIntrLength >= maxIntrLength)
				{
					ptnIntrLength = maxIntrLength;
					ptnIntrLengthIncrease = ptnIntrLength - IntrLength;
					IntrLength = ptnIntrLength;
				}
				else
				{
					//Z redundant code
					//ptnIntrLengthIncrease = ptnIntrLength - IntrLength;
					IntrLength = ptnIntrLength;
				}
			}

			//Z under enlongation, consider the mass and nitrogen can be assigned
			ptnIntrMassIncrease = __max(IntrLength * stl - IntrMass, 0.0);
			ptnIntrNitrogenMassIncrease = __max(IntrLength * stl * MAX_N_PCT * 10.0 - IntrNitrogenMass, 0.0);
		}
		else
		{
			//Z Attention
			//  so far we do not have "senescent" for internode, 
			//  therefore when that internode is initiated, it always has the tend to grow, until the whole tiller is dead/killed by external stress 
			//
			//  may change the else condition to "if(!aging)" if internode senescence is included
			ptnIntrMassIncrease = __max(IntrLength * stl - IntrMass, 0.0);
			ptnIntrNitrogenMassIncrease = __max(IntrLength * stl * MAX_N_PCT * 10.0 - IntrNitrogenMass, 0.0);
		}
	}
}

/*Z living factor adjustment
	after the "single internode" morphology computation (geometry)
	use the living fraction to determine the "statistical/expectation" of the internode morphology and mass
*/
//Z we use "Rep" as "representative" or "statistically representative"
//  used for field scale mass or geometry computations
void RyeInterNode::InterNodeLivingFractorAdjustment()
{
	/*Z living fraction or living fraction old may NOT be necessary started at 1
        that is because if the tiller's living fraction is from, e.g. 0.5, 
        then even the internode is fully expanded, its living fraction is at most 0.5
		*/

	/*Z the computation is different from leaf
	    for internode is living or dead, do not include a "sencesent" period
	    */

	//Z first reset the representative mass demand and output
	ptnIntrMassIncrease_Rep = 0.0;
	ptnIntrNitrogenMassIncrease_Rep = 0.0;
	IntrNitrogenRelease_Rep = 0.0;

	if (livingFrac < livingFrac_old)
	{
		//Z livingfrac must be smaller than livingfrac old
		double livingDiff = livingFrac_old - livingFrac;

		//-------------------------
		/*Z adjust living internode geometry(length), massand potential mass required for leaf growth
		    1. first define an internode loss factor for internode length and mass, as an cumulative quantity
		    2. second define the living portion of internode length and mass
		    3. Then the living and dead portions are separated
		*/
		/*Z Why we have to do this?
		    similar to leaf, internode can dead due to tiller death at any time, so at any length
			therefore we cannot simply take "InternodeLength*livingfrac" as the representative value
			becaues some of the tiller/internode may die when they are pretty young and small
		*/
		/*Z 
		    ptnIntrLengthDecrease is reserved for internode dying, now it is always 0.
		    maybe in future put internode senecent?
		*/
		IntrDeadLength_Rep += (ptnIntrLengthDecrease * livingFrac + IntrLength * livingDiff);
		IntrLivingLength_Rep = IntrLength * livingFrac;
		IntrLength_Rep = IntrLivingLength_Rep + IntrDeadLength_Rep;

		IntrDeadMass_Rep += (DeadIntrMass * livingFrac + IntrMass * livingDiff);
		IntrLivingMass_Rep = (IntrMass - DeadIntrMass) * livingFrac;
		IntrMass_Rep = IntrLivingMass_Rep + IntrDeadMass_Rep;

		//-------------------------
		/*Z adjust internode nitrogen mass and their contribution to nitrogen releasing
		    1. mass computation is simple, since we assume only dead portion totally recycle nitrogen
		*/

		IntrNitrogenMass_Rep = IntrNitrogenMass * livingFrac;

		//-------------------------
		//Z adjust internode nitrogen mass
		//  1. note that for the N decline to dead tissues due to livingfraction changes, still need to redistribute the N
		//  2. for releasing, the living portion contribute instantaneous releasing, will the dropped one contribute the total nitrogen
		//  3. so far, "IntrNitrogenDecline" and "IntrNitrogenRelease" is reserved for internode senescence
		DeadIntrNitrogenMass_Rep += IntrNitrogenDecline * livingFrac;
		IntrNitrogenRelease_Rep = IntrNitrogenRelease * livingFrac;
		
		//  need to redistribute those N due to livingfraction changes
		double IntrNitrogenFree_temp = IntrNitrogenMass * livingDiff;
		if (IntrNitrogenContent <= IntrNitrogenReleaseThreshold)
		{
			DeadIntrNitrogenMass_Rep += IntrNitrogenFree_temp;
		}
		else
		{
			double NitrogenReleaseFrac = (IntrNitrogenContent - IntrNitrogenReleaseThreshold) / IntrNitrogenContent;
			NitrogenReleaseFrac = __max(__min(NitrogenReleaseFrac, IntrNitrogenMaxReleasePtge / 100.0), 0.0);
			double NitrogenDeclineFrac = 1.0 - NitrogenReleaseFrac;

			IntrNitrogenRelease_Rep += IntrNitrogenFree_temp * NitrogenReleaseFrac;
			DeadIntrNitrogenMass_Rep += IntrNitrogenFree_temp * NitrogenDeclineFrac;
		}


		//-------------------------
		//Z adjust potential internode growth
		//  only living portion can grow, so this part is very simple
		//ptnIntrLengthIncrease_Rep = ptnIntrLengthIncrease * livingFrac;
		ptnIntrMassIncrease_Rep = ptnIntrMassIncrease * livingFrac;
		ptnIntrNitrogenMassIncrease_Rep = ptnIntrNitrogenMassIncrease * livingFrac;

		if (!living)
		{
			livingFrac_old = 0.0;
			livingFrac = 0.0;
			IntrNitrogenMass = 0.0;

			ptnIntrMassIncrease_Rep = 0.0;
			ptnIntrNitrogenMassIncrease_Rep = 0.0;
		}

		//Z finally update the living fraction numbers for this Internode
		livingFrac_old = livingFrac;
	}
	else
	{
		//Z basically do nothing here, but still need to count for the "Rep" values
		//  1. if "livingFrac=livingFrac_old = 0.0", protection path, i.e., ensure all the data is kept or zeroed
		//  2. if "livingFrac=livingFrac_old != 0.0", count for the "Rep" values
		if ((livingFrac == 0.0) && (livingFrac_old == 0.0))
		{
			/*Z the tiller is dead
			    1. internode does not have its own senescent function, so it will die unless the tiller is dead
				2. the internode is dead, then no need for the cumulative computation
			*/
			//IntrDeadLength_Rep += (ptnIntrLengthDecrease * livingFrac + IntrLength * livingDiff);
			IntrLivingLength_Rep = 0.0;
			IntrLength_Rep = IntrDeadLength_Rep;

			//IntrDeadMass_Rep += (DeadIntrMass * livingFrac + IntrMass * livingDiff);
			IntrLivingMass_Rep = 0.0;
			IntrMass_Rep = IntrDeadMass_Rep;

			//-------------------------
			/*Z adjust internode nitrogen mass and their contribution to nitrogen releasing
				internode is dead, then all nitrogen is recycled and no future nitrogen release can be assumed
			*/
			IntrNitrogenMass_Rep = 0.0;
			IntrNitrogenRelease_Rep = 0.0;

			//-------------------------
			//Z adjust potential internode growth
			//  internode is dead, then no potential/actual growth
			//ptnIntrLengthIncrease_Rep = 0.0;
			ptnIntrMassIncrease_Rep = 0.0;
			ptnIntrNitrogenMassIncrease_Rep = 0.0;

			//make sure the N part should be reused
			IntrNitrogenMass = 0.0;
		}
		else
		{
			//-------------------------
			/*Z adjust living internode geometry(length), massand potential mass required for leaf growth
				1. first define an internode loss factor for internode length and mass, as an cumulative quantity
				2. second define the living portion of internode length and mass
				3. Then the living and dead portions are separated
			*/
			/*Z Why we have to do this?
				similar to leaf, internode can dead due to tiller death at any time, so at any length
				therefore we cannot simply take "InternodeLength*livingfrac" as the representative value
				becaues some of the tiller/internode may die when they are pretty young and small
			*/
			IntrDeadLength_Rep += (ptnIntrLengthDecrease * livingFrac);
			IntrLivingLength_Rep = IntrLength * livingFrac;
			IntrLength_Rep = IntrLivingLength_Rep + IntrDeadLength_Rep;

			IntrDeadMass_Rep += (DeadIntrMass * livingFrac);
			IntrLivingMass_Rep = (IntrMass - DeadIntrMass) * livingFrac;
			IntrMass_Rep = IntrLivingMass_Rep + IntrDeadMass_Rep;

			//-------------------------
			/*Z adjust internode nitrogen mass and their contribution to nitrogen releasing
				1. mass computation is simple, since we assume only dead portion totally recycle nitrogen
				2. for releasing, the living portion contribute instantaneous releasing, will the dropped one contribute the total nitrogen
			*/
			IntrNitrogenMass_Rep = IntrNitrogenMass * livingFrac;
			IntrNitrogenRelease_Rep = IntrNitrogenRelease * livingFrac;
			DeadIntrNitrogenMass_Rep += IntrNitrogenDecline * livingFrac;

			//-------------------------
			//Z adjust potential internode growth
			//  only living portion can grow, so this part is very simple
			//ptnIntrLengthIncrease_Rep = ptnIntrLengthIncrease * livingFrac;
			ptnIntrMassIncrease_Rep = ptnIntrMassIncrease * livingFrac;
			ptnIntrNitrogenMassIncrease_Rep = ptnIntrNitrogenMassIncrease * livingFrac;

			//Z finally update the living fraction numbers for this internode
			livingFrac_old = livingFrac;
			if (!living)
			{
				livingFrac_old = 0.0;
				livingFrac = 0.0;
				IntrNitrogenMass = 0.0;

				ptnIntrMassIncrease_Rep = 0.0;
				ptnIntrNitrogenMassIncrease_Rep = 0.0;
			}
		}
	}
}


//Z only one internode grow at one time, based on the tiller leaf number 
//  this is only for the internode mass, should be called recursively when biomass assimilation is computed in Plant class
void RyeInterNode::RyeInterNodeMassUpdate(double biomassIncomeRate, double nitrogenIncomeRate)
{
	//Z this initialization is a protection
	//  all the potential values are recommended to be reset to 0 before use and assignment
	IntrMassIncrease = 0.0;
	IntrNitrogenIncrease = 0.0;
	if ((living) && (livingFrac > 0.0))
	{
		//Z potential_rep growth of biomass and N is the only judgement for internode growth
		//  two mass filling period in this model
		//  1. when internode is initiated, fill the biomass and N to satisfy the 0.1 cm default internode length
		//  2. when elongation starts, the target elonged internode should receive biomass and N

		double biomassIncomeRate_adj = biomassIncomeRate / livingFrac;
		double nitrogenIncomeRate_adj = nitrogenIncomeRate / livingFrac;

		//Z under elongation, assign biomass and N to internode (stem)
		if (ptnIntrMassIncrease_Rep > 0.0)
		{
			IntrMassIncrease = ptnIntrMassIncrease_Rep * biomassIncomeRate_adj;		// mass in g
			IntrMass = IntrMass + IntrMassIncrease;									//Z update internode stem mass, N mass
		}
		if (ptnIntrNitrogenMassIncrease_Rep > 0.0)
		{
			IntrNitrogenIncrease = ptnIntrNitrogenMassIncrease_Rep * nitrogenIncomeRate_adj;	// mass in mg
			IntrNitrogenMass = IntrNitrogenMass + IntrNitrogenIncrease;							//Z update internode stem mass, N mass
		}
		if ((ptnIntrMassIncrease_Rep > 0.0) || (ptnIntrNitrogenMassIncrease_Rep > 0.0))
		{
			//Z update N content
			IntrNitrogenContent = 0.10 * IntrNitrogenMass / IntrMass; // 0.1=100/1000, 100 converts fraction to % and 1000 convert mg to g
		}

	}
}

/*Z internode death function, this is a manually call internode death function
	put internode to death and collect the nitrogen as released N
	Note that internode death is computed for "single internode" rather than representative value
*/
void RyeInterNode::RyeInterNodeDeath()
{
	ptnIntrLengthDecrease = 0.0;
	IntrNitrogenRelease = 0.0;
	IntrNitrogenDecline = 0.0;

	//Z prevent to kill the internode twice
	if (living)
	{
		//Z "force_to_death_current_step" seems not useful for internode
		//  however, compare to leaf, it is used to block the ordinary senescence process
		//  because "naturally senescence -> death" and force_to_death due to tiller death are Mutually Exclusive
		//Z if there is a "naturally senescence -> death" function for internode,
		//  "force_to_death_current_step" will be useful
		living = false;
		force_to_death_current_step = true;

		//Z first handle the geometry of internode death
		ptnIntrLengthDecrease = IntrLength;
		DeadIntrLength = IntrLength;

		//Z second handle the internode mass death
		DeadIntrMass = IntrMass;

		//Z little bit tricky for N
		if (IntrNitrogenContent <= IntrNitrogenReleaseThreshold)
		{
			IntrNitrogenRelease = 0.0;
			IntrNitrogenDecline = IntrNitrogenMass;

			DeadIntrNitrogenMass += IntrNitrogenDecline;
		}
		else
		{
			double NitrogenReleaseFrac = (IntrNitrogenContent - IntrNitrogenReleaseThreshold) / IntrNitrogenContent;
			NitrogenReleaseFrac = __max(__min(NitrogenReleaseFrac, IntrNitrogenMaxReleasePtge / 100.0), 0.0);
			double NitrogenDeclineFrac = 1.0 - NitrogenReleaseFrac;

			IntrNitrogenRelease = IntrNitrogenMass * NitrogenReleaseFrac;
			IntrNitrogenDecline = IntrNitrogenMass * NitrogenDeclineFrac;

			DeadIntrNitrogenMass += IntrNitrogenDecline;
		}

	}

	//Z DO NOT change "livingFrac" here
	//  DO NOT set "IntrNitrogenMass" to 0
	//  That is because we need those values to calculate N return and convert greenlf to droplf in the "LivingFractionAdjustment" function
	//  The correctness can be viewed by computing N mass fraction from the output.
	//  The N mass% ranges 0.25% to 4.0%

	//Z Setting necessary values back to 0 can be done in "LivingFractionAdjustment" function
	//  under the condition that "drop==ture".
}

//Z tmpr effects on growth from K.Paff
//  design for leaf but temporoarily use here
double RyeInterNode::TmprEffect(double Tmpr)
{
	// These parameters were estimated for various spring wheat from the International Heat Stress Genotype Experiment 
	// for the "AlphaGrowth", Z just compute it based on the leaf temperature
	// not sure if there is a function for rye/wheat internode elongation and temperture effects

	double AlphaGrowth = 1.0;
	double Growth_Temperature_Effect = 0.0;
	//If the average daily temperature is outside of the acceptable range.
	if (Tmpr < 0.0 || Tmpr > 20.0) {
		Growth_Temperature_Effect = 0.0;
	}
	else {
		double temp1 = pow(Tmpr, AlphaGrowth);
		double temp2 = pow(10.0, AlphaGrowth);
		Growth_Temperature_Effect = (2.0 * temp1 * temp2 - pow(temp1, 2.0)) / pow(temp2, 2.0);
	}
	Growth_Temperature_Effect = __max(Growth_Temperature_Effect, 0.0);
	Growth_Temperature_Effect = __min(Growth_Temperature_Effect, 1.0);
	return Growth_Temperature_Effect;
}
