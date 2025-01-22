#pragma once
#ifndef _RYE_DEVELOPMENT_H_
#define _RYE_DEVELOPMENT_H_
#include "weather.h"
#include "initinfo.h"
#include "RyeThermalTime.h"
#include <iostream>
#include <string>
using namespace std;

struct TEvent
{
public:
	TEvent() { daytime = 0.0;  done = false; }
	double daytime;
	bool done;
};

class RyeDevelopment
{
public:
	RyeDevelopment(const TInitInfo&);
	~RyeDevelopment();
	void RyeDelpUpdate(const TWeather& wthr);
	void PlantDie(double frac);
	double EmCurv(double xval, int phase);
	double ReCurv(int phase, double yval);
	double LWPeffect(double predawn_psi_bars, double threshold);

	//Z IO functions
	RyeThermalTime get_RyeTTd() { return RyeTTd; }
	double get_Tcur() { return T_cur; }		//Z get current temperature (maybe airT or soilT)
	double get_Tair() { return T_air; }		//Z get air temperature
	double get_Tave() { return T_ave; }
	double get_TimeStep() { return TimeStep; }
	double get_ColdTime() { return Cold_Time; }
	double get_ColdTimeRatioPlant() { return Cold_Time_Ratio_Plant; }
	double get_ColdTimeRatioJoint() { return Cold_Time_Ratio_Joint; }

	double get_TTd() { return Cur_TTd; }	//Z incremental of Gdd for this time step
	double get_TTd_Plant() { return TTd_plant; }
	double get_TTd_Joint() { return TTd_joint; }
	double get_TTd_Elong() { return TTd_2_elongation; }
	double get_TTd_FlagLf_min() { return TTd_FlagLf_min; }  //Z the min flag leaf complete gdd, to compute elongation internodes
	double get_TTd_sinceSingleRidge() { return TTd_since_singleRidge; }

	//Z water effects group
	double get_PredawnLWP() { return PredawnLWP; }
	double get_tillerWaterEffect() { return tiller_water_effect; }
	double get_tillerWaterEffectAve() { return tiller_water_effect_ave; }
	double get_leafWaterEffectExpand() { return leaf_water_effect_expand; }
	double get_leafWaterEffectSenesc() { return leaf_water_effect_senesc; }
	double get_intrnodeWaterEffectGrowth() { return intrnode_water_effect_growth; }
	
	double get_shadeEffect() { return ShadeEffect; }
	double get_plantLivingFraction() { return PlantLivingFraction; }
	//Z plant tiller in development must be "single plant" rather than "representative plant"
	//  why? If we say "growing stage" for "representative plant", it is vague, because "representative plant" stands for multiple plants statistically
	//  while I reach certain "growing stage" is single plant issue
	int get_plantTillerNum() { return PlantTillerNumber; }

	int get_doy() { return doy; }

	//Yes-No stage functions
	bool is_germinInit() { return germinInit.done; }	//Z if germination starts, emergence may come before germination completes.
	bool is_germination() { return germination.done; }  //Z if germinated, yes then start root initialization
	bool is_emerge() { return emergence.done; }			//Z if emerged, yes then start plant growth
	bool is_singleRidge() { return singleRidge.done; }	//Z if single-ridge stage is reached, yes then there will be some tiller forced to die (leaf number <4)
	bool is_startEnlongation() { return elongationStart.done; }	//Z if the elongation starts
	bool is_accel() { return acceleration.done; }		//Z if single ridge is reached and growth acceleration occurs
	bool is_grainFillBegan() { return grainFillBegan.done; }
	bool is_mature() { return maturity.done; }
	bool is_dead() { return death.done; }

	//Set stage functions
	void set_maturity(bool d, double daytime) { maturity.done = d; maturity.daytime = daytime; }
	void set_death(bool d, double daytime) { death.done = d; death.daytime = daytime; }
	void set_shadeEffect(double x) { ShadeEffect = x; }
	void set_plantTillerNum(int x) { PlantTillerNumber = x; }

private:

	RyeThermalTime RyeTTd;

	//****** GDD numbers and growing stage **************
	double Cur_TTd;
	//Z germination start > 0.01 (small number)
	TEvent germinInit;
	//Z germination rate > 0.5, the whole germination will be 14 days and overlap the emergence time
	TEvent germination;
	//Z emergence rate > 0.01 (small number)
	TEvent emergence;
	//Z single ridge
	TEvent singleRidge;
	//Z elongation stage start, for the whole plant
	TEvent elongationStart;
	//Z grainFillingBegin
	TEvent grainFillBegan;
	//Z maturity
	TEvent maturity;
	//Z death
	TEvent death;
	//Z acceleration is not a stage, but we treat it as one since it triggers faster growth
	TEvent acceleration;

	//****** Environment Conditions **************
	double T_cur, T_air, T_ave;
	//Z this predawn leaf water potential is stored but not processed in this class, save and offer it to leaf class
	double PredawnLWP;
	//Z also compute a moving average for the predawn LWP for tiller water effect
	double tiller_water_effect;
	double tiller_water_effect_ave;
	double tiller_water_effect_Record[120];
	double tiller_water_effect_Record_Number;
	int tiller_water_effect_Current_Position;

	double leaf_water_effect_expand;
	double leaf_water_effect_senesc;
	double intrnode_water_effect_growth;


	//Z shaded effects, compute in plant class, assigned to development and stored, used in leaf class
	double ShadeEffect;
	//Z time step
	double TimeStep;

	//****** Germination & Emergence Group **************

	//soil water condition
	int temcon;				//Z current water content status;
	int condit;				//Z soil condition class
	int	oldcon;				//Z old soil condition class
	int doy;				//Z record doy for future usage;
	int doyRecord;          //Z record the DOY s.t. ready to update daily soil water condition

	//Z germination rate;
	double SeedGerminationRate;
	//Z emergence fraction and its real (truncated) fraction to be used outside;
	double EmergFrac, FracEmerg_Real;	
	//Z germination fraction and its real (truncated) fraction to be used outside;
	double GerminFrac, FracGermin_Real;	
	//Z seed number should be input, since there exist emergence fraction, so actual plant number changes;
	int SeedNum, GerminateSeedNum, PlantNum;
	double PlantLivingFraction;
	//Z save a tiller number for one single plant, which is related to the usage of shaded factor in tiller emerging
	//Z intuitively, single plant is used because we talk about one plant response here
	int PlantTillerNumber;
	//Z GDD for germination and emergence, TTd_seed will be from sowing.
	double TTd_seed, Cur_Germin_TTd, Cur_Emerg_TTd;
	
	//Z germination and emergence truncation limits
	double upcutg[3], locutg[3], upcute[3], locute[3]; 
	//Z mean elong rate(gdd / cm) * planting depth(cm) + mean germination rate(gdd) = memerg
	//     sgerm, semerg, std of germination and emergence gdd
	//     emercv, the covariance coefficients
	double melrat[3], mgerm[3], memerg[3], sgerm[3], semerg[3], emercv;

	//****** Growth Stage Till Jointing **************
	//Z gdd from doy==1 to till Jointing
	double TTd_joint;
	//Z gdd for plant growth, from emergence
	double TTd_plant;
	//Z gdd min for the earliest time that the flag leaf will complete its growth (MINGDF in shootgro)
	//  This is used in determining exactly when the flag leaf will complete growth (GDDF)
	//  and which will be the last internode to undergo accelerated elongation(ILAST, the peduncle).
	double TTd_FlagLf_min;
	double TTd_FlagLf_min_aux;  //Z auxiliary variable, the last leaf will appear within .5 phyllochron of boot, and will thus complete growth no earlier than .5 phyllochron after boot.
	//Z gdd to single ridge based on joint
	double TTd_2_singleridge;
	//Z gdd to elongation based on joint
	double TTd_2_elongation;
	//Z gdd since single ridge start
	double TTd_since_singleRidge;
	
	//****** Plant Death Factors **************
	//     4 fractions to mark the kill percentage range
	//     Cumulate Present value to avoid 
	double LowKill_Frac, HighKill_Frac, PresentKill_Frac, TotalKill_Frac;

	//****** cumulate the cold temperature time that the plant experienced
	double Cold_Time;
	double Cold_Time_Ratio_Plant;
	double Cold_Time_Ratio_Joint;

};
#endif