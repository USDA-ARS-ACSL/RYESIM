#include "stdafx.h"
#include <cmath>
#include "RyeDevelopment.h"
#include <iostream>
#include <string>
#include <algorithm>
//#using <mscorlib.dll>
#define PHYLLOCHRON 106.0
#define MINUTESPERDAY 1440.0
#define DAYPERMINUTES 0.00069444444
#define DAYPERHOUR 0.041666666667

using namespace std;

RyeDevelopment::RyeDevelopment(const TInitInfo& info)
{
	//Z GDD numbers and thermal time
	Cur_TTd = 0.0;
	//Z time step
	TimeStep = info.timeStep;
	RyeTTd.initialize((double)TimeStep, info.gdd_base, info.gdd_opt, info.gdd_max);

	// initialization
	//****** Germination & Emergence Group **************

	temcon = 0;
	condit = -1;			//Z -1 marks before the simulation
	oldcon = -1;			//Z -1 marks before the simulation
	doy = -1;
	doyRecord = -1;			//Z just need to make sure this marks the starting of the code, should be some nonsense number here

	SeedGerminationRate = 0.0;
	GerminFrac = 0.0, FracGermin_Real = 0.0;
	EmergFrac = 0.0, FracEmerg_Real = 0.0;
	SeedNum = info.plantDensity, GerminateSeedNum = 0, PlantNum = 0;
	PlantTillerNumber = 1.0;
	TTd_seed = 0.0, Cur_Germin_TTd = 0.0, Cur_Emerg_TTd = 0.0;

	//Z plant ratio
	//  cereal rye is different from maize, that the plant number will change over time
	//  possible reason: 1. germination/emergence have a certian fraction
	//                   2. plant may die in winter
	//Z "PlantLivingFraction" = plant number / initialized number
	//  can be used to adjust any variables like plant density, popslab ......
	PlantLivingFraction = 1.0;

	//emercv : Coefficient of variation for emergence curves.
	emercv = 0.20;
	//mgerm(1:3) : Mean degree-days needed for germination for each of the three soil water condition values.
	mgerm[0] = 80.0, mgerm[1] = 90.0, mgerm[2] = 110.0;
	//melrat(1:3) : Mean elongation rate for each of the three soil water condition values (growing degree-days/cm).
	melrat[0] = 20.0, melrat[1] = 25.0, melrat[2] = 30.0;
	//Z seedDepth (cm) in the "info" structure, emergence gdd depends on seed depth
	for (int ii = 0; ii < 3; ii++)
	{
		memerg[ii] = melrat[ii] * info.seedDepth + mgerm[ii];
		sgerm[ii] = emercv * mgerm[ii];
		semerg[ii] = emercv * memerg[ii];
		upcutg[ii] = mgerm[ii] + 3.0 * sgerm[ii];
		locutg[ii] = mgerm[ii] - 3.0 * sgerm[ii];
		upcute[ii] = memerg[ii] + 3.0 * semerg[ii];
		locute[ii] = memerg[ii] - 3.0 * semerg[ii];
	}

	//****** Growth Stage Till Jointing **************
	//Z gdd from doy==1 to till Jointing
	TTd_joint = 0.0;
	TTd_plant = 0.0;
	TTd_2_singleridge = PHYLLOCHRON * info.srpa;
	TTd_2_elongation = PHYLLOCHRON * (info.srpa + info.drpa + info.iepa);
	TTd_FlagLf_min = 0.0;
	TTd_FlagLf_min_aux = PHYLLOCHRON * (info.srpa + info.drpa + info.iepa + info.jtpa + info.bootpa + 0.5);
	TTd_since_singleRidge = 0.0;

	//****** Plant Death Factors **************
	//rye death should enable a cumulation computation
	LowKill_Frac = 0.0, HighKill_Frac = 0.0, PresentKill_Frac = 0.0, TotalKill_Frac = 0.0;

	T_cur = 0.0;
	T_air = 0.0;
	T_ave = 0.0;

	Cold_Time = 0.0;
	Cold_Time_Ratio_Plant = 1.0;
	Cold_Time_Ratio_Joint = 1.0;

	//Z this predawn leaf water potential (bar) is stored but not processed in this class, save and offer it to leaf class
	PredawnLWP = -0.05;

	//Z Water effects seems to be the same among tiller, leaf and internode organs.
	//  all of the water effects are driven by the predawn 
	//  Therefore, I put it in development module to reduce the computing load
	//Z Tiller water effects
	tiller_water_effect = 1.0;
	tiller_water_effect_ave = 1.0;
	for (int ii = 1; ii < 120; ii++) { tiller_water_effect_Record[ii] = 0.0; }
	tiller_water_effect_Current_Position = 0;
	tiller_water_effect_Record_Number = 0.0;

	//Z Leaf/Internode Growing water effects
	leaf_water_effect_expand = 1.0;
	intrnode_water_effect_growth = 1.0;

	//Z Leaf Senesc water effects
	leaf_water_effect_senesc = 1.0;

	//Z shaded effects, compute in plant class, assigned to development and stored, used in leaf class
	ShadeEffect = 1.0;
	
}

RyeDevelopment::~RyeDevelopment()
{
}

void RyeDevelopment::RyeDelpUpdate(const TWeather& wthr)
{
	//Z get necessary climate data for this module.
	T_cur = __max(0.0, wthr.airT);
	T_air = wthr.airT;
	doy = wthr.doy;
	PredawnLWP = wthr.PredawnLWP;

	//Z first update the thermal time computation for the whole plant
	// RyeTTd.update(T_cur, wthr.airTDmax, wthr.airTDmin, wthr.dayLength, wthr.soilT, elongationStart.done);
	RyeTTd.update(T_cur, wthr.airTDmax, wthr.airTDmin, wthr.dayLength, wthr.soilT, singleRidge.done);

	Cur_TTd = RyeTTd.get_currentTT();
	T_ave = RyeTTd.get_movingAveTmpr();

	//Z Germination and Emergence
	//  The criterion is if the final emergence fraction = 1.0, i.e., the whole emergence is complete
	if(FracEmerg_Real<1.0)
	{	
		//Z water saturation, need 3 cut-off values
		if (wthr.wfps >= 0.35)			//Z wfps > 0.35 (averaged degree of saturation)
		{
			temcon = 0;
		}
		else if (wthr.wfps >= 0.25)		//Z wfps = 0.25 - 0.35
		{
			temcon = 1;
		}
		else if (wthr.wfps >= 0.15)		//Z wfps = 0.15 - 0.25
		{
			temcon = 2;
		}
		else							//Z wfps < 0.15
		{
			temcon = 3;
		}

		//Z first call to set up a germination rate
		if (GerminateSeedNum == 0 && temcon != 3)
		{
			//Z make germination as a continous function w.r.t. wfps
			if (wthr.wfps >= 0.35)
			{
				SeedGerminationRate = 1.0;
			}
			else if (wthr.wfps >= 0.25)
			{
				SeedGerminationRate = 0.9 + (wthr.wfps - 0.25);
			}
			else if (wthr.wfps >= 0.15)
			{
				SeedGerminationRate = 0.8 + (wthr.wfps - 0.15);
			}
			else if (wthr.wfps >= 0.05)
			{
				SeedGerminationRate = 0.7 + (wthr.wfps - 0.05);
			}
			else
			{
				SeedGerminationRate = 0.7;
			}
			//SeedGerminationRate = 1.0 - (double)temcon * 0.1;
			GerminateSeedNum = (int)(SeedGerminationRate * (double)SeedNum + 0.5);
			PlantLivingFraction = double(GerminateSeedNum / SeedNum);
		}

		//Z start the germination and emergence procedure
		/*
		Use today and previous soil conditions to determine whether development occurs today and the algorithm to calculate it.
        Four general situations must be handled, in shootgro, temcon = 1,2,3,4 and condit = 0, while in c++, we use 0,1,2,3 and -1 (always minus 1)
			Condit = yesterday's soil condition index
			temcon = today's soil condition index

		(1) If soil conditions are too dry for any development (temcon = 3), 
		    and this is either the day of planting (condit = -1) 
            or not the first day of the dry period (condit = 3), 
            do nothing except reset condit to today's value of temcon in preparation for tomorrow. 
            
        (2) If this is the first dry day after a period of conditions wet enough to allow growth (temcon = 3, condit = 0,1,2)
            kill any plants that have germinated but not yet emerged, 
            then set oldcon to yesterday's soil condition value (condit) for
            comparison when conditions improve.

        (3) If soil conditions are adequate for seed activity (temcon = 0,1,2) 
            and there has been a change in the soil condition value since the previous day on which development occurred,
            either with (condit = 3, oldcon != temcon) or without an intervening dry period (condit = 1,2, or 3)
            recalculate progress using a new growth curve, then perform the calculations associated with growth.

        (4) If soil conditions are adequate for seed activity (temcon = 0,1,2)
            and this is the first day of activity (oldcon = -1) either with (condit = 3) or without (condit = -1) an initial dry spell, 
            or if soil conditions are the same as on the previous day on which any activity occurred, 
			either with (oldcon = temcon) or without (condit = temcon) an intervening dry spell,
            perform the calculations associated with growth (there is NO need to update the growth curve, 
			OTHERWISE will be no normal germination or emergence).
		*/

		//Z  If soil water conditions are not too dry for growth to occur today, case (3) and (4)
		//   no need to perform curve recalculation for (4)
		if (temcon != 3)
		{
			/*
			If soil conditions have changed since the previous day on which seed activity occurred
            (either with case (1) or without case (2) an intervening severely dry spell)
            recalculate the germination and combined emergence curves.
            
			Gfrac and efrac are the current cumulative fractions of plants that have germinated or emerged, respectively.
            The iterative function RECURV calculates the GDD values corresponding to gfrac and efrac on the curve associated with the current soil condition (temcon).
  
            The germination and emergence curves that result from changes in soil conditions will lie between the curves for condition 0 and 2.
			Can be shown by plotting gfrac and efrac  against gdds, with gdds on the same scale as bgddg and bgdde.  
			(Gdds is the cumulative GDD since planting, not including days at condition 3.)
			The first parameter passed to RECURV indicates whether the germination curves (1) or the emergence curves (2) are to be used.
			*/

			//Z General case 3, intervening dry spell:
			if (condit != 3 && oldcon != temcon && oldcon != -1)
			{
				if (GerminFrac < 1.0)
				{
					//first argument = 1 germination
					Cur_Germin_TTd = ReCurv(1, GerminFrac);
				}
				if (EmergFrac < 1.0)
				{
					//first argument !=1 emergence
					Cur_Emerg_TTd = ReCurv(2, EmergFrac);
				}
			}
			//Z General case 3, no intervening dry spell :
			if (condit != temcon && condit != -1 && condit != 3)
			{
				if (GerminFrac < 1.0)
				{
					//first argument = 1 germination
					Cur_Germin_TTd = ReCurv(1, GerminFrac);
				}
				if (EmergFrac < 1.0)
				{
					//first argument !=1 emergence
					Cur_Emerg_TTd = ReCurv(2, EmergFrac);
				}
			}
			
			TTd_seed += Cur_TTd;
			Cur_Germin_TTd += Cur_TTd;
			Cur_Emerg_TTd += Cur_TTd;

			if (Cur_Germin_TTd > upcutg[temcon])
			{
				GerminFrac = 1.0;
				FracGermin_Real = 1.0;
			}
			else
			{
				GerminFrac = EmCurv(Cur_Germin_TTd, 1);
				FracGermin_Real = GerminFrac;
				if (Cur_Germin_TTd < locutg[temcon])
				{
					FracGermin_Real = 0.0;
				}
			}

			if (Cur_Emerg_TTd > upcute[temcon])
			{
				EmergFrac = 1.0;
				FracEmerg_Real = 1.0;
			}
			else
			{
				EmergFrac = EmCurv(Cur_Emerg_TTd, 2);
				FracEmerg_Real = EmergFrac;
				if (Cur_Emerg_TTd < locute[temcon])
				{
					FracEmerg_Real = 0.0;
				}
			}

			//Z set up the germination and emergence marker
			//Z first germination initiation/start, use a small number here, but same as the one for emergence
			//  germination initiation must be earlier than emergence
			if (!germinInit.done && FracGermin_Real >= 0.01)
			{
				germinInit.daytime = wthr.daytime;
				germinInit.done = true;

			}
			if (!germination.done && FracGermin_Real >= 0.5)
			{
				germination.daytime = wthr.daytime;
				germination.done = true;

			}
			if (!emergence.done && FracEmerg_Real >= 0.01)
			{
				emergence.daytime = wthr.daytime;
				emergence.done = true;
			}

			if (TTd_seed >= 350)
			{
				// exceed the target emergence day, kill all the plants
				PlantDie(1.0);
			}
		}
		// If it is the first day of the dry period, but not the day of planting execute the general case (2) :
		else if (condit != -1 && condit != 3)
		{
			//Z the ambient water condition is very pool
			//  all germinated plants want to die
			PlantDie(FracGermin_Real);
			if (doy > doyRecord)
			{
				//Z doyRecord mark day + 1, where the water condition should be updated
				oldcon = condit;
				doyRecord = doy;
			}
		}

		//Z daily update the previous day water condition
		if (doy > doyRecord) condit = temcon;

		//Z based on the emergence rate
		//  a certain number of plant will grow
		if (FracEmerg_Real < 0.01)
		{
			PlantNum = 0.0;
			PlantLivingFraction = 0.001;
		}
		else if (FracEmerg_Real < 1.0)
		{
			PlantNum = (int)((FracEmerg_Real - TotalKill_Frac) * (double)GerminateSeedNum + 0.5);
			PlantNum = __max(PlantNum, 1.0);
			PlantLivingFraction = double(PlantNum) / double(SeedNum);
		}
		else 
		{
			PlantNum = (int)((1.0 - TotalKill_Frac) * (double)GerminateSeedNum + 0.5);
			PlantNum = __max(PlantNum, 1.0);
			PlantLivingFraction = double(PlantNum) / double(SeedNum);
		}
		
	}

//Z use the growth gdd (not joint gdd) to determine leaf emergence
//Z use the gdd to joint to determine start of elongation
	if (emergence.done)
	{
		//Z gdd since emergence
		TTd_plant = TTd_plant + Cur_TTd;
		if (doy < 213)
		{
			//Z in shootgro, TTd_joint start after vernalization, aka, it starts on Jan. 1 for the following year
			TTd_joint = TTd_joint + Cur_TTd;
		}

		//Z judge if the single ridge stage is reached
		if (!singleRidge.done)
		{
			if (TTd_joint >= TTd_2_singleridge)
			{
				singleRidge.daytime = wthr.daytime;
				singleRidge.done = true;
			}
		}

		//Z judge if growth accerlation should occur under the condition that the single ridge is reached
		if (!acceleration.done)
		{
			if (TTd_joint >= TTd_2_singleridge)
			{
				acceleration.daytime = wthr.daytime;
				acceleration.done = true;
			}
		}

		//Z judge if the elongation starts
		if (!elongationStart.done)
		{
			if (TTd_joint >= TTd_2_elongation)
			{
				elongationStart.daytime = wthr.daytime;
				elongationStart.done = true;
				TTd_FlagLf_min = TTd_plant - TTd_joint + TTd_FlagLf_min_aux;
			}
		}

		//Z low temperature conditions for the tiller
		double base_Tmpr = 5.0;
		if (doy < 50 || doy > 300) {
			if (T_ave < base_Tmpr) 
			{ 
				Cold_Time += DAYPERHOUR; 
				//Cold_Time_Ratio_Plant = __min(Cold_Time / TTd_plant, 1.0);
				//Cold_Time_Ratio_Joint = __min(Cold_Time / TTd_joint, 1.0);
			}
		}
		if (singleRidge.done)
		{
			TTd_since_singleRidge = TTd_since_singleRidge + Cur_TTd;
		}
		else 
		{
			Cold_Time_Ratio_Plant = __min(Cold_Time / TTd_plant, 1.0);
			Cold_Time_Ratio_Joint = __min(Cold_Time / TTd_joint, 1.0);
		}
		
		//Z start to compute water effect 
		//Z tiller, (with moving average for all the tillers)
		tiller_water_effect = 1.0;
		const double tiller_psi_threshold_bars = -0.8657;
		tiller_water_effect = LWPeffect(PredawnLWP, tiller_psi_threshold_bars);
		tiller_water_effect = __min(tiller_water_effect, 1.0);
		tiller_water_effect = __max(tiller_water_effect, 0.01);

		if (tiller_water_effect_Record_Number <= 0.0)
		{
			tiller_water_effect_Record[0] = tiller_water_effect;
			tiller_water_effect_ave = tiller_water_effect;
			tiller_water_effect_Current_Position = 1;
			tiller_water_effect_Record_Number = 1.0;
		}
		else
		{
			double water_effect_kick_out = tiller_water_effect_Record[tiller_water_effect_Current_Position];
			double water_effect_sum = tiller_water_effect_ave * tiller_water_effect_Record_Number;
			water_effect_sum = water_effect_sum - water_effect_kick_out + tiller_water_effect;
			tiller_water_effect_Record[tiller_water_effect_Current_Position] = tiller_water_effect;

			tiller_water_effect_Record_Number = __min(tiller_water_effect_Record_Number + 1.0, 120.0);
			tiller_water_effect_ave = water_effect_sum / tiller_water_effect_Record_Number;

			tiller_water_effect_Current_Position = tiller_water_effect_Current_Position + 1;
			if (tiller_water_effect_Current_Position == 120)
			{
				tiller_water_effect_Current_Position = 0;
			}
		}

		//Z Leaf expansion and Internode growth
		leaf_water_effect_expand = 1.0;
		intrnode_water_effect_growth = 1.0;
		const double psi_threshold_bars_expand = -0.8657;
		leaf_water_effect_expand = LWPeffect(PredawnLWP, psi_threshold_bars_expand);
		leaf_water_effect_expand = __min(leaf_water_effect_expand, 1.0);
		leaf_water_effect_expand = __max(leaf_water_effect_expand, 0.01);
		intrnode_water_effect_growth = leaf_water_effect_expand;

		//Z Leaf Senesc water effects
		leaf_water_effect_senesc = 1.0;
		//threshold predawn leaf water potential (in bars) below which water stress triggers senescence, needs to be substantiated with lit or exp evidence, SK
	    // This is the water potential at which considerable reduction in leaf growth takes place in corn, sunflower, and soybean in Boyear (1970)
		const double psi_threshold_bars_senesc = -4.0;
		leaf_water_effect_senesc = LWPeffect(PredawnLWP, psi_threshold_bars_senesc);
		leaf_water_effect_senesc = __min(leaf_water_effect_senesc, 1.0);
		leaf_water_effect_senesc = __max(leaf_water_effect_senesc, 0.01);
	}
}

//Z this compute plant dead fraction, from 0.0 to 1.0
void RyeDevelopment::PlantDie(double top)
{
	double lobnd = 0.0, upbnd = 1.0;
	LowKill_Frac = __max(FracEmerg_Real, HighKill_Frac);
	HighKill_Frac = top;
	if (HighKill_Frac < upbnd)
	{
		PresentKill_Frac = __min(HighKill_Frac - LowKill_Frac, HighKill_Frac - lobnd);
	}
	else
	{
		PresentKill_Frac = __min(upbnd - LowKill_Frac, upbnd - lobnd);
	}
	PresentKill_Frac = __max(PresentKill_Frac, 0.0);
	TotalKill_Frac = TotalKill_Frac + PresentKill_Frac;
	TotalKill_Frac = __min(TotalKill_Frac, 1.0);

	if ((1.0 - TotalKill_Frac) <= 0.000001)
	{
		//Z plant is totally killed, need to stop this program.
		TotalKill_Frac = 1.0;
	}
}

//Z germinaiton and emergence gdd curves
//  input gdd and obtain a fraction
double RyeDevelopment::EmCurv(double xval, int phase)
{
	double b1 = 0.319381530, b2 = -0.356563782, b3 = 1.781477937, b4 = -1.821255978, b5 = 1.330274429, p = 0.2316419, pi = 3.1415927;
	double dev = 0.0;
	if (phase == 1)
	{
		// germination
		dev = (xval - mgerm[temcon]) / sgerm[temcon];
	}
	else
	{
		// emergence
		dev = (xval - memerg[temcon]) / semerg[temcon];
	}
	double abdev = fabs(dev);
	double t1 = 1.0 / (1.0 + p * abdev);
	double yval = 1.0 - 1.0 / sqrt(2.0 * pi) * exp(-pow(abdev, 2.0) / 2.0) * t1 * (b1 + t1 * (b2 + t1 * (b3 + t1 * (b4 + t1 * b5))));

	double EmcurvRes = 0.0;
	if (dev >= 0.0)
	{
		EmcurvRes = yval;
	}
	else
	{
		EmcurvRes = 1.0 - yval;
	}
	return EmcurvRes;
}

//Z germinaiton and emergence gdd INVERSE curves
//  input a fraction and output gdd
//  with weather condition changes, given a germination/emergence fraction
//  may need to adjust the gdd needed to reach germination/emergence
double RyeDevelopment::ReCurv(int phase, double yval)
{
	int iterno = 0;
	double highx = 0.0, lowx = 0.0, medx = 0.0, tempy = 0.0, tol = 0.00001;
	highx = upcute[temcon] + 10.0;
	double RecurvRes = 0.0;
	while (true)
	{
		iterno = iterno + 1;
		medx = (highx + lowx) / 2.0;
		tempy = EmCurv(medx, phase);
		if (abs(yval - tempy) < tol)
		{
			RecurvRes = medx;
			break;
		}
		else if (yval < tempy)
		{
			highx = medx;
		}
		else
		{
			lowx = medx;
		}
	}
	return RecurvRes;
}

//Z water potential effects on growth from Boyer (1970) and Tanguilig et al (1987) YY
//   design for leaf but temporoarily use here for tiller
double RyeDevelopment::LWPeffect(double predawn_psi_bars, double threshold)
{
	double psi_f = -1.4251; // -1.0, was -1.4251   later changed to -2.3;
	double s_f = 0.4258; // was 0.4258 0.5;
	double psi_th = threshold; // threshold wp below which stress effect shows up
	double effect;
	effect = __min(1.0, (1.0 + exp(psi_f * s_f)) / (1.0 + exp(s_f * (psi_f - (predawn_psi_bars - psi_th)))));
	return effect;
}