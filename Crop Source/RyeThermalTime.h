#pragma once
#ifndef _RYE_THERMALTIME_H_
#define _RYE_THERMALTIME_H_

class RyeThermalTime
{
public:
	RyeThermalTime(void);
	~RyeThermalTime(void);
	void initialize(double timestep, double gddbase, double gddopt, double gddmax);
	//Z: IO for thermal time
	double get_Tcur() { return airTmprCur; }			//Z get current temperature
	double get_sum_plant() { return sum_plant; }		//Z get cumulated thermal time "observed by plant"
	double get_sum_ambient() { return sum_ambient; }	//Z get cumulated thermal time purely based on weather
	double get_currentTT() { return TTd_Cur; }			//Z get current thermal time for the current hour
	double get_movingAveTmpr() { return Tave; }

	//Z: gdd computation
	void add(double x, double y, bool z);
	void update(double airTmpr, double airTmax, double airTmin, double Daylength, double soilTmpr, bool elongation);
	void TmprMovingAve(double airTmpr);
	
	void FactorPhotoperiod(double Daylength);
	void FactorVernalisation(double airTmax, double airTmin);

private:
	double airTmprCur;		//Z current temperature
	double soilTmprCur;
	bool elongation;

	double Tbase;		//Z based temperature for gdd computation
	double Topt;		//Z optimal temperature for gdd computation, plant growth
	double Tmax;		//Z max temperature for gdd computation, beyond that, plant growth stops and gdd=0

	double T_diff_bo;
	double T_diff_mo;

	double TTd_Cur;		//Z current gdd at this hour period
			
	double dTmpr;		//Z "apparent tmpr" based on tmpr, opt tmpr, ...,
	double timeStep;	//Z timestep, should be 1/24 day
	double decreasingBeta; //Z Topt<Tcur<Tmax, the slope for gdd decreasing

	//Z sum of the gdd based on our "peak function" for thermal time
	double sum_plant;   // the sum "observed" by the plants, including photoperiod and vernalisation factor
	double sum_ambient; // the gdd based on weather

	//Z add daily max and min temperature for Vernalisation
	//Z add daylength for photoperiod
	double TDmax, TDmin, DayL;
	//Z add photoperiod and vernalisation factors
	double fD, fV, VernalTotal;

	//Z do a 5-day moving average of airtmpr
	double Tave;
	double TRecord[120];
	int Tcurrentposition;
	double TRecordNumber;
};
#endif