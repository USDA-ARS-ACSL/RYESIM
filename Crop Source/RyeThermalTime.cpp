
#include "stdafx.h"
#include "RyeThermalTime.h"

#include <iostream>
using namespace std;

#define DAYPERMINUTES 0.00069444444

RyeThermalTime::RyeThermalTime(void)
{
	airTmprCur = 25.0;
	soilTmprCur = 25.0;
	elongation = false;

	Tbase = 0.0;
	Topt = 22.0;
	Tmax = 30.0;

	T_diff_bo = Topt - Tbase;
	T_diff_mo = Tmax - Topt;

	TDmax = 0.0;
	TDmin = 0.0;
	DayL = 0.0;
	sum_plant = 0.0;
	sum_ambient = 0.0;
	dTmpr = 0.0;
	timeStep = 60.0;
	VernalTotal = 0.0;
	fD = fV = 0.0;
	TTd_Cur = 0.0;
	decreasingBeta = Topt / (Tmax - Topt);

	//Z do a 5-day moving average of air tmpr
	Tave = 0.0;
	for (int ii = 1; ii < 120; ii++) { TRecord[ii] = 0.0; }
	Tcurrentposition = 0;
	TRecordNumber = 0.0;
}

RyeThermalTime::~RyeThermalTime(void)
{
}

void RyeThermalTime::initialize(double dt, double gddbase, double gddopt, double gddmax)
{
	airTmprCur = 0.0;
	soilTmprCur = 0.0;
	elongation = false;

	TDmax = 0.0;
	TDmin = 0.0;
	DayL = 0.0;
	timeStep = dt;
	sum_plant = 0.0;
	sum_ambient = 0.0;

	//Z redo the setting for the three critical tmpr
	//  this is from the input file
	Tbase = gddbase;
	Topt = gddopt;
	Tmax = gddmax;
	decreasingBeta = Topt / (Tmax - Topt);

	add(airTmprCur, soilTmprCur, elongation);
}

void RyeThermalTime::add(double x, double y, bool z)
{
	airTmprCur = x;
	soilTmprCur = y;
	elongation = z;

	//Z timestep in develop in minutes, convert it to days
	double dD = timeStep * DAYPERMINUTES;

	//Z TTd cumulated based on air tmpr, use for comparing with weather station data
	if (airTmprCur <= Tbase)
	{
		dTmpr = 0.0;
	}
	else if (airTmprCur > Tbase && airTmprCur <= Topt)
	{
		dTmpr = airTmprCur - Tbase;
	}
	else if (airTmprCur > Topt && airTmprCur <= Tmax)
	{
		dTmpr = airTmprCur - Tbase;
	}
	else
	{
		dTmpr = 0.0;
	}
	sum_ambient += dTmpr * dD;

	//Z TTd cumulated based on air/soil tmpr, observed by plant
	double Tcur = airTmprCur;
	//if (!elongation) Tcur = __max(__max(airTmprCur, soilTmprCur), 0.0);
	if (!elongation) Tcur = __max((0.75 * airTmprCur + 0.25 * soilTmprCur), 0.0);

	if (Tcur <= Tbase)
	{
		dTmpr = 2.0;
	}
	else if (Tcur > Tbase && Tcur <= Topt)
	{
		dTmpr = __max(Tcur - Tbase, 2.0);
	}
	else if (Tcur > Topt && Tcur <= Tmax)
	{
		dTmpr = decreasingBeta * (Tmax - Tcur);
	}
	else
	{
		dTmpr = 0.0;
	}

	TTd_Cur = (dTmpr * dD) * __min(fD, fV);
	sum_plant += TTd_Cur;

}

void RyeThermalTime::FactorPhotoperiod(double Daylength)
{
	DayL = Daylength;
	fD = 1.0 - 0.002 * 1.5 * (20.0 - DayL) * (20.0 - DayL);
	fD = max(fD, 0.0);
	fD = min(fD, 1.0);
}

void RyeThermalTime::FactorVernalisation(double airTmax, double airTmin)
{
	TDmax = airTmax;
	TDmin = airTmin;
	if (TDmax <= 20 && TDmin <= 10)
	{
		double DeltaT = TDmax - TDmin;
		double aa = 1.4 - 0.0778 * airTmprCur;
		double bb = 0.5 + 13.44 * airTmprCur / (DeltaT + 3.0) / (DeltaT + 3.0);
		VernalTotal += (__min(aa, bb) / 24.0);
	}

	if (TDmax > 20 && VernalTotal <= 10)
	{
		double VT = VernalTotal;
		double aa = 0.5 * (TDmax - 20);
		VernalTotal -= (__min(aa, VT) / 24.0);
	}

	fV = 1.0 - (0.0054545 * 2.0 + 0.0003) * (50.0 - VernalTotal);
	fV = __max(fV, 0.0);
	fV = __min(fV, 1.0);

//	cout << "VernalTotal: " << VernalTotal << " fV: " << fV << '\n';
}

void RyeThermalTime::TmprMovingAve(double airTmpr)
{
	//Z do a 5-day moving average of air tmpr
	if (TRecordNumber <= 0.0)
	{
		TRecord[0] = airTmpr;
		Tave = airTmpr;
		Tcurrentposition = 1;
		TRecordNumber = 1.0;
	}
	else
	{
		double tmpr_kick_out = TRecord[Tcurrentposition];
		double T_sum = Tave * TRecordNumber;
		T_sum = T_sum - tmpr_kick_out + airTmpr;
		TRecord[Tcurrentposition] = airTmpr;

		TRecordNumber = __min(TRecordNumber + 1.0, 120.0);
		Tave = T_sum / TRecordNumber;

		Tcurrentposition = Tcurrentposition + 1;
		if (Tcurrentposition == 120) 
		{ 
			Tcurrentposition = 0; 
		}
	}
}


void RyeThermalTime::update(double airTmpr, double airTmax, double airTmin, double Daylength, double soilTmpr, bool elongation)
{
	FactorPhotoperiod(Daylength);
	FactorVernalisation(airTmax, airTmin);
	add(airTmpr, soilTmpr, elongation);
	TmprMovingAve(airTmpr);
}


