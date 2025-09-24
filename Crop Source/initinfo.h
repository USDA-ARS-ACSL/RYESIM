#pragma once
#ifndef _INITINFO_H_
#define _INITINFO_H_
#define MINUTESPERDAY (24.0*60.0)
#ifndef FLOAT_EQ
#define EPSILON 0.001   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#endif

struct TInitInfo
{
public:
	TInitInfo()
	{
		char description[255] = "\0";
		char cultivar[255] = "\0";
		GDD_rating = 1331;
		genericLeafNo = 15;
		latitude = 38.0; longitude = 0.0; altitude = 50.0;
		sowingDay = 150;
		beginDay = 1; endDay = 365;
		year = 2004;
		timeStep = 5.0;
		plantDensity = 8.0;
		CO2 = 370.0;
		Rmax_LIR = .0978;
		Rmax_LTAR = 0.53;
		DayLengthSensitive = true;
		PhyllochronsToSilk = 8;
		PhyllochronsToTassel = 1;
		stayGreen = 4.5;
		LM_min = 125.0;

		//Z make the following variables in this structure for ryesim
		fpibpa = srpa = drpa = tspa = iepa = jtpa = bootpa = headpa = antspa = antepa = matpa = 0.0;
		seedDepth = 0.0;

		gdd_base = 0.0;
		gdd_opt = 26.0;
		gdd_max = 34.0;

	};
	char description[255];
	char cultivar[255];
	int GDD_rating; // GDD or GTI rating of the cv, see Stewart 1999 for conversion between MRMR and other ratings
	int genericLeafNo; // leaf number at the end of juvenile phase independent of environmental ques of leaf initiation
	double plantDensity;
	double latitude, longitude, altitude;
	int sowingDay, beginDay, endDay;
	double CO2;
	int year;
	double timeStep;
	bool DayLengthSensitive; //1 if daylength sensitive
	double Rmax_LIR, Rmax_LTAR; //  Maximum Leaf tip initiation and appearance rates
	double stayGreen;  // staygreen trait of the hybrid (originally 4.5)
	double LM_min; //Length of largest leaf
                                  // stay green for this value times growth period after peaking before senescence begins
                                  // An analogy for this is that with no other stresses involved, it takes 15 years to grow up, stays active for 60 years, and age the last 15 year if it were for a 90 year life span creature.
	                              //Once fully grown, the clock works differently so that the hotter it is quicker it ages
	double PhyllochronsToSilk; //number of phyllochrons from tassel initiation for 75% silking.
	double PhyllochronsToTassel; // number of phyllochrons past tassel initiation when tassels are fully emerged. (not input yet)
	//todo these above 2 variables are also in development - need to remove them from there.
	//check units

	// ***************************************************
	//Z make the following variables in this structure for ryesim

	//Z phyllochrons numbers for cereal phenology, the last one will be based on gdd
	double fpibpa;	//number of phyllochrons after double ridge that flower primordium initiation begins.
	double srpa;	//number of phyllochrons from vernalization to single ridge.
	double drpa;	//number of phyllochrons between singleridge and double ridge.
	double tspa;	//number of phyllochrons between double ridge and terminal spikelet growth stages.
	double iepa;	//number of phyllochrons between double ridge and start of internode elongation growth stages.
	double jtpa;	//number of phyllochrons between beginning of internode elongation and jointing.
	double bootpa;	//number of phyllochrons between the beginning of jointing and booting (i.e., flag leaf fully grown).
	double headpa;	//number of phyllochrons between the beginning of booting and heading.
	double antspa;	//number of phyllochrons between the beginning of heading and beginning of anthesis.
	double antepa;	//number of phyllochrons between the beginning of anthesis and end of anthesis (duration of anthesis).
	double matpa;	//number of growing degree-days between the beginning of anthesis and physiological maturity.

	double seedDepth;  //seeding depth

	double gdd_base; //base tmpr for computing gdd
	double gdd_opt;  //optimal tmpr for computing gdd
	double gdd_max;  //max tmpr for computing gdd, beyond that, no gdd and plant stop morph growth

};
#endif