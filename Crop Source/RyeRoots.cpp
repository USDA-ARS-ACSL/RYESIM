#include "stdafx.h"
#include "RyeRoots.h"
#include "weather.h"
#include "initinfo.h"
#include <cmath>
#include <algorithm>
#define DAYPERMINUTES 0.00069444444
#define MAXLEAFNUM 20
#define PHYLLOCHRON 106.0
#define CO2_MW 44.0098
#define C_MW 12.011
#define CH2O_MW 30.03

using namespace std;

//Z constructor
RyeRoots::RyeRoots()
{
	//Z root initializaiton
	rootInitialization = false;

	//Z root age
	RootsAge = 0.0;
	physAge = 0.0;

	//Z root mass managements
	RootMass = 0.0;
}

RyeRoots::~RyeRoots() {}