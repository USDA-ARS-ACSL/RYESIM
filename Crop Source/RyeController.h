
//Class Controller for rye
//
// RyeController.h
//
//Based on CPM simulation_controller
#pragma once
#ifndef _RYE_CONTROLLER_H_
#define _RYE_CONTROLLER_H_
#include "timer.h"
#include "RyeDevelopment.h"
#include "RyePlant.h"
#include "weather.h"
#include "initinfo.h"
#include "gas_ex_species_param.h"
#ifndef FLOAT_EQ
#define EPSILON 0.001   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#define MINUTESPERDAY (24*60);
#endif

class RyeController
{
public:
	RyeController(const char*, const char*, const char*, TInitInfo);
	~RyeController();

	//Z read files and prepare the initilization
	void initialize();

	//Z read files from 2dsoil module
	void readWeatherFrom2DSOIL(const TWeather& wthr);

	//Z output to the crop results g01
	void outputToCropFile();

	int run(const TWeather& wthr);

//Z IO functions
	RyePlant* get_plant(void) { return ryePlant; }
	TInitInfo* get_initInfo(void) { return &initInfo; } //Z alhtough initInfo is defined as a structure, we return a pointer to avoid copy the whole structure

//Z other IO functions
	int get_sowingDay() { return SowingDay; }

//Z reserve for error catching
	int get_errStatus() { return 0; }

private:
	TInitInfo initInfo;
	TGasExSpeciesParam  GasExParam;
	RyePlant* ryePlant;
	Timer* time;
	TWeather* weather;

	//Z INPUT, name of the crop variety file
	char varietyFile[256];
	//Z OUTPUT, name of the crop output file as g01
	char cropFile[256];
	//Z OUTPUT, name of the debug file
	char DebugFile[256];
	//Z OUTPUT, name of the leaf output file
	char LeafFile[256];
	char outputFile[133];
	char logFile[256];

	//Z current record number, used to index weather files
	int iCur;

	//Z simulation time recorder
	int firstDayOfSim;
	int lastDayOfSim;
	int SowingDay;

	//Z error flag, maybe not that 
	int errorFlag;

	//Z actual supply of water to plant, mol (water) m-2 (leaf) s-1
	double ET_supply;   

	double RootWeightFrom2DSOIL;
	float MaxRootDepth;
	float AvailableWater;


};
#endif