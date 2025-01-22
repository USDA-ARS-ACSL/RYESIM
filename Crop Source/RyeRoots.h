#pragma once
#ifndef _RYE_ROOTS_H_
#define _RYE_ROOTS_H_

//Z I personally think this is an unnecessary class
//  it does not provide useful information
// 

class RyeRoots
{
//Z We build this class for completeness
//  This class is totally unnecessary
//  Plant.class will assign mass to the Root module in 2DSOIL directly
public:
	RyeRoots(void);
	~RyeRoots(void);

	bool get_Initialization() { return rootInitialization; }
	void set_Initialization() { rootInitialization = true; }
	void inputBiomass(double biomass) { RootMass = biomass; }

private:

	//Z initialize of not
	bool rootInitialization;

	//Z root has its own age like other organs
	double RootsAge;		// chronological ordinary time, day
	double physAge;		// gdd time (in reference to endGrowth and lifeSpan, days)

	//Z root mass managements
	double RootMass;

};
#endif