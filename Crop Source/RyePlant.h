#pragma once
#ifndef _RYE_PLANT_H_
#define _RYE_PLANT_H_

#include "RyeDevelopment.h"
#include "RyeTiller.h"
#include "RyeRoots.h"
#include "gas_exchange.h"
#include "initinfo.h"
#include "weather.h"
#include <iostream>
#include <string>

class RyePlant
{
public:
	RyePlant(const TInitInfo&, TGasExSpeciesParam&);
	~RyePlant();
	void RyePlantUpdate(const TWeather& weather);
	void RyeRecursiveUpdateSummary();

	//Z for each rye plant, this only return the mainstem
	RyeTiller* get_mainstem() { return mainstem; };
	//Z get rye roots
	RyeRoots* get_roots() { return roots; }
	//Z get rye development
	RyeDevelopment* get_develop() { return develop; }
	//Z return two pointers that point to a sunlit and a shaded leaf
	CGasExchange* get_sunlit() { return this->sunlit; }
	CGasExchange* get_shaded() { return this->shaded; }

//********** Z Plant manupulate tiller update and summary **********************
//Z tiller (including leaf/internode) will not do anything until plant call them to do it
//  Recursive UPDATE and SUMMARY of leaf mass and internode mass, call by plant (upper class) and act on tiller (lower class)
//  Plant class has authority on tiller class
	void PltCallTillerRecursiveUpdate(RyeTiller* ryeTiller);
//Z mass distribution among the tillers
//  plant class (upper class) will call and authorize such action
//  then the tiller "named" by the plant, following the linked list, will call its own mass distribution function for growth of leaves and internodes
	void PltCallTillerRecursiveMassDistribution(RyeTiller* ryeTiller);

	//Z IO function for GET
//********** Z reorganize the code ****************
	double get_leafNum() { return LeafNumPlt; }
	double get_greenLfNum() { return GreenLfNumPlt; }
	double get_tillerNum() { return TillerNumPlt; }						//Z living tiller number
	double get_plantLivingFrac() { return PlantLivingFraction; }		//Z not all seeds are germinated or emerged, there is a fraction
	double get_leafArea() { return LeafAreaPlt; }						//Z plant scale leaf area in cm^2
	double get_greenLeafArea() { return GreenLfAreaPlt; }				//Z return green leaf area, leaf area = green + aged (not dropped)
	double get_LAI() { return LAI; }
	double get_SeedMass() { return SeedMass; }

	double get_plantMass() { return RyeMass; }							//Z rye plant mass,   g/plant
	double get_plantTotalNitrogen() { return NitrogenMassPlt; }			//Z Total N in plant, mg/plant
	double get_grossPhotosynthesis() { return photosynthesis_gross; }	//Z Gross Photosynthesis, grams biomass per plant per hour
	double get_netPhotosynthesis() { return photosynthesis_net; }		//Z Net Photosynthesis, grams biomass per plant per hour
	double get_ET() { return transpiration; }							//Z transpiration, gr per plant per hour 
	double get_ET_Old() { return transpirationOld; }					//Z previous transpiration, gr per plant per hour 
	double get_plantTmpr() { return temperature; }						//Z plant temperature
	double get_MaintenanceRespiration() { return maintRespiration; }	//Z maintenance respiration g biomass/plant/hour
	double get_stemMass() { return IntrNodeMassPlt; }					//Z stem mass is the internode mass, since there may be multiple tillers, the existing part at this time step, not include potential growth
	double get_leafMass() { return LeafMassPlt; }						//Z leaf mass per plant
	double get_shootMass() { return ShootMass; }						//Z leaf mass + sheath mass + internode mass
	double get_rootMass() { return RootMass; }							//Z root mass, need to be assigned from 2DSOIL

	double get_leafPart() { return leafPart; }							//Z leaf biomass assignment g/plant
	double get_actualRootBiomassAssignment_PCRL() { return ActuralRootBiomassAssignment_PCRL; }	//Z root biomass assignment g/plant, already adjusted with the BiomassRoot pool, different from maizsim (which done in crop)
	double get_actualShootBiomassAssignment() { return ActuralShootBiomassAssignment; }			//Z shoot biomass assignment g/plant
	
	double get_conductance() { return conductance; }					//Z leaf surface conductance, output from photosynthesis
	double get_VPD() { return VPD; }									//Z leaf surface VPD, output from photosynthesis
	double get_leafNitrogenMass() { return LeafNitrogenMassPlt; }       //Z plant scale leaf nitrogen mass mg/plant

	double get_cumulativeNitrogenDemand() { return CumulativeNitrogenDemand; }			//Z cumulated N demand to be assigned externally, "gN"/plant, 
	double get_cumulativeNitrogenSoilUptake() { return CumulativeNitrogenSoilUptake; }	//Z cumulated N uptake to be assigned externally, "gN"/plant, 
	double get_droppedLfArea() { return DropLfAreaPlt; }				//Z plant scale drop leaf area, cm^2/plant
	//double get_ptnLfIncrease() { return PtnLeafAreaIncreasePlt; }		//Z plant scale leaf area increase at this time step, cm^2/plant

	//Z photosynthesis variables to be output
	double get_sunlit_LAI() { return sunlit_LAI; }
	double get_shaded_LAI() { return shaded_LAI; }
	double get_sunlit_PFD() { return sunlit_PFD; }
	double get_shaded_PFD() { return shaded_PFD; }
	double get_sunlit_A_net() { return sunlit_A_net; }
	double get_shaded_A_net() { return shaded_A_net; }
	double get_sunlit_A_gross() { return sunlit_A_gross; }
	double get_shaded_A_gross() { return shaded_A_gross; }
	double get_sunlit_gs() { return sunlit_gs; }
	double get_shaded_gs() { return shaded_gs; }
	double get_C2Effect() { return C2_effect; }
	double get_SunlitRatio() { return SunlitRatio; }
	double get_MaintRespiration() { return maintRespiration; }

	//Z IO function for biomass pools
	//  these function should be used for model evaluation purpose, maybe not that useful for model application
	//  biomass pools are more like internal variables
	double get_biomassRootPool() { return BiomassRoot; }				//Z root biomass storage, partially release to root biomass at current time g/plant
	double get_biomassReserve() { return BiomassReserve; }				//Z long term biomass reservior before allocation. 
	double get_nitrogenPool() { return NitrogenPool; }
	double get_HourlyNitrogenSoilUptake() { return HourlyNitrogenSoilUptake; }
	double get_LeafNitrogenRelease() { return LeafNitrogenReleasePlt; }
	double get_SheathNitrogenRelease() { return SheathNitrogenReleasePlt; }
	double get_IntrNodeNitrogenRlease() { return IntrNodeNitrogenReleasePlt; }
	double get_PltShootNitrogenAssignment() { return ShootNitrogenAvailiableAllocation; }
	double get_PltRootNitrogenAssignment() { return RootNitrogenAvailiableAllocation; }
	double get_BiomassRootAllocationDuringGermination() { return BiomassRootAllocationDuringGermination; }	//Z cumulative biomass allocated to roots during germination g/plant
	

	//cZ IO function for set
	void set_HourlyNitrogenDemand(double x) { HourlyNitrogenDemand = x; }
	void set_CumulativeNitrogenDemand(double x) { CumulativeNitrogenDemand = x; }
	void set_HourlyNitrogenSoilUptake(double x) { HourlyNitrogenSoilUptake = x; }
	void set_CumulativeNitrogenSoilUptake(double x) { CumulativeNitrogenSoilUptake = x; }
	void set_NitrogenRatio(double x) { NitrogenRatio = x; }
	void add_BiomassRootAllocationDuringGermination(double x) { BiomassRootAllocationDuringGermination += x; }	//Z set cumulative biomass allocated to roots during germination g/plant

	void calcGasExchange(const TWeather& weather, const TGasExSpeciesParam& photoparam);
	void calcMaintRespiration(const TWeather&);
	void calcBiomassAllocation(const TWeather&);
	void calcNitrogenAllocation();
	void calcSetMass();
	void calcRed_FRedRatio(const TWeather&);

	double get_averagedBiomassLeftover() { return BiomassLeftover; }

private:
	TInitInfo initInfo;
	TGasExSpeciesParam gasExparam;
	//Z for rye plant, it only trace the mainstem
	RyeTiller* mainstem;
	//Z rye roots
	RyeRoots* roots;
	//Z rye development object
	RyeDevelopment* develop;
	//Z declare two pointers that point to a sunlit and a shaded leaf
	CGasExchange* sunlit;
	CGasExchange* shaded;
	
	//***** Plant Mass Group (g biomass) *********
	double RyeMass;
	double SeedMass, SeedNitrogenMass;
	double ShootMass;
	double RootMass;
	double SeedRootMass;
	double PlantLivingFraction;    //Z used to adjust plant density, popslab ... due to not all seeds become living plant

	//***** Plant Leaf Area Group (cm^2 or g) *********
	//Z Dropped (dead) leaves must be Senescent leaves, so SenescentLfAreaPlt include dropped parts
	//  Senescent parts may still stay on the plant, so Dropped leaf will not be all the senescent leaf
	//  Thus, DropXXX is part of the total XXX
	double LeafNumPlt;
	double GreenLfNumPlt;
	double DropLfNumPt;

	double LeafAreaPlt;
	double GreenLfAreaPlt;
	double GreenLfLengthPlt;
	double GreenLfWidthPlt;
	double SenescentLfAreaPlt;
	double DropLfAreaPlt;			 //Z leaf area dropped must be senecent leaf, but senecent (part) of the leaf may not drop
	double LAI; //Leaf area Index
	
	double LeafMassPlt;
	double SheathMassPlt;
	double DropLeafMassPlt;
	double DropSheathMassPlt;
	double LeafNitrogenMassPlt;
	double SheathNitrogenMassPlt;
	double LeafNitrogenReleasePlt;
	double SheathNitrogenReleasePlt;

	//Z for dead plant nitrogen mass, in dead organ and cannot be recycled, we further partition that into
	//  "Yellow": N in the senecent portion but not dropped from the stem
	//  "Dropped" N in the dropped leaf, so already in the residue
	//  Note that for nitrogen: 
	//     Dead = Yellow + Dropped
	double LeafNitrogenDeadMassPlt;
	double SheathNitrogenDeadMassPlt;
	double LeafNitrogenYellowMassPlt;
	double SheathNitrogenYellowMassPlt;
	double LeafNitrogenDroppedMassPlt;
	double SheathNitrogenDroppedMassPlt;

	double PtnLeafMassIncreasePlt;			//Z ptn=potential 
	double PtnLeafNitrogenMassIncreasePlt;

	//***** Plant Internode Group (cm or g) *********
	int TillerNumSingle;
	double TillerNumPlt;
	double IntrNodeLengthPlt;
	double IntrNodeMassPlt;
	double IntrNodeNitrogenMassPlt;
	double IntrNodeNitrogenReleasePlt;

	double IntrNodeNitrogenDeadMassPlt;

	double PtnIntrNodeMassIncreasePlt;		
	double PtnIntrNodeNitrogenMassIncreasePlt;

	//**** Plant Root N
	double RootNitrogenMassPlt;
	double NitrogenRootGrowth_ptn;

	//***** Plant Based N group *********
	double NitrogenMassPlt;					// TotalNitrogen per plant, mg
	double NitrogenRatio;					// optimal N ratio according to N Dilution ratio	
	double LeafNitrogenContentPlt;			// averaged leaf N content

	//***** PhotoSynthesis Group *********
	double sunlit_LAI, shaded_LAI;			// sunlit and shaded LAI values
	double sunlit_PFD, shaded_PFD;			// sunlit and shaded PFD (umol m-2 s) PFD Photon Flux Density
	double sunlit_A_gross, shaded_A_gross;	// leaf area based gross assmilation for one unit sunlit and shaded LAI
	double sunlit_A_net, shaded_A_net;		// leaf area based net assmilation for one unit sunlit and shaded LAI
	double assimilate;						// assimilation flux, gCO2/plant/timstep (hour)
	double photosynthesis_gross;			// gross photosynthesis (biomass version of the assimilation), unit converted from "umol CO2 m-2 s-1" to "g Biomass plant^-1 hour^-1"
	double photosynthesis_net;				// net photosynthesis, unit converted from "umol CO2 m-2 s-1" to "g Biomass plant^-1 hour^-1"
	double transpiration;					// current and previous values of transpiration, g/plant/hr
	double transpirationOld; 
	double temperature;						// leaf temperature in the photosynthesis model
	double VPD;								// vapor pressure deficit
	double sunlit_gs, shaded_gs;			// sunlit and shaded stomatal conductance
	double conductance;						// averaged stomatal conductance

	//***** Carbon/Nitorgen Partitioning and Stroage (C in g; N in mg) *********

	double BiomassReserve;					//Z long term biomass pool in g per plant
	double BiomassPool;						//Z short term biomass pool in g per plant
	double BiomassSupply;					//Z biomass supply to plant growth in g per plant
	double BiomassDemand;					//Z biomass demand for REPRODUCTIVE GROWTH in g per plant
	double BiomassLeftover;                 //Z BiomassPool+BiomassReserve and averaged over the green leaf area
	double BiomassRoot;						//Z biomass root pool in g per plant
	double BiomassLeafGrowth_ptn;			//Z potential (ideal) leaf mass increases g per plant
	double BiomassIntrNodeGrowth_ptn;		//Z potential (ideal) internode mass increases g per plant
	double BiomassShootGrowth_ptn;			//Z potential (ideal) leaf and internode mass increases g per plant
	double BiomassPltGrowth_ptn;			//Z potential shoot + root mass increases g per plant
	double ActuralRootBiomassAssignment_PCRL;
	double ActuralShootBiomassAssignment;
	double BiomassRootAllocationDuringGermination;  // used to proportion more biomass to roots as seeds germinate. this is the cumulative amount

	double NitrogenPool;					//Z nitrogen pool in mg per plant

	double shootPart;						//g per plant Carbohydrate partitioined to shoot
	double rootPart;						//g per plant Carbohydrate partitioned to root
	double shootPart_old;					//shoot/root in the previous time step
	double rootPart_old;
	double leafPart;						//g per plant biomass partitioned to leaf (in ryesim, this include sheath)
	double internodePart;

	double leafPartNitrogen;
	double internodePartNitrogen;

	double LeafBiomassAssignmentRate;
	double IntrNodeBiomassAssignmentRate;
	double LeafNitrogenAssignmentRate;
	double IntrNodeNitrogenAssignmentRate;
	double maintRespiration;

	double HourlyNitrogenDemand;			// Nitrogen demand in mg N plant-1
	double CumulativeNitrogenDemand;		// cumulativeNitrogen demand in mg N plant-1
	double HourlyNitrogenSoilUptake;		// Nitrogen uptake from the soil mg N plant-1
	double CumulativeNitrogenSoilUptake;	// Nitrogen uptake from the soil mg N plant-1
	double NitrogenLeafGrowth_ptn;
	double NitrogenIntrNodeGrowth_ptn;
	double NitrogenShootGrowth_ptn;			// Based on the potential increase of the plant mass (leaf and internode), there should be a potential N requirement in mg N plant^-1
	double ShootNitrogenAvailiableAllocation;
	double RootNitrogenAvailiableAllocation;

	//***** Plant Information *********
	double sowingDay;
	double age;
	double emerge_gdd;						//records thermal time needed for plant to emergy
	double SunlitRatio;						// daily ratio of sunlit leaf area to use for scaling leaf expansion due to carbon stress
	double C2_effect;

};
#endif