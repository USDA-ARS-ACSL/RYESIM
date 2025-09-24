#include "stdafx.h"
#include "RyeController.h"
#include "initinfo.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#ifndef FLOAT_EQ
#define EPSILON 0.001   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#endif
#define comma ","
#ifdef _WIN32
const char* pathSymbol = "\\";
#else
const char* pathSymbol = "/";
#endif
// const a = 17.27; b = 237.7; //constant in deg C
inline double E_sat(double T) { return 0.6105 * exp(17.27 * T / (237.7 + T)); }
//#using <mscorlib.dll>
using namespace std;


RyeController::RyeController(const char* filename, const char* outfile, const char* LFile, TInitInfo iniInfo)
{
	time = NULL;
	weather = NULL;
	ryePlant = NULL;

	char* pch = (char*)calloc(256, sizeof(char));
	char* next = (char*)calloc(256, sizeof(char));
	const char* ext_dbg = ".dbg";
	char* temp = (char*)calloc(256, sizeof(char));

	//Z copy varitey input filename to the "varietyFile" buffer
	//  strcpy_s is the thread-safe version of strcpy
	strcpy_s(varietyFile, filename);

	//Z copy g01 output filename to the "cropFile" buffer
	strcpy_s(cropFile, outfile);

	// find value of cropFile(G01) strip off the path and filename in order to create a dbg file
	std::string basePath;
	std::size_t found;
	std::string cropFileAsString = cropFile;
	// find last path separator and break path from file name
	found = cropFileAsString.find_last_of(pathSymbol);
	basePath = cropFileAsString.substr(0, found);
	// now get filename to use as a root
	std::string fileName = cropFileAsString.substr(found + 1);
	std::string root = fileName.substr(0, fileName.length() - 3);
	// create summFile and debug file  names
	// first have to determine if we are linux or windows
	DebugFile = basePath.append("/" + root + ext_dbg);



	//Z copy leaf output file to the "LeafFile" buffer
	strcpy_s(LeafFile, LFile);

	//Z a local version of the initial information file
	initInfo = iniInfo;

	//Z iterator for the weather
	iCur = 0;

	//Z simulation time recorder
	firstDayOfSim = 0;
	lastDayOfSim = 365;

	//Z initialize and read files
	this->initialize();

	//Z initialize the error flag
	errorFlag = 0;

}

RyeController::~RyeController()
{
	if (time != NULL)
	{
		delete time;
		time = NULL;
	}
	if (weather != NULL)
	{
		delete[] weather;
		weather = NULL;
	}
	if (ryePlant != NULL)
	{
		delete ryePlant;
		ryePlant = NULL;
	}
	//	if ( output != NULL )
	//		delete output;
}

//Z initialize the simulation
//  mainly initialized the parameters and read fromt the input files
void RyeController::initialize()
{
	//Z define a timer to convert dates
	//  present days based on mm/dd/yy
	Timer dConvert;
	int mm, dd, yy; 
	cout << "Initializing Controller object...." << endl << endl << endl;
	//Z setiosflags is an io function to control the alignment
	cout << setiosflags(ios::left) << endl
		<< " ***********************************************************" << endl
		<< " *      RYESIM: for Cereal Rye Winter Cover Crop           *" << endl
		<< " *                 VERSION  0.3 2024                       *" << endl
		<< " *   USDA-ARS, Adaptive Cropping Sysems Laboratory         *" << endl
		<< " ***********************************************************" << endl
		<< endl << endl;

//Z read parameters first
//  read variety file here and fill in data structures
	//Z to hold duumy strings from variety file
	//  unimportant thing goes to BUFFER
	char* Buffer = (char*)calloc(256, sizeof(char));
	try
	{
		ifstream cfs(varietyFile, ios::in);
		if (!cfs)
		{
			throw "Variety File not found.";
		}
		cfs.getline(initInfo.description, sizeof(initInfo.description), '\n');
		string cult = initInfo.description;
		//Z bypass two description lines, 
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> initInfo.gdd_base >> initInfo.gdd_opt >> initInfo.gdd_max;
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> initInfo.fpibpa >> initInfo.srpa >> initInfo.drpa >> initInfo.tspa
			>> initInfo.iepa >> initInfo.jtpa >> initInfo.bootpa >> initInfo.headpa
			>> initInfo.antspa >> initInfo.antepa >> initInfo.matpa;;
		initInfo.GDD_rating = 1900;
		// end reading cultivar specific data from variety file
		// now read species specific data at the end of the file
		// loop until we find the location in the file
		string location = "[Gas_Exchange Species Parameters]";
		int result = -1;
		char strTest[255];
		string strTest2;
		do
		{
			//Z loop until finding the "[Gas_Exchange Species Parameters]"
			cfs.getline(Buffer, 256, '\n');
			strcpy_s(strTest, Buffer);
			strTest2.assign(strTest);
			result = strTest2.find(location);

		} while (result == -1);

		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> GasExParam.EaVp >> GasExParam.EaVc >> GasExParam.Eaj >> GasExParam.Hj >>
			GasExParam.Sj >> GasExParam.TPU25 >> GasExParam.Vcm25 >> GasExParam.Jm25 >>
			GasExParam.Rd25 >> GasExParam.Ear >> GasExParam.g0 >> GasExParam.g1;
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> GasExParam.f >> GasExParam.scatt >> GasExParam.Kc25 >> GasExParam.Ko25 >>
			GasExParam.Eac >> GasExParam.Eao;
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> GasExParam.sf >> GasExParam.phyf >> GasExParam.stomaRatio;
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs.getline(Buffer, 256, '\n');
		cfs >> GasExParam.internalCO2Ratio >> GasExParam.SC_param >> GasExParam.BLC_param;

		//Z finish reading, close the variety file
		if (cfs.eof()) cfs.close();
	}
	catch(const char* message)
	{
		cerr << message << "\n";
		exit(1);
	}

	firstDayOfSim = initInfo.beginDay;
	lastDayOfSim = initInfo.endDay;
	SowingDay = initInfo.sowingDay;

	//Z Timer class gets stepsize in hours
	dConvert.caldat(firstDayOfSim, mm, dd, yy);
	time = new Timer(dd, mm, yy, initInfo.timeStep / 60.0);

	//Z Weather class for the whole
	//  counting total records of weather data
	int dim = (int)(((lastDayOfSim + 1) - firstDayOfSim) * ((MINUTESPERDAY) / initInfo.timeStep)); 
	weather = new TWeather[dim];

	//Z initialize the rye plant here
	ryePlant = new RyePlant(initInfo, GasExParam);

//Z set heading for output (see function "output to crop file")
//  TODO may need remove the "setiosflags(ios::left)"
	{
		ofstream cropOut(cropFile, ios::out);
		cropOut << setiosflags(ios::left)
			<< setiosflags(ios::fixed)
			<< setw(12) << "date,"
			<< setw(12) << "jday,"
			<< setw(10) << "time,"
			<< setw(12) << "LeafNum,"			// number per representative plant
			<< setw(12) << "LeafArea,"			// cm2 per representative plant
			<< setw(12) << "GreenLfNum,"		// number per representative plant
			<< setw(12) << "GreenLfArea,"		// cm2 per representative plant
			<< setw(10) << "LAI,"               //Leaf Area Index
			<< setw(12) << "TillerNum,"			// number per representative plant
			<< setw(13) << "ShootMass,"			// grams per representative plant
			<< setw(14) << "RootMass,"			// grams per representative plant
			<< setw(13) << "NitroMass,"			// mg per representative plant
			<< setw(13) << "GrossPhotosyn,"		// mg biomass per representative plant per hour
			<< setw(13) << "NetPhotosyn,"		// mg biomass per representative plant
			<< setw(15) << "CumuTTd(plant),"
			<< setw(15) << "CumuTTd(climate),"
			<< setw(15) << "PltLivingFrac,"
			<< setw(15) << "BiomassReserve,"     //g per plant
			<< setw(15) << "MaintRespir (mg/p)," //mg per plant - model is grams per plant but it is multiplied by 1000
			<< setw(15) << "NitrogenAssignment,"
			<< setw(15) << "LeafNRelease,"
			<< setw(15) << "SheathNRelease,"
			<< setw(15) << "InternodeNRelease,"
			<< setw(15) << "HourlyNitrUp,"
			<< setw(15) << "NPool,"
			<< setw(15) << "Sunlit_PFD,"  //umol m-2 s-1
			<< setw(15) << "Shaded_PFD,"
			<< setw(15) << "Sunlit_A_gross,"
			<< setw(15) << "Shaded_A_gross,"
			<< setw(15) << "LeftOverC,"      //grams per plant
			<< setw(10) << "IntercepPAR,"  //intercepted PAR, umol m-2 s-1
			<< setw(9) << "Note"
			
			<< endl;
	}
}

//Z read files from 2dsoil module
//  weather data format
void RyeController::readWeatherFrom2DSOIL(const TWeather& wthr)
{
	weather[iCur] = wthr;
	//Z actual supply of water to plant, mol (water) m-2 (leaf) s-1
	ET_supply = wthr.ET_supply;
	weather[iCur].daytime = wthr.jday + wthr.time;
}

//Z main engine of the rye model
//  weather is from soil, and gas exchange parameters are read in this controller
//
int RyeController::run(const TWeather& wthr)
{
	//Z first read the weather file
	readWeatherFrom2DSOIL(wthr);

	//Z use the weather data to call the plant variable once,
	//  then process the rye growth computation.
	if (weather[iCur].jday >= initInfo.sowingDay && weather[iCur].jday <= lastDayOfSim)
	{
		RootWeightFrom2DSOIL = wthr.TotalRootWeight;
		ryePlant->RyePlantUpdate(weather[iCur]);
		MaxRootDepth = wthr.MaxRootDepth;
		AvailableWater = wthr.ThetaAvail;
		outputToCropFile();

		iCur++;
		time->step();
	}
	return 0;
}

//Z output to the crop results g01
void RyeController::outputToCropFile()
{
	int mm, id, iyyy;
	string DateForOutput;
	double av_gs = ryePlant->get_conductance();
	if (!ryePlant->get_develop()->is_emergenceStart())
	{
		av_gs = 0;
	}
	double vpd = ryePlant->get_VPD();
	if (vpd < 0)
	{
		vpd = 0;
	}
	time->caldat(weather[iCur].jday, mm, id, iyyy);
#if 0
	DateForOutput.Format("%.2d/%.2d/%4i", mm, id, iyyy);
#else
	char DateForOutputBuff[16];
	sprintf(DateForOutputBuff, "%.2d/%.2d/%4i", mm, id, iyyy);
	DateForOutput = DateForOutputBuff;
#endif
	RyeTiller* myStem= ryePlant->get_mainstem();
	string s = "";

// First handle terminal states that should override others
if (ryePlant->get_develop()->is_dead()) { 
    s = "Inactive"; 
}
else if (ryePlant->get_develop()->is_mature()) { 
    s = "Matured"; 
}
// For developmental stages that can occur together, remove the else
else {
    // Check each state independently
    if (ryePlant->get_develop()->is_startEnlongation()) { 
        s += "Elong/"; 
    }
    if (ryePlant->get_develop()->is_singleRidge()) { 
        s += "SRidge/"; 
    }
    if (ryePlant->get_develop()->is_accel()) { 
        s += "Accel"; 
    }
    // Basic states that shouldn't combine with others
    if (s.empty()) {  // Only set these if no other states were set
		if (ryePlant->get_develop()->is_emergenceEnd()) {
			s = "Emerge End";
		}
        else if (ryePlant->get_develop()->is_emergenceStart()) { 
            s = "Emerge Start"; 
        }
        else if (ryePlant->get_develop()->is_germinationEnd()) { 
            s = "Germinate End"; 
        }
        else if (ryePlant->get_develop()->is_germinationStart()) { 
            s = "Germinate start"; 
        }
        else { 
            s = "none"; 
        }
    }
}

// Optional: Remove trailing slash if present
if (!s.empty() && s.back() == '/') {
    s.pop_back();
}

	int accelerated=0, singleRidge=0, growing=0, mature = 0;
	int leafNum = ryePlant->get_leafNum();
	int i;	

	// get intercepted light - use 4.6 moles light per MJ quantum flux data
	double PAR = ryePlant->get_shaded_PFD() + ryePlant->get_sunlit_PFD();
	double cover = 1.0 - exp(-0.79 * ryePlant->get_LAI());
	ofstream ostr(cropFile, ios::app);
	ostr << setiosflags(ios::right)
		<< setiosflags(ios::fixed)
		<< setw(9) << DateForOutput << comma
		<< setw(6) << weather[iCur].jday << comma
		<< setw(8) << setprecision(0) << weather[iCur].time * 24.0 << comma
		<< setw(11) << setprecision(2) << ryePlant->get_leafNum() << comma
		<< setw(11) << setprecision(2) << ryePlant->get_leafArea() << comma
		<< setw(11) << setprecision(2) << ryePlant->get_greenLfNum() << comma
		<< setw(12) << setprecision(2) << ryePlant->get_greenLeafArea() << comma
		<< setw(6)  << setprecision(1) << ryePlant->get_LAI()     << comma
		<< setw(12) << setprecision(2) << ryePlant->get_tillerNum() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_shootMass() << comma
		<< setw(13) << setprecision(4) << RootWeightFrom2DSOIL << comma //Z it is interesting that root weight per plant from 2DSOIL is set as a component of weather
		<< setw(13) << setprecision(6) << ryePlant->get_plantTotalNitrogen() << comma
		<< setw(13) << setprecision(4) << ryePlant->get_grossPhotosynthesis()*1000. << comma
		<< setw(13) << setprecision(4) << ryePlant->get_netPhotosynthesis()*1000. << comma
		<< setw(13) << setprecision(2) << ryePlant->get_develop()->get_RyeTTd().get_sum_plant() << comma
		<< setw(13) << setprecision(2) << ryePlant->get_develop()->get_RyeTTd().get_sum_ambient() << comma
		<< setw(13) << setprecision(4) << ryePlant->get_develop()->get_plantLivingFraction() << comma
		<< setw(13) << setprecision(2) << ryePlant->get_biomassReserve() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_MaintenanceRespiration()*1000.0 << comma
		<< setw(13) << setprecision(6) << ryePlant->get_PltShootNitrogenAssignment() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_LeafNitrogenRelease() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_SheathNitrogenRelease() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_IntrNodeNitrogenRlease()<< comma
		<< setw(13) << setprecision(6) << ryePlant->get_HourlyNitrogenSoilUptake() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_nitrogenPool() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_sunlit_PFD() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_shaded_PFD() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_sunlit_A_gross() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_shaded_A_gross() << comma
		<< setw(13) << setprecision(6) << ryePlant->get_averagedBiomassLeftover() <<comma
		<< setw(13) << setprecision(3) << PAR*cover << comma
		<< setw(25) << setiosflags(ios::skipws) << "\"" + s + "\""
		<< endl;
	ostr.close();
}