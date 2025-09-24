/* unit Solar;
{Unit to calculate solar geometry including solar elevation, declination,
azimuth etc using TSolar class. Data are hidden. 03/15/00 SK
- 1st Revision 10/10/00: Changed to dealing only upto the top of the canopy. Radiation transter within the canopy is now a separate module.
- added functins to calculate global radiation, atmospheric emissivity, etc as in Spitters et al. (1986), 3/18/01, SK
24Dec03, SK
- added a function to calculate day length based on Campbell and Norman (1998) p 170,
- added overloads to SetVal
- added direct and diffuse separation functions according to Weiss and Norman (1985), 3/16/05
Translated from pascal to c++ by Dennis Timlin, October, 2008
CSolar  destructor - doesn't pass any values but consider updating to pass location constants
Such as lattitude etc

Note on measurements and measured values:
Typically measurements are made of PAR (visible light) using sensors such as light bars from Li-Cor or Apogee.
These are measured in units of Watts m-2 or quanta of photons m-2 s-2. Pyranometers provide measurements
of the total light spectra. In this case NIR is inferred from the two values. The model here assumes all values
are in Watts m-2. If PAR is available, PFD (photon flux density) is calculated directly using cPFD (4.6 - see
Solar.h). If only PAR is available, the total solar radiation (PAR+NIR) is calculated using ratios of potential
PAR/NIR. This ratio is usually near 0.5. In the case where PAR and Pyronometer data are both available the
ratio can be calculated directly. For the purose of the model, PFD is always input as PAR and calculated.

  this may have to change to utilize multiple days
  SetVal int Day, double Time, double Lat,  double Longi, double tauin, double SolRad0,
					double PFD0 ,double Alti
  CalcSolarNoon - calculates solar noon - no arguments
  CalcDaylength
  Sunrise  uses solarnoon
  Sunset   uses solarnoon
  CalcDeclination = calculates declinatino of sun
  CalcElev - calculates the elevation of the sun from the horizon in radians
    (from CalcSinElevation)
  CalcSinElevation - calculates the sine of the elevation angle from a formula based
    on declination and time of day
    CSolar::CosTheta() the same as Sin_Elev, but theta is the zenith angle (solar elevation
	 measured from zenith or top of sky)
	 Cos_Elev is the cosine of solar elevation (from the horizon)
	 CalcAzimuth The solar azimuth angle is the angular distance between due South and the
projection of the line of sight to the sun on the ground.
STotal() Total solar radiation SDiffuse() + SDirect(); //eqn 11.9 Campbell and Norman, 1998
SDiffuse() diffuse  and  SDirect() direct components of solar radiation
press() // atmospheric pressure in kPa
m() // the optical air mass number Campbell and Norman, pg 173
setPARFraction() // calculate PAR fraction
PotentialPARTotal() - total potential PAR
PotentialNIRTotal() total of par dif and nir
setPARFraction() set PAR fraction
setPARFraction(double fr)use when both PFD and SolRad measurements are available
PotentialPARDirect()
SolRad includes both PAR and NIR in Wm-2.
PFD is in flux density which will vary between 4.2-4.6 umol/J depending on cloudiness.

Values of direct and diffuse PAR and NIR are calculated but only used for calculating fractions of these components
To be consistent, total solar radiation is calculated one way and then fractioned into diffuse, direct, par and nir
using fractions.

*/
#include <cmath>
#include <algorithm>
#include "solar.h"

using namespace std;

inline double DegToRad(double Deg) { return (Deg/180.0)*PI; }
inline double RadToDeg(double Rad) { return (Rad*180.0)/PI; }

//// **** Initialization of the CSolar class ***** ////

double CSolar::SDERP[9] = {
		0.3964E-0,0.3631E1,0.3838E-1,0.7659E-1,0.0000E0,-0.2297E2,-0.3885E0,-0.1587E-0,-0.1021E-1
	};
CSolar::CSolar()
{
	//Z make initialization row by row for the private variables
	JDay = 0;
	Time = Latitude = Longitude = altitude = 0.0, tau = 0.75;
	DayLength = 12.0, Declination = SolarNoon = SinElevation = Azimuth = 0.0;
	Sunrise = Sunset = 0.0, HalfDay = 6.0, Elevation = CosElevation = 0.0;
	CosTheta = 0.0;

	PAR = -99, SolarRadiation = -99, PotentialSolarTotal = PotentialSolarDirect = PotentialSolarDiffuse = 0.0;
	NIRTotal = NIRDiffuse = NIRDirect = PARTotal = PARDiffuse = PARDirect = 0.0;
	PFD = -99, NIR = 0.0, PARFraction = 0.5, NIRFraction = 0.0;
	PotentialPARDiffuse = PotentialPARDirect = PotentialPARTotal = 0.0;
	PotentialNIRDiffuse = PotentialNIRDirect = PotentialNIRTotal = 0.0;
	SolarDirect = SolarDiffuse = FracDiffuse = PFDDirect = PFDDiffuse = 0.0;
	FracPARDirect = FracPARDiffuse = FracNIRDirect = FracNIRDiffuse = FracSolarDiffuse = FracSolarDirect = 0.0;
	FracPFDDirect = FracPFDDiffuse = 0.0;
	RowAzimuth=0.0;
	
	useObs = false; useTau = false;

}
CSolar::~CSolar()
{
}

// Altitude available
void CSolar::SetVal(int Day, double Time, double Lat, double Longi, double Alti, double SolRad0)
{
	// We assume only radiation is available (MJ) PAR can be calculated from radiation
	// usetau - flag to use measured tau for calculation of the PAR fraction.
	// here we use the potential values of PAR and NIR to calculate the PAR fraction
	JDay = Day;
	this->Time = Time * 24.0;
	Latitude = DegToRad(Lat);		//Z Latitude in rad
	Longitude = Longi;				//Z Longitude in deg
	useTau = false;
	useObs = true;
	altitude = max(50.0,Alti);	//Z Assume minimum altitude is 50 m from the sea level
	SetDeclination();				//Z declination in rad
	SetDayLength();					//Z day length in hour
	SetSolarNoon();					//Z solar noon in hour
	SetSolarElevation();
	SetAzimuth();
// set up potential radiation values
	SetPotentialSolar();			//Z Potential solar						W m^-2
    SetPotentialPAR();				//Z Potential solar photosynthesis rad  W m^-2
	SetPotentialNIR();				//Z Potential solar near infrared rad   W m^-2

//by default consider both solar radiation and PFD available
//using the parfraction to partition solar rad is reasonable, see Wiess and Norman, 1985, Eq 7,8
//   https://www.sciencedirect.com/science/article/pii/0168192385900206
	SetPARFraction();				//Z uses tau to set PAR fraction (NIR and PAR)
	SolarRadiation = SolRad0;
	PAR = SolarRadiation * PARFraction;
	NIR = SolarRadiation * (1.0 - PARFraction);

//Z from here, we have any unnessary calling, circular argument, "NIR_total+PAR_total = solarradiation"?
//  temporarily remove these two lines for fast computation
//	SetPARFraction(PAR/SolarRadiation);
//	SetPotentialNIR();

//Z for setting PAR or NIR, please follow this order
//  since the total amount comes first and then partitioning can be done

	SetFracNIRDirect();				//Z Fraction of NIR in direct beam
	SetNIRTotal();
	SetNIRDirect();
	SetNIRDiffuse();

	SetFracPARDirect();				//Z Fraction of PAR in direct beam
    SetPARTotal();
    SetPARDirect();
	SetPARDiffuse();

//	PFD = PAR * cPFD;
// get light components (diffuse and direct)
}


//// ***** solar geometry calculations ***** ////

//Z Declination (delta) in astronomy is comparable to geographic latitude, i.e., the latitude of straight solar rad
//  From GLYCIM
void CSolar::SetDeclination()
{
	double Ang1;
	int i, n, j;

    Declination = SDERP[0];
	for (i = 2; i <= 5; i++)
	{
          n = i - 1;
          j = i + 4;
		  Ang1 = n * 0.01721 * JDay;
		  Declination = Declination + SDERP[i - 1] * sin(Ang1) + SDERP[j - 1] * cos(Ang1);
	}
    Declination = DegToRad(Declination);
}

void CSolar::SetDayLength()
{
	/*from Jeffrey Amthor, Calculation of Daylength. Bioinformatics, 13:479-480
       Also see campbell and Norman, 1998, chapter 11, eqn 11.6
	   Daylength(N)= 24*hs/PI where hs varies from 0 to pi and is the hour angle of the sun (angular distance from the
	   meridian of a site in radians) at sunset
	   hs is obtained by rearranging the geometric equation for sin solar elevation

	*/
	double LatRad, D1, D2, D3, hs;
	LatRad = Latitude;

	D1 = sin(LatRad) * sin(Declination);
	D2 = cos(LatRad) * cos(Declination);
	D3 = D1 + D2;

	hs = acos((-0.014544 - D1) / D2);	//- 0.14544 is a small negative number to indicate twilight is considered
	DayLength = hs * 24.0 / PI;			//Z hs=\pi/2=90deg means 12 hour (geometrical) day length
	HalfDay = DayLength / 2.0;
}

void CSolar::SetSolarNoon()
// only works now for positive longitude
// DT-March 1 - 2013 changed program to be universally applicable
 {
	double TimeZone;
	double divisor = 15.0;
	double LC, EqTime, Epsil; 

//Z --------- old code I ----------------------
/*
// LC is longitude correction for solar noon, Wohlfart et al, 2000; Campbell & Norman 1998
// find offset in degrees from meridian of standard time
//  LC is a correction for each degree east of a standard meridian
	LC = (75-Longitude)*1/15;  // standard meridIn for pacific time zone is 120 W, Eastern Time zone : 75W

// this should be more general
	LC=fmod(Longitude,15);
	LC=(LC*-1.0)*1/15;; //+4 minutes for each degree east of the standard meridian (0, 15 30, etc) and -4 minutes for each degree west
	// east longitude is positive, west longitude is negative. North latitude is positive
	Epsil = DegToRad(279.575 + 0.9856*JDay);
	EqTime = (-104.7*sin(Epsil)+596.2*sin(2*Epsil)+4.3*sin(3*Epsil)
		-12.7*sin(4*Epsil) - 429.3*cos(Epsil)-2.0*cos(2*Epsil)
		+ 19.3*cos(3*Epsil))/3600;  // Calculating Equation of Time Eq 11.4 in Campbell and Norman (1998)
	// fixed errors in original version (some components adjusted for radians twice)

	// EqTime = 0; //Ignore for computing speed
	SolarNoon=12.0 - LC - EqTime;
	Sunrise=SolarNoon - HalfDay;
	Sunset=SolarNoon + HalfDay;
*/

//Z --------- old code II ----------------------
/*	LC=fmod(Longitude,divisor);
	LC=(LC*-1.0) / divisor; //+4 minutes (1/15 hour) for each degree east of the standard meridian (0, 15 30, etc) and -4 minutes for each degree west
	// east longitude is positive, west longitude is negative. North latitude is positive
	Epsil = DegToRad(279.575 + 0.9856*JDay);
	EqTime = (-104.7*sin(Epsil)+596.2*sin(2*Epsil)+4.3*sin(3*Epsil)
		-12.7*sin(4*Epsil) - 429.3*cos(Epsil)-2.0*cos(2*Epsil)
		+ 19.3*cos(3*Epsil))/3600;  // Calculating Equation of Time Eq 11.4 in Campbell and Norman (1998)
	// fixed errors in original version (some components adjusted for radians twice)

	// EqTime = 0; //Ignore for computing speed
	SolarNoon=12.0 - LC - EqTime;
	Sunrise=SolarNoon - HalfDay;
	Sunset=SolarNoon + HalfDay;
*/

	TimeZone = floor((Longitude + 7.5) / divisor);	//Z determine the time zone using floor function
	LC = Longitude - TimeZone * divisor;			//Z residue longtitude within that time zone, positive (east, early) or negative (west, late)
	LC = LC / divisor;								//Z residue longtitude -> time, "-" indicate earlier or later

	// Calculating "Equation of Time" Eq 11.4 in Campbell and Norman (1998)
	Epsil = DegToRad(279.575 + 0.9856 * JDay);
	EqTime = (-104.7 * sin(Epsil) + 596.2 * sin(2 * Epsil) + 4.3 * sin(3 * Epsil)
		- 12.7 * sin(4 * Epsil) - 429.3 * cos(Epsil) - 2.0 * cos(2 * Epsil)
		+ 19.3 * cos(3 * Epsil)) / 3600;

	SolarNoon = 12.0 - LC - EqTime;

	//Z because the daylength includes the twilight, so sunrise/sunset includes that correction
	Sunrise = SolarNoon - HalfDay;
	Sunset = SolarNoon + HalfDay;
}

//Z solar altitude angle, the angle of the sun relative to the earth's horizon
//  We compute elevation here, while zenith (theta) is "\pi/2-elevation"
//  for beam radiation:
//		Elevation: angle with earth surface
//		Theta: angle with vertical line (zenith angle)
void CSolar::SetSolarElevation()
{
	SinElevation = sin(Latitude) * sin(Declination) + cos(Latitude) * cos(Declination) * cos(DegToRad(15.0 * (Time - SolarNoon)));
	Elevation = asin(SinElevation);
	//Z macro "FDIV_GUARD" is a small number that prevent the result from being 0.0
	Elevation = max(FDIV_GUARD, Elevation);
	// identical to sin_elev, theta is the zenith angle
	CosTheta=max(FDIV_GUARD, SinElevation);	//Z CosTheta=SinElevation, theta is the zenith angle, and Theta+Elevation=pi/2
	//if (fabs(CosTheta) < FDIV_GUARD)
	//{
	//	CosTheta = 0.0;
	//}
	CosElevation = sqrt(1.0 - SinElevation * SinElevation);
}


/* {Obtained from http://www.susdesign.com/sunangle/index.html
The solar azimuth angle is the angular distance between due South and the
projection of the line of sight to the sun on the ground.
View point from south, morning: -, afternoon: +} (+ is clockwise from south)
This is similar to equations for this caclulation but time of day is left out
See eqn 11.1 in campbell and norman, 1998. Calculations are also in glycim. this should
be updated for time of day and other dependent equations modified.
*/

//Z south is the center 0 rad
//  clockwise, so east is -90deg or -1.57 rad
//                west is 90deg or 1.57 rad
void CSolar::SetAzimuth()
{
	if (Time < SolarNoon)
	{
		Azimuth = -acos((SinElevation * sin(Latitude) - sin(Declination)) / (cos(Latitude) * CosElevation));
	}
	else if (Time > SolarNoon)
	{
		Azimuth = acos((SinElevation * sin(Latitude) - sin(Declination)) / (cos(Latitude) * CosElevation));
	}
	else Azimuth = 0.0;
}

//Z atmospheric pressure in kPa
//  see example 11.2 in Normann and Campbell 1998
double CSolar::press() 
{
	return 101.325 * exp(-altitude / 8200.0); // campbell and Norman (1998), p 41
}

//Z the optical air mass number
//  see example 11.2 in Normann and Campbell 1998
double CSolar::m() 
{
	return press() / (101.325 * CosTheta); // campbell and Norman (1998), p 173
}

//// **** solar radiation calculations for total and components
///        potential values are those calculated for the conditions
///    of no cloud cover and value of tau - transmissivity***** ////

//Z STotal is total solar radiation, W m^-2 (based on the solar constant unit)
void CSolar::SetPotentialSolar()
{
	if (SinElevation >= 0.0)
	{
		// note that the adjustment factor (1+0.033*Cos...) is a very small number that varies from 0.97 to 1.03 with a larger magnitude in the winter.
		// This is a correction for the elliptical shape of the earth that results in a slighlty different path of light
		// 
		// Uses solar elevation rather than zenith
        // so sin function is appropriate
        // compare to eq 11.11 in campbell norman book.

		//Z D11 is just a temperoary value
		double Spo = SolarConst * SinElevation * (1.0 + 0.033 * cos(2.0 * PI * (JDay - 10.0) / 365.0));
		PotentialSolarDirect = pow(tau, m()) * Spo;									// tau: atmospheric transmissivity, HGJones 1991
		PotentialSolarDiffuse = 0.3 * (1.0 - pow(tau, m())) * Spo * CosTheta;		// campbell and Norman (1998), 11.13
		//campbell and norman eqns (11.8 and 11.11)
	}
	else
	{
		PotentialSolarDirect = 0.0;
		PotentialSolarDiffuse = 0.0;
	}

	PotentialSolarTotal = PotentialSolarDiffuse + PotentialSolarDirect; //eqn 11.9 Campbell and Norman, 1998
}

//(beam radiation - direct radiation on a surface perpendicular to the beam)

/////***** PFD calculations ******/////

//Z potential PAR calculations (visible light), unit in W m^-2
//  According to Weiss and Norman (1985)
void CSolar::SetPotentialPAR()
{
	PotentialPARDirect = 600.0 * exp(-0.185 * m()) * CosTheta;				//Z Eq.1 in Weiss and Norman (1985) 
	PotentialPARDiffuse = 0.4 * (600.0 - PotentialPARDirect) * CosTheta;	//Z Eq.3 in Weiss and Norman (1985) 
	PotentialPARTotal = PotentialPARDirect + PotentialPARDiffuse;
}

// //***** NIR calculations ***** /////
//Z potential total Near-Infrared Light (NIR, W m-2) high wavelength energy
//  From Weiss and Norman, 1985
void CSolar::SetPotentialNIR()
{
	// potential diffuse, direct and total NIR (W m-2)

	double x;
	double w;			//Z water absorption in the near infrared for 10 mm of precipitable water
	x = 1 / CosTheta;	//Z optical airmass (m) as in Weiss and Norman (1985),
						// note: m in Campbell and Norman (1989) includes pressure correction, see above
	w = -1.1950 + 0.4459 * log10(x) - 0.0345 * log10(x) * log10(x);
	w = 1320.0 * pow(10.0, w);		//water apsorption, 1320 is solar constant excluding UV

	PotentialNIRDirect = (720.0 * exp(-0.06 * m()) - w) * CosTheta; // m here is actually (P/P0)*m in Weiss and Norman
	PotentialNIRDiffuse = 0.6 * (720.0 - PotentialNIRDirect - w) * CosTheta;
	PotentialNIRTotal = PotentialNIRDirect + PotentialNIRDiffuse;
}

// calculate PAR fraction
void CSolar::SetPARFraction()
{
	double tmp;
	tmp = PotentialPARTotal / (PotentialPARTotal + PotentialNIRTotal);
	if (useTau)
	{
		if (tau >= 0.7) { tmp = 0.45; }         // clear sky (tau >= 0.7): 45% is PAR      Goudriaan and van Laar (1994)
		else
		{
			if (tau <= 0.3) { tmp = 0.55; }		// cloudy sky (<= 0.3): 55% is PAR
			else { tmp = 0.625 - tau * 0.25; }
		}
	}
	PARFraction = min(max(FDIV_GUARD, tmp), (1.0 - FDIV_GUARD));
}

//use when both PFD and SolRad measurements are available
void CSolar::SetPARFraction(double Fraction)
{
	PARFraction = min(max(FDIV_GUARD, Fraction), (1.0 - FDIV_GUARD));
}

void CSolar::SetFracPARDirect() //Fraction of PAR in direct beam
{
	const double A = 0.9, B = 0.7;
	double ratio;
//	double ratio2, R1, R2, R3;
	//Z SolarRadiation is the measured one, while potential is computed one
	ratio = SolarRadiation / (PotentialPARTotal + PotentialNIRTotal);
//	ratio2=A-FDIV_GUARD;
//	R1=A-min(A-FDIV_GUARD,ratio);
//	R2=R1/B;
//	R3=pow(R2,2.0/3.0);
	FracPARDirect = max(0.0, (PotentialPARDirect / PotentialPARTotal)
		* (1 - pow((A - min(A - FDIV_GUARD, ratio)) / B, (2.0 / 3.0))));
}

void CSolar::SetPARTotal() // for measured
{
	//PARTotal= SolarRadiation*GetPARFraction();
	PARTotal = SolarRadiation * PARFraction;
}

void CSolar::SetPARDirect() // for measured
{
	//
	PARDirect = PARTotal * FracPARDirect;
}

void CSolar::SetPARDiffuse()
{
	//PARDiffuse= PARTotal*GetFracPARDiffuse();
	FracPARDiffuse = 1.0 - FracPARDirect;
	PARDiffuse = PARTotal * FracPARDiffuse;
}


void CSolar::SetFracNIRDirect() //Fraction of NIR in direct beam
{
	const double C = 0.88, D = 0.68;
	double ratio;
	ratio = SolarRadiation / (PotentialPARTotal + PotentialNIRTotal);
	FracNIRDirect = max(0.0, (PotentialNIRDirect / PotentialNIRTotal) * (1.0 - pow((C - min(C - FDIV_GUARD, ratio)) / D, (2.0 / 3.0))));
}

void CSolar::SetNIRTotal() // for measured
{
	//NIRTotal= SolarRadiation*GetNIRFraction();
	NIRFraction= 1.0 - PARFraction;
	NIRTotal = SolarRadiation * NIRFraction;
}

void CSolar::SetNIRDirect() // for measured
{
	NIRDirect = NIRTotal * FracNIRDirect;
}

void CSolar::SetNIRDiffuse()
{
	//NIRDiffuse= NIRTotal*GetFracNIRDiffuse();
	FracNIRDiffuse = 1.0 - FracNIRDirect;
	NIRDiffuse = NIRTotal * FracNIRDiffuse;
}