/*unit CanopyRadTrans;
{Basic canopy architecture parameters, 10/10/00 S.Kim
 modified to represent heterogeneous canopies
 Uniform continuouse canopy: Width2 = 0
 Hedgerow canopy : Width1 = row width, Width2 = interrow, Height1 = hedge height, Height2 = 0
 Intercropping canopy: Height1, Width1, LA1 for Crop1, and so on
 Rose bent canopy: Height1=Upright canopy, Height2 = bent portion height, 10/16/02 S.Kim}

 This is the interface into the Solar Class but calculates transmission coefficients. To use it you must
 initialize the Solar object. See Solar.cpp for more information.

 This routine adds calculations that require use of LAI and LAF (Leaf Angle Factor)
 You must pass parameters to the Solar Class to initialize that class and use the methods therin.
 The LeafAngleFactor can be numeric in which case it is used directly in one equation to calculate
 Kd as a function of zenith angle (eqn 15.4 in Campbell and Norman). In this case use 'IsLeafAngleFactorUsed' as true
 Otherwise you can use shapes as an enumeration.
*/
// dateutils ;

#include <cmath>  // need for math functions
#include <algorithm>  // need for max min functions
#include "radtrans.h"

using namespace std;

inline double cot(double a) { return 1/tan(a); }
inline double sqr(double a) { return (a * a); }


// include Math unit for Real mathmatical doubles

CRadTrans::CRadTrans()
{
//Z LeafAngle needs to be updated, while the "IsLeafAngleFactorUsed" indicates it is not used here
//Z still recommend to be updated, since is updated later.
//  directly use "LeafAngleFactor", see Table 15.1 in Campbell and Normann 1998, see "ryegrass"
	LeafAngle = Spherical;
	LeafAngleFactor = 1.035; //Z average leaf angle factor for rye grass
	IsLeafAngleFactorUsed = false;
	absorp = 0.85;   //leaf absorptivity for PAR
	clump = 1.0;
	rho_soil = 0.10; // soil reflectivity for PAR band

//Z get other variables initial values, such that
	IrradianceDirect = IrradianceDiffuse = LAI = Elev = LeafAngleFactor = KbVal = KdVal = 0.0;

}

CRadTrans::~CRadTrans()
{
}


 void CRadTrans::SetVal(CSolar Irradiance, double SLAI, double LeafAngleFactorIn)
 {

	 IrradianceDirect = Irradiance.GetPFDDirect();		//Z PFD (umol m-2 s) direct from solar irradiance
	 IrradianceDiffuse = Irradiance.GetPFDDiffuse();	//Z PFD (umol m-2 s) diffuse from solar irradiance
// the transmittance values obtained from Day and Bailey (1999), chap 3, 
// ecosystems of the world 20: greenhouse ecosystem, page 76
     LAI = SLAI;
     Elev =  Irradiance.GetSolarElevation();			//Z solar elevation angle in rad
     LeafAngleFactor = LeafAngleFactorIn;				//Z input leaf angle factor
     IsLeafAngleFactorUsed = true;
	 Kb(GetZenith());					//Z KbVal: light interception by canopies, beam radiation extinction coefficient 
	 Kd(LAI);							//Z KdVal: extinction coefficient for black leaf in diffusive radiation

//  GDiffuse = S;
}

// need to move this to CSolar
//Z can be moved to CSolar but not necessary: in CSolar, use "Theta" angle as Zenith
//  Zenith in rad
double CRadTrans::GetZenith()  
{
	double zenith;
	// need to constrain zenith to not quite reach PI/2 when elevation is zero
	// i.e., the sun is near the horizon.
	zenith = abs(PI / 2.0 - Elev);
	zenith = min(zenith, 1.56);
	return zenith;
}

double CRadTrans::Reflect()
 {
	//Z for private variable, no need to use get function, can use this->
	// return (1-sqrt(absorp))/(1+sqrt(absorp))*(2*GetKb()/(GetKb()+GetKd()));
	return (1 - sqrt(absorp)) / (1 + sqrt(absorp)) * (2 * KbVal / (KbVal + KdVal));
 }


//Z Kb (a ratio for radiation)
//  Campbell, p 251, Eq. 15.4, Ratio of projected area to hemi-surface area for an ellisoid
//  x is a leaf angle distribution parameter
void CRadTrans::Kb(double theta) 
{
	double x, tmp;
	tmp = 0.5;

/*Z seems this section will be overwritten by the following section with "IsLeafAngleFactorUsed == true"
	switch (LeafAngle)
	{
	case Spherical:
		tmp = 0.5/sin(Elev); // (0.5/sin_elev); When Lt accounts for total path length, division by sin(elev) isn't necessary
		break;
	case Horizontal:
		tmp = 1/sin(Elev); //1
		break;
	case Vertical:
		tmp = (2/cot(Elev))/PI;
		break;
	case Diaheliotropic:
		tmp = 1; // 1/sin_elev;
		break;
	case Empirical:
		tmp = 0.667; //(0.5)*KDiffuse/(0.8*sqrt(1-scatt)); //(0.5/sin_elev)*KDiffuse/(0.8*sqrt(1-scatt));
		break;
	default:
		tmp= 0.5/sin(Elev);
	}
*/

	if (IsLeafAngleFactorUsed == true) { x = LeafAngleFactor; }
	else
	{
		switch (LeafAngle)
		{
		case Spherical:
			x = 1; //
			break;
		case Horizontal:
			x = 10; //1
			break;
		case Vertical:
			x = 0;
			break;
		case Corn:
			x = 1.37;
			break;
		default:
			x =1;
		}
	}
	//   if Sin(Elev) > sin(5) then
	//   tmp =  0.5/sin(Elev)

	tmp = sqrt(sqr(x) + sqr(tan(theta))) / (x + 1.774 * pow(x + 1.182, -0.733));

	//   else tmp = 0.5/sin(5);

	KbVal = tmp * clump;
}

//Z Kd, extinction coefficient 
//  for black leaf in diffusive radiation
void CRadTrans::Kd(double LA)
{
	const double gauss3[3] = { -0.774597,0,0.774597 };			// abscissas
	const double weight3[3] = { 0.555556,0.888889,0.555556 };

	double K, FDiffuse, tmp, angle, x;
	if (IsLeafAngleFactorUsed == true) { x = LeafAngleFactor; }
	else
	{
		switch (LeafAngle)
		{
		case Spherical:
			x = 1; //
			break;
		case Horizontal:
			x = 10; //1
			break;
		case Vertical:
			x = 0;
			break;
		case Corn:
			x = 1.37;
			break;
		default:
			x = 1;
		}
	}	//end else
	
	FDiffuse = 0.0;
	for (int ii = 0; ii < 3; ii++)  //diffused light ratio to ambient, itegrated over all incident angles from -90 to 90
	{
		angle = (PI / 2.0) / 2.0 * (gauss3[ii]) + (PI / 2.0) / 2.0;
		tmp = sqrt(sqr(x) + sqr(tan(angle))) / (x + 1.774 * pow(x + 1.182, -0.733));
		FDiffuse = FDiffuse + (PI / 2.0) / 2.0 * (2.0 * exp(-tmp * LA) * sin(angle) * cos(angle)) * weight3[ii];
	}
	if (LA <= 0.0) 
	{
		K = 0.0;
	}
	else { K = -log(FDiffuse) / LA; }
	KdVal = K * clump;
  }


double CRadTrans::Irradiancetot() //total irradiance at the top of the canopy, passed over from either observed PAR or TSolar or TIrradiance
{
	//IrradianceDirect: beam radiation at top of canopy, IrradianceDiffuse: diffuse radiation at top.
	// 
	//Z To be exact:
	//Z IrradianceDirect = Irradiance.GetPFDDirect()	 PFD (umol m-2 s) direct from solar irradiance
	//Z IrradianceDiffuse = Irradiance.GetPFDDiffuse() 	 PFD (umol m-2 s) diffuse from solar irradiance
	
	return (IrradianceDirect + IrradianceDiffuse);  
}

double CRadTrans::Qtot(double L) //total irradiance (Direct + Diffuse) at depth L, simple empirical approach
{
	//return Irradiancetot() * exp(-sqrt(absorp) * (KdVal + KbVal) / 2.0 * L);

	//Z separate the direct and diffusive radiation and then sum together
	//  refer to Campbell and Norman, 1998 P258 Eqs 15.15-17
	return IrradianceDirect * exp(-sqrt(absorp) * KbVal * L) + IrradianceDiffuse * exp(-sqrt(absorp) * KdVal * L);
}

double CRadTrans::Qbt(double L) // total beam radiation at depth L
{
	return IrradianceDirect * exp(-sqrt(absorp) * KbVal * L);
}

double CRadTrans::Qd(double L) // net diffuse flux at depth of L within canopy
{
	return IrradianceDiffuse * exp(-sqrt(absorp) * KdVal * L);
}

//Z weighted average absorved diffuse flux over depth of L within canopy accounting for exponential decay
//  the return value is the averaged photon flux density on an "one-LAI" representative surface

//Z For example, suppose your LAI = 10, and you think that it is a 10-layer structure (each layer has LAI=1). 
//  Choose 1 layer (LAI=1) to represent that 10-layer canopy, 
//  then those equations are calculating the flux density received on that one representative layer with LAI=1.
//  so there is no difference between leaf and ground area numerically.

//Z the computed value can be used as leaf-area based computations, by multiplying correct LAI fractions
//  Campbel and Normal called it "leaf-hemi-surface area". 
double CRadTrans::Qdm() 
{
	if (LAI <= 0.0) { return 0.0; }
	else {
		//  Integral Qd / Integral L
		//Z see Campbell and Norman section 15.9 P261 for an example
		return IrradianceDiffuse * (1.0 - exp(-sqrt(absorp) * KdVal * LAI)) / (sqrt(absorp) * KdVal * LAI);
	}
 }

double CRadTrans::Qb(double L) // unintercepted beam (direct beam) flux at depth of L within canopy
{
	return IrradianceDirect * exp(-KbVal * L);
}

//Z mean flux density on sunlit leaves: 
//  two components, one is direct PFD with extinction factor, one is diffusive PFD
//  because "sunlit", use GetKb() or KbVal is sufficient

//Z add absorp=0.85 as the leaf absorped photon flux density
//  Qsh() should already have absorp=0.85 in its own function
double CRadTrans::Qsl()
{
	return absorp * KbVal * IrradianceDirect + Qsh();
}

//Z add absorp=0.85 as the leaf absorped photon flux density
//  Qsh(L) should already have absorp=0.85 in its own function
double CRadTrans::Qsl(double L) // flux density on sunlit leaves at depth L (depth is marked via LAI incremental)
{
	return absorp * KbVal * IrradianceDirect + Qsh(L);
}

//Z mean flux density on shaded leaves over LAI, 
//  just one diffusive PFD component (but has three parts)

//Z add absorp=0.85 as the leaf absorped photon flux density
double CRadTrans::Qsh()
{
	return absorp * (Qdm() + Qsc() + Qsoilm());  // include soil reflection
}

double CRadTrans::Qsh(double L) // diffuse flux density on shaded leaves at depth L
{
	return absorp * (Qd(L) + Qsc(L) + Qsoilm());   // include soil reflection
}

// weighted average of Soil reflectance over canopy accounting for exponential decay

//Z the computed value can be used as leaf-area based computations, by multiplying correct LAI fractions
//  Campbel and Normal called it "leaf-hemi-surface area". 
//  That means it is the average PFD value on a representative, LAI=1, intersection surface
//  similar to Qdm()
double CRadTrans::Qsoilm() 
{
	if (LAI <= 0.0) { return 0.0; }
	else {
		//  Integral Qd / Integral L
		//Z to understand this equation
		//   "Qsoil()" is the radiation reach to soil surface
		//   "Qsoil() * rho_soil" is the radiation reflected by soil surface
		//   "whole equation" is the absorption of the reflected equation took by the canopy when the radiation is transporting upwards
		return Qsoil() * rho_soil * (1.0 - exp(-sqrt(absorp) * KdVal * LAI)) / (sqrt(absorp) * KdVal * LAI); 
	}
 }

// weighted average scattered radiation within canopy
double CRadTrans::Qsc() 
{
	double totBeam, nonscatt;
	if (LAI == 0.0) { return 0.0; }
	else
	{
		//total beam including scattered absorbed by canopy
		totBeam = IrradianceDirect * (1.0 - exp(-sqrt(absorp) * KbVal * LAI)) / (sqrt(absorp) * KbVal);
		//non scattered beam absorbed by canopy
		nonscatt = IrradianceDirect * (1.0 - exp(-KbVal * LAI)) / (KbVal);
		//mean scattered flux density
		return (totBeam - nonscatt) / LAI; 
	}
}

double CRadTrans::Qsc(double L) // scattered radiation at depth L in the canopy
{
	return Qbt(L) - Qb(L); // total beam - nonscattered beam at depth L
}

double CRadTrans::Qsoil() // total PFD at the soil sufrace under the canopy
{
	return Qtot(LAI);
}

double CRadTrans::LAIsl() // sunlit LAI assuming closed canopy; thus not accurate for row or isolated canopy
{
	if (Elev <= 0.01) { return 0.0; }
	else { return (1.0 - exp(-KbVal * LAI)) / KbVal; }
}

double CRadTrans::LAIsh()// shaded LAI assuming closed canopy
{
    return LAI - LAIsl();
}

// sunlit fraction of current layer
double CRadTrans::Fsl(double L) 
{
	if (Elev <= 0.01) { return 0.0; }
	else { return exp(-KbVal * L); }
}

// shaded fraction of current layer
double CRadTrans::Fsh(double L)
{
	return 1 - Fsl(L);
}

//Z this is for testing, total PAR absorption (umol m^-2 ground sec^-1), but could provide useful information
//  total radiation absorbed by the plant canopy
double CRadTrans::get_RadAbsorbTot()
{
	if (LAI <= 0.0)
	{
		return 0.0;
	}
	else
	{
		return (Qsl() * LAIsl() + Qsh() * LAIsh());
	}
}