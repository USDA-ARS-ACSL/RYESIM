#include "solar.h"
class CRadTrans
{
private:

	enum CLeafAngle {Spherical, Horizontal, Vertical, Diaheliotropic, Empirical, Ellipsoidal, Corn};
	//type TCover = (Glass, Acrylic, polyethyl, doublepoly, whitewashed, NoCover);
	double      absorp;		//leaf absorptivity for PAR
	double      clump ;
	double      rho_soil ;	//soil reflectivity for PAR band
	double IrradianceDirect, IrradianceDiffuse, LAI, Elev, LeafAngleFactor, KbVal, KdVal;
	CLeafAngle LeafAngle;	//Solar elevation of a beam and cumulative LAI at the layer, diffused fraction (fdf)
	bool IsLeafAngleFactorUsed;

	double Qbt(double L);	//Z total beam radiation at depth L
	double Qb(double L);	//Z unintercepted beam (direct beam) flux at depth of L within canopy
	double Qd(double L);	//Z net diffuse flux at depth of L within canopy
	double Qsoil();			//Z total PFD at the soil surface under the canopy
	void   Kb(double theta);//Z Ratio of projected area to hemi-surface area for an ellisoid 
	void   Kd(double LA);	//Z diffused light ratio to ambient, integrated over all incident angles from - 90 to 90

public:
	CRadTrans(void);		//sdf: diffused fraction of solar radiation
	~CRadTrans(void);

	void SetVal(CSolar Irradiance, double LAI, double leafAngleFactor);
	double Qsc();			//Z weighted average scattered radiation within canopy
	double Qsc(double L);	//Z scattered radiation at depth L in the canopy
	double Qtot(double L);	//Z total irradiance (dir + dif) at depth L, simple empirical approach
	double Irradiancetot(); //Z total PAR at top of the canopy
	double Qsl();			//Z mean flux density on sunlit leaves
	double Qsh();			//Z mean flux density on shaded leaves over LAI
	double Qsl(double L);	//Z flux density on sunlit leaves at delpth L
	double Qsh(double L);	//Z diffuse flux density on shaded leaves at depth L
	double Qdm() ;			//Z weighted average absorved diffuse flux over depth of L within canopy accounting for exponential decay
	double Qsoilm();		//Z weighted average of Soil reflectance over canopy accounting for exponential decay
	double GetZenith();
	double Reflect();
	double LAIsl();			//Z sunlit LAI assuming closed canopy; thus not accurate for row or isolated canopy
	double LAIsh();			//Z shaded LAI assuming closed canopy
	double Fsl(double L);	//Z sunlit fraction of current layer
	double Fsh(double L);	//Z shaded fraction of current layer, "1.0 - sunlit"


	double GetKb(){return KbVal;}    // extiction coefficient assuming spherical leaf dist
	double GetKd(){return KdVal;}

	double get_RadAbsorbTot();	//Z this is for testing, total PAR absorption (umol m^-2 ground sec^-1), but could provide useful information
};
