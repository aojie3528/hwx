/***************************************************************************

 Orbit computation using two body theory


 Reference: IAU SOFA EPV00 for Model 1
  
***************************************************************************/

#ifndef TWOBODY_H
#define TWOBODY_H

#include <windows.h>
#include <math.h>

#include "constant.h"

const double dfEarthGMIOD = sqrt( g_dfEarthGM / pow( g_dfEarthSemiMajor, 3.0 ) * 3600.0 );
const double dfMuIOD = 1.0;

struct TwoBodyOrbitElement
{
    double semi_major;				// meters
    double eccentricity;
    double inclination;				// radians
    double perigee;					// radians
    double Longi_RisingNode;		// radians
    double MeanAnomaly;				// radians
    double referencetime_jd;		// julian date

	double referencetime_jd_int;
	double referencetime_jd_fra;

	int GM;							// 1 semi-major in meters
									// 0 semi-major in earth radius
};

class cTwoBody
{
    TwoBodyOrbitElement m_stSatellite;
    double m_dfX, m_dfY, m_dfZ, m_dfXdot, m_dfYdot, m_dfZdot;

	double m_dfTransCoe;

public:
    cTwoBody();
    ~cTwoBody();

    void SetElement( TwoBodyOrbitElement &stSatellite );
    void ComputeInertialXYZ( double dfTimeJD );
    void ComputeInertialXYZ( double dfIntJD, double dfFraJD );
    void ComputeECEFXYZ( double dfTimeJD );
	void TransformInertialToECEF( double dfTimeJD );
    void GetXYZ( double &dfX, double &dfY, double &dfZ );
    void GetXYZ( double *pdfPos );
    void GetXYZdot( double &dfXdot, double &dfYdot, double &dfZdot );
    void GetXYZdot( double *pdfVel );	
	double GetSemiMajor();
	double GetInclination();
	int GetGM() { return m_stSatellite.GM; }
	double GetReferenceTime() {	return m_stSatellite.referencetime_jd; }

	void GetElements( double &dfSM, double &dfEcc, double &dfInc,
		              double &dfRAAN, double &dfPerigee, double &dfAM ); 

	void SetPosVel( double dfX,    double dfY,    double dfZ,
		            double dfXdot, double dfYdot, double dfZdot, double dfTimeJD );
	void SetPosVel( double *pdfPos, double *pdfVel, double dfTimeJD );

	void SetPosVel( double dfX,    double dfY,    double dfZ,
		            double dfXdot, double dfYdot, double dfZdot, 
					double dfTimeIntJD, double dfTimeFraJD );
	void SetPosVel( double *pdfPos, double *pdfVel, 
		            double dfTimeIntJD, double dfTimeFraJD );

	void ComputeElementFromPosVel();

	static BOOL ComputeElementsFromPosVel( double *pdfPos, double *pdfVel, double dfGM,
		                                   double *pdfElements,
										   BOOL bPartial = FALSE, double *pdfPartial = NULL );
	static BOOL ComputePosVelFromElements( double dfJDRef, double dfJD,
		                                   double *pdfElements, double dfGM, 
		                                   double *pdfPos, double *pdfVel );
	static BOOL MeanAnomaly2TrueAnomaly( double a, double e, double &M, double &f, double &E );

	static bool computePartials( double *pos, double *vel, double GM, double *element, bool E2X, double *partial );

	static void VectorProduct( double *pdfVec1, double *pdfVec2, double *pdfVec3 );

	static double ComputeMeanMotion( double a, double dfGM );

	BOOL deltaVEffect( double dxdot, double dydot, char *fileName );
};

#endif