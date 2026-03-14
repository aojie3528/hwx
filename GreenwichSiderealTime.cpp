/***************************************************************************************

 compute Greenwich sidereal time


 Reference: IERS TECH NOTE 21

***************************************************************************************/

#include <math.h>
#include "constant.h"
#include "GreenwichSiderealTime.h"

cGreenwichST::cGreenwichST()
{

}

cGreenwichST::~cGreenwichST()
{

}

/***************************************************************************************

  Compute GST at the given time, the equinox equation is ignored.
 
  Time should be in UT1 system.

***************************************************************************************/

double cGreenwichST::ComputeGST( double dfTJD )
{
	double dfDt;		
	
	double dfJD = ( ( int ) dfTJD ) * 1.0;
	
	double dfFday = dfTJD - dfJD;

	if( dfFday >= 0.5 )
	{
		dfJD += 0.5;
		dfFday -= 0.5;
	}
	else
	{
		dfJD -= 0.5;
		dfFday += 0.5;
	}
	
	// against 2000 January 1st 12 UT1 (JD=2451545.0), p.21 IERS tech note 21
	dfDt = dfJD - g_dfJ2000;	
	double dfRA = 100.4606184 + 0.9856473663 * dfDt + 2.908e-13 * dfDt * dfDt;

	dfRA *= g_dfDEG2RAD;		// in radians
	
	double dfR = 1.002737909350795 + 5.9006e-11 / 36525.0 * dfDt;
	double dfDSRT = dfR * dfFday * g_dfTWOPI;		// in radians

	dfRA += dfDSRT;		// consider earth rotation from midnight
	
	dfRA = fmod( dfRA, g_dfTWOPI );
	
	if( dfRA < 0.0 ) dfRA += g_dfTWOPI;

	return dfRA;
}

/***************************************************************************************

  Compute GST at the given time
 
  Time is in UT1 system

***************************************************************************************/

BOOL cGreenwichST::ComputeGST( double dfIntJD, double dfFractionJD, double dfEquinoxEquation,
                               double &dfGST )
{
	double dfJD = dfIntJD, dfDay = dfFractionJD;

	if( dfDay >= 0.5 )
	{
		dfJD += 0.5;
		dfDay -= 0.5;
	}
	else
	{
		dfJD -= 0.5;
		dfDay += 0.5;
	}
	
	double dfDT = dfJD - g_dfJ2000;

	dfGST = 100.460618375 + 0.98564736628633356 * dfDT + 
            2.90788e-13 * dfDT * dfDT - 5.3e-22 * dfDT * dfDT * dfDT;	// in deg
	dfGST *= g_dfDEG2RAD;	// in radians
 
	dfDT /= 36525.0;

	double dfR = 1.002737909350795 + 5.9006e-11 * dfDT - 5.9e-15 * dfDT * dfDT;

	dfGST += ( dfR * dfDay * g_dfTWOPI );

	dfGST += dfEquinoxEquation;

    dfDT = ( dfJD + dfDay - g_dfJ2000 ) / 36525.0;

	double dfOmega = 2.1824391966 - 33.7570446126362 * dfDT + 3.62262478e-05 * dfDT * dfDT;
    dfOmega = fmod( dfOmega, g_dfTWOPI );
    if( dfOmega < 0.0 ) dfOmega += g_dfTWOPI;
	 
	dfGST += ( 1.279908e-8 * sin( dfOmega ) + 3.054326e-10 * sin( 2.0 * dfOmega ) );
	
    dfGST = fmod( dfGST, g_dfTWOPI );

	return TRUE;
}

/***************************************************************************************

  Compute GST at the given time
 
  Time is in UT1 system

***************************************************************************************/

BOOL cGreenwichST::ComputeGST( double dfIntJD, double dfFractionJD, double dfEquinoxEquation,
                               double &dfGST,  double &dfSinGST,    double &dfCosGST )
{
	if( !ComputeGST( dfIntJD, dfFractionJD, dfEquinoxEquation, dfGST ) ) return FALSE;

	dfSinGST = sin( dfGST );
	dfCosGST = cos( dfGST );

	return TRUE;
}

/***************************************************************************************

  Compute GST matrix from ECEF to TOD at given time 

  Time is in UT1 system

***************************************************************************************/

BOOL cGreenwichST::ComputeGSTMatrixECEF2TOD( double dfIntJD, double dfFractionJD, double dfEquinoxEquation,
                                             double *pdfGSTMatrix )
{
	double dfGST, dfSinGST, dfCosGST;

	if( !ComputeGST( dfIntJD, dfFractionJD, dfEquinoxEquation, dfGST, dfSinGST, dfCosGST ) ) return FALSE;

	pdfGSTMatrix[ 0 ] = pdfGSTMatrix[ 4 ] = dfCosGST;
	pdfGSTMatrix[ 2 ] = pdfGSTMatrix[ 5 ] = pdfGSTMatrix[ 6 ] = pdfGSTMatrix[ 7 ] = 0.0;
	pdfGSTMatrix[ 1 ] = -dfSinGST;
	pdfGSTMatrix[ 3 ] = dfSinGST;
	pdfGSTMatrix[ 8 ] = 1.0;


	return TRUE;
}

/***************************************************************************************

  Compute GST matrix from TOD to ECEF at given time 

***************************************************************************************/

BOOL cGreenwichST::ComputeGSTMatrixTOD2ECEF( double dfIntJD, double dfFractionJD, double dfEquinoxEquation,
                                             double *pdfGSTMatrix )
{
	if( !ComputeGSTMatrixECEF2TOD( dfIntJD, dfFractionJD, dfEquinoxEquation, pdfGSTMatrix ) )
		return FALSE;

	pdfGSTMatrix[ 1 ] *= (-1.0);
	pdfGSTMatrix[ 3 ] *= (-1.0);

	return TRUE;
}

