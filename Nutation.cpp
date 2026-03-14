/***************************************************************************************

 compute nutation parameters


***************************************************************************************/

#include <cmath>

#include "Constant.h"
#include "Nutation.h"
#include "AstroArgument.h"
#include "Matrix.h"
#include "EOP.h"

cNutation::cNutation( int nModel )
{
	m_nModel = nModel;
}

cNutation::~cNutation()
{
}

/***************************************************************************************

  Fortran Subroutine: compute_nutation_components_single_epoch

  Compute the nutation parameters at the given time

***************************************************************************************/

bool cNutation::ComputeParameters( double dfJD, double &dfPSI, double &dfEPS,
								   double &dfOblm, int nModel )
{
	switch( m_nModel )
	{
		case 1980:
			if( !IAU1980( dfJD ) )
			{
				return false;
			}
			break;

		case 1996:
			if( !IERS1996( dfJD ) )
			{
				return false;
			}
			break;
	}

	ComputeMeanObliquity( dfJD );

	dfPSI = m_dfPSI;
	dfEPS = m_dfEPS;
	dfOblm = m_dfOblm;

	return true;
}

/***************************************************************************************

  Compute nutation matrix at given time

***************************************************************************************/

bool cNutation::Nutation( double dfJD, bool bTrue2Mean, int nModel )
{
	m_nModel = nModel;

	switch( m_nModel )
	{
		case 1980:
			if( !IAU1980( dfJD ) )
			{
				return false;
			}
			break;

		case 1996:
			if( !IERS1996( dfJD ) )
			{
				return false;
			}
			break;
	}

	ComputeMeanObliquity( dfJD );
	if( bTrue2Mean ) ComputeNutationMatrixTrue2Mean();
	else ComputeNutationMatrixMean2True();

	return true;
}

/***************************************************************************************

  Compute nutation matrix at given time

  pdfNutationMatrix is vectorized by line

***************************************************************************************/

bool cNutation::Nutation( double dfJD, bool bTrue2Mean,
						  double *pdfNutationMatrix, int nModel )
{
	if( !Nutation( dfJD, bTrue2Mean, nModel ) )
	{
		return false;
	}

	if( bTrue2Mean ) return GetNutationMatrixTrue2Mean( pdfNutationMatrix );
	else return GetNutationMatrixMean2True( pdfNutationMatrix );
}

/***************************************************************************************

  Compute mean obliquity at given time

***************************************************************************************/

void cNutation::ComputeMeanObliquity( double dfJD )
{
	double dfT, dfT2, dfT3;

	dfT = ( dfJD - g_dfJ2000 ) / 36525.0;
	dfT2 = dfT * dfT;
	dfT3 = dfT2 * dfT;

	// lieske's expression for the obliquity
	// P.22 IERS TECHNICAL NOTE 21, IERS CONVENTIONS 1996, D.D McCARTHY
	m_dfOblm = ( 84381.448 - 46.8150 * dfT - 0.00059 * dfT2 +
		         0.001813 * dfT3 ) * g_dfSEC2RAD; // in radian

}

/***************************************************************************************

  Compute mean obliquity at given time

***************************************************************************************/

void cNutation::ComputeMeanObliquity( double dfJD, double &dfOblm )
{
	double dfT, dfT2, dfT3;

	dfT = ( dfJD - g_dfJ2000 ) / 36525.0;
	dfT2 = dfT * dfT;
	dfT3 = dfT2 * dfT;

	// lieske's expression for the obliquity
	// P.22 IERS TECHNICAL NOTE 21, IERS CONVENTIONS 1996, D.D McCARTHY
	dfOblm = m_dfOblm =
		( 84381.448 - 46.8150 * dfT - 0.00059 * dfT2 +
		  0.001813 * dfT3 ) * g_dfSEC2RAD; // in radian
}

/***************************************************************************************

  Get the nutation parameters

***************************************************************************************/

bool cNutation::GetParameters( double &dfPSI, double &dfEPS )
{
	dfPSI = m_dfPSI;
	dfEPS = m_dfEPS;

	return true;
}

/***************************************************************************************

  Compute the nutation matrix from True to mean by given the nutation parameters

  pdfNutationMatrix is vectorized by line

***************************************************************************************/

bool cNutation::GetNutationMatrixTrue2Mean( double dfPSI, double dfEPS, double dfOblm,
	                                        double *pdfNutationMatrix )
{
	m_dfPSI = dfPSI;
	m_dfEPS = dfEPS;
	m_dfOblm = dfOblm;

	ComputeNutationMatrixTrue2Mean();

	return GetNutationMatrixTrue2Mean( pdfNutationMatrix );
}

/***************************************************************************************

  Compute the nutation matrix from True to mean by given the nutation parameters

  pdfNutationMatrix is vectorized by line

***************************************************************************************/

bool cNutation::GetNutationMatrixMean2True( double dfPSI, double dfEPS, double dfOblm,
		                                    double *pdfNutationMatrix )
{
	m_dfPSI = dfPSI;
	m_dfEPS = dfEPS;
	m_dfOblm = dfOblm;

	ComputeNutationMatrixMean2True();

	return GetNutationMatrixMean2True( pdfNutationMatrix );

}

/***************************************************************************************

  Get nutation matrix from true to mean

  pdfNutationMatrix is vectorized by line

***************************************************************************************/

bool cNutation::GetNutationMatrixTrue2Mean( double *pdfNutationMatrix )
{
	for( int i = 0; i < 9; i++ )
		pdfNutationMatrix[ i ] = m_pdfNutationMatrixTrue2Mean[ i ];

	return true;
}

/***************************************************************************************

  Get nutation matrix from mean to true

   pdfNutationMatrix is vectorized by line

***************************************************************************************/

bool cNutation::GetNutationMatrixMean2True( double *pdfNutationMatrix )
{
	for( int i = 0; i < 9; i++ )
		pdfNutationMatrix[ i ] = m_pdfNutationMatrixMean2True[ i ];

	return true;
}

/***************************************************************************************

  compute nutation matrix from true to MEAN at epoch

  see p.21 IERS TECH NOTE 21

  N = R1(-eA ) R3(dPSI) R1( eA + dE )

  m_pdfNutationMatrixTrue2Mean is vectorized by line

***************************************************************************************/

void cNutation::ComputeNutationMatrixTrue2Mean()
{
	double dfCosOBM,
		   dfSinOBM,
		   dfCosOBT,
		   dfSinOBT,
		   dfCosPSI,
		   dfSinPSI;

	double dfObliquity = m_dfOblm + m_dfEPS;

	dfCosOBM = cos( m_dfOblm );
	dfSinOBM = sin( m_dfOblm );
	dfCosOBT = cos( dfObliquity );
	dfSinOBT = sin( dfObliquity );
	dfCosPSI = cos( m_dfPSI );
	dfSinPSI = sin( m_dfPSI );

	m_pdfNutationMatrixTrue2Mean[ 0 ] =   dfCosPSI;
	m_pdfNutationMatrixTrue2Mean[ 3 ] = - dfSinPSI * dfCosOBM;
	m_pdfNutationMatrixTrue2Mean[ 6 ] = - dfSinPSI * dfSinOBM;
	m_pdfNutationMatrixTrue2Mean[ 1 ] =   dfSinPSI * dfCosOBT;
	m_pdfNutationMatrixTrue2Mean[ 4 ] =   dfCosPSI * dfCosOBM * dfCosOBT + dfSinOBM * dfSinOBT;
	m_pdfNutationMatrixTrue2Mean[ 7 ] =   dfCosPSI * dfSinOBM * dfCosOBT - dfCosOBM * dfSinOBT;
	m_pdfNutationMatrixTrue2Mean[ 2 ] =   dfSinPSI * dfSinOBT;
	m_pdfNutationMatrixTrue2Mean[ 5 ] =   dfCosPSI * dfCosOBM * dfSinOBT - dfSinOBM * dfCosOBT;
	m_pdfNutationMatrixTrue2Mean[ 8 ] =   dfCosPSI * dfSinOBM * dfSinOBT + dfCosOBM * dfCosOBT;
}

/***************************************************************************************

  compute nutation matrix from mean to true at epoch

  m_pdfNutationMatrixMean2True is vectorized by line

***************************************************************************************/

void cNutation::ComputeNutationMatrixMean2True( )
{
	ComputeNutationMatrixTrue2Mean();
	g_Matrix.Vec_Assign
		( m_pdfNutationMatrixTrue2Mean, m_pdfNutationMatrixMean2True, 9 );
	g_Matrix.Mat_Transpose( m_pdfNutationMatrixMean2True, 3, 3 );
}


/***************************************************************************************

  IAU 1980 nutation theory

***************************************************************************************/

bool cNutation::IAU1980( double dfJD )
{
	try
	{
		int i, j;
		double dfL, dfLp, dfF, dfD, dfOm, dfArg;
		double dfT = ( dfJD - g_dfJ2000 ) / 36525.0;

		m_dfPSI = 0.0;
		m_dfEPS = 0.0;

		m_AstroArg.ComputeFundamentals( dfT, dfL, dfLp, dfF, dfD, dfOm );

		for( i = 0; i < 106; i++ )
		{
			j = i * 9;

			dfArg =  dfX1980[ j ] * dfL    +
					 dfX1980[ j + 1 ] * dfLp +
					 dfX1980[ j + 2 ] * dfF  +
					 dfX1980[ j + 3 ] * dfD  +
					 dfX1980[ j + 4 ] * dfOm;

			dfArg = fmod( dfArg, g_dfTWOPI );

			m_dfPSI += ( dfX1980[ j + 5 ] + dfX1980[ j + 6 ] * dfT ) * sin( dfArg ); // IN 0".0001
			m_dfEPS += ( dfX1980[ j + 7 ] + dfX1980[ j + 8 ] * dfT ) * cos( dfArg ); // IN 0".0001
		}

		m_dfPSI /= 10000.0;   // IN arc seconds
		m_dfEPS /= 10000.0;   // IN arc seconds

		m_dfPSI *= g_dfSEC2RAD;	// in radian
		m_dfEPS *= g_dfSEC2RAD;

		// correct for the dpsi and deps given in the IERS bulletin

		double dfDPSI, dfDEPS;

		if( g_EOP.GetDPSIDEPS( dfJD, dfDPSI, dfDEPS, true ) )
		{
			m_dfPSI += dfDPSI;
			m_dfEPS += dfDEPS;
		}

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  IERS 1996 nutation theory

***************************************************************************************/

bool cNutation::IERS1996( double dfJD )
{
	try
	{
		int i, j;
		double dfL, dfLp, dfF, dfD, dfOm, dfArg, dfCosArg, dfSinArg;
		double dfT = ( dfJD - g_dfJ2000 ) / 36525.0;

		m_dfPSI = 0.0;
		m_dfEPS = 0.0;

		m_AstroArg.ComputeFundamentals( dfT, dfL, dfLp, dfF, dfD, dfOm );

		for( i = 0; i < 263; i++ )
		{
			j = i * 11;

			dfArg =  dfX1996[ j ] * dfL    +
					 dfX1996[ j + 1 ] * dfLp +
					 dfX1996[ j + 2 ] * dfF  +
					 dfX1996[ j + 3 ] * dfD  +
					 dfX1996[ j + 4 ] * dfOm;

			dfArg = fmod( dfArg, g_dfTWOPI );

			dfCosArg = cos( dfArg );
			dfSinArg = sin( dfArg );

			m_dfPSI +=
				( dfX1996[ j + 5 ] + dfX1996[ j + 6 ] * dfT ) * dfSinArg +
				  dfX1996[ j + 9 ] * dfCosArg; // IN 0.001 mas
			m_dfEPS +=
				( dfX1996[ j + 7 ] + dfX1996[ j + 8 ] * dfT ) * dfCosArg +
				  dfX1996[ j + 10 ] * dfSinArg; // IN 0".0001
		}

		m_dfPSI *= 1.0e6;   // IN arc seconds
		m_dfEPS *= 1.0e6;   // IN arc seconds

		m_dfPSI *= g_dfSEC2RAD;	// in radian
		m_dfEPS *= g_dfSEC2RAD;

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


