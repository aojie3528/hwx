/***************************************************************************************

 Precession computation


***************************************************************************************/

#include <cmath>

#include "Constant.h"
#include "Precession.h"
#include "Matrix.h"

cPrecession::cPrecession()
{

}


cPrecession::~cPrecession()
{

}


/***************************************************************************************

  Fortran Subroutine: compute_precession_components_single_epoch

  compute the precession parameters from J2000 to epoch of interest JD

  see IERS Tech Note 21, p. 22

***************************************************************************************/

bool cPrecession::ComputeParameters( double dfJD, double &dfZeta,
									 double &dfZee, double &dfTheta )
{
	double dfT = ( dfJD - g_dfJ2000 ) / 36525.0,
		   dfT2 = dfT * dfT, dfT3 = dfT2 * dfT;

	dfZeta  = ( 2306.2181 * dfT + 0.30188 * dfT2 + 0.017998 * dfT3 ) * g_dfSEC2RAD;		// in rad
	dfZee =   ( 2306.2181 * dfT + 1.09468 * dfT2 + 0.018203 * dfT3 ) * g_dfSEC2RAD;		// in rad
	dfTheta = ( 2004.3109 * dfT - 0.42665 * dfT2 - 0.041833 * dfT3 ) * g_dfSEC2RAD;		// in rad

	return true;
}

/***************************************************************************************

  Compute precession matrix at given time

***************************************************************************************/

bool cPrecession::Precession( double dfJD, bool bTOR2J2000 )
{
	try
	{
		double dfT = ( dfJD - g_dfJ2000 ) / 36525.0,
			dfT2 = dfT * dfT, dfT3 = dfT2 * dfT;

		m_dfZeta  = ( 2306.2181 * dfT + 0.30188 * dfT2 + 0.017998 * dfT3 ) * g_dfSEC2RAD;		// in rad
		m_dfZee =   ( 2306.2181 * dfT + 1.09468 * dfT2 + 0.018203 * dfT3 ) * g_dfSEC2RAD;		// in rad
		m_dfTheta = ( 2004.3109 * dfT - 0.42665 * dfT2 - 0.041833 * dfT3 ) * g_dfSEC2RAD;		// in rad

		if( bTOR2J2000 ) ComputePrecessionMatrixTOR2J2000();
		else ComputePrecessionMatrixJ20002TOR();

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Compute precession matrix at given time

  dfPrecessionMatrix is vectorized by line

***************************************************************************************/

bool cPrecession::Precession( double dfJD, bool bTOR2J2000, double *dfPrecessionMatrix )
{
	try
	{
		if( !Precession( dfJD, bTOR2J2000 ) ) return false;

		if( bTOR2J2000 ) return GetPrecessionMatrixTOR2J2000( dfPrecessionMatrix );
		else return GetPrecessionMatrixJ20002TOR( dfPrecessionMatrix );
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

 Compute the precession parameters from epoch JD1 to epoch JD2

 see Section 3.21, Explanatory supplement to the astronomical almanac

 pdfPrecessionMatrix is the precession matrix from JD1 to JD2

***************************************************************************************/

bool cPrecession::Precession( double dfJD1, double dfJD2, double *pdfPrecessionMatrix )
{
	try
	{
		double dfT0 = ( dfJD1 - g_dfJ2000 ) / 36525.0, dfT02 = dfT0 * dfT0;
		double dfT = ( dfJD2 - dfJD1 ) / 36525.0, dfT2 = dfT * dfT, dfT3 = dfT2 * dfT;

		m_dfZeta = ( ( 2306.2181 + 1.39656 * dfT0 - 0.000139 * dfT02 ) * dfT +
					( 0.30188 - 0.000344 * dfT0 ) * dfT2 +
					0.017998 * dfT3 ) * g_dfSEC2RAD;		// in rad

		m_dfZee = ( ( 2306.2181 + 1.39656 * dfT0 - 0.000139 * dfT02 ) * dfT +
	  				( 1.09468 + 0.000066 * dfT0 ) * dfT2 +
					0.018203 * dfT3 ) * g_dfSEC2RAD;		// in rad

		m_dfTheta = ( ( 2004.3109 - 0.85330 * dfT0 - 0.000217 * dfT02 ) * dfT +
					( - 0.42665 - 0.000217 * dfT0 ) * dfT2 -
					0.041833 * dfT3 ) * g_dfSEC2RAD;	// in rad

		ComputePrecessionMatrixJ20002TOR();

		g_Matrix.Vec_Assign( m_pdfPrecessionMatrixJ20002TOR, pdfPrecessionMatrix, 9 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the precession parameters

***************************************************************************************/

bool cPrecession::GetParameters( double &dfZeta, double &dfZee, double &dfTheta )
{
	dfZeta = m_dfZeta;	// in rad
	dfZee = m_dfZee;
	dfTheta = m_dfTheta;

	return true;
}

/***************************************************************************************

  Get the precession matrix from TOR to J2000

  dfPrecessionMatrix is vectorized by line

***************************************************************************************/

bool cPrecession::GetPrecessionMatrixTOR2J2000( double dfZeta, double dfZee,
											    double dfTheta,
		                                        double *pdfPrecessionMatrix )
{
	try
	{
		m_dfZeta = dfZeta;	// in rad
		m_dfZee = dfZee;
		m_dfTheta = dfTheta;

		ComputePrecessionMatrixTOR2J2000();

		return GetPrecessionMatrixTOR2J2000( pdfPrecessionMatrix );
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the precession matrix from J2000 to TOR

  dfPrecessionMatrix is vectorized by line

***************************************************************************************/

bool cPrecession::GetPrecessionMatrixJ20002TOR( double dfZeta, double dfZee,
											    double dfTheta,
		                                        double *pdfPrecessionMatrix )
{
	try
	{
		m_dfZeta = dfZeta;	// in rad
		m_dfZee = dfZee;
		m_dfTheta = dfTheta;

		ComputePrecessionMatrixJ20002TOR();

		return GetPrecessionMatrixJ20002TOR( pdfPrecessionMatrix );
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the precession matrix from TOR to J2000

  dfPrecessionMatrix is vectorized by line

***************************************************************************************/

bool cPrecession::GetPrecessionMatrixTOR2J2000( double *pdfPrecessionMatrix )
{
	try
	{
		for( int i = 0; i < 9; i++ )
			pdfPrecessionMatrix[ i ] = m_pdfPrecessionMatrixTOR2J2000[ i ];

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the precession matrix from TOR to J2000

  dfPrecessionMatrix is vectorized by line

***************************************************************************************/

bool cPrecession::GetPrecessionMatrixJ20002TOR( double *pdfPrecessionMatrix )
{
	try
	{
		for( int i = 0; i < 9; i++ )
			pdfPrecessionMatrix[ i ] = m_pdfPrecessionMatrixJ20002TOR[ i ];

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Compute precession matrix from TOR to J2000

  see IERS Tech Note 21, p. 21
  Precession Rotation Matrix
  P = R3( Zeta0 ) * R2( -Theta ) * R3( Zee ) from epoch of interest to J2000

  m_pdfPrecessionMatrixTOR2J2000 is vectorized by line

***************************************************************************************/

void cPrecession::ComputePrecessionMatrixTOR2J2000()
{
	try
	{
		double dfCosZeta = cos( m_dfZeta );
		double dfSinZeta = sin( m_dfZeta );

		double dfCosZee = cos( m_dfZee );
		double dfSinZee = sin( m_dfZee );

		double dfCosTheta = cos( m_dfTheta );
		double dfSinTheta = sin( m_dfTheta );

		m_pdfPrecessionMatrixTOR2J2000[ 0 ] =   dfCosZeta * dfCosTheta * dfCosZee -
												dfSinZeta * dfSinZee;
		m_pdfPrecessionMatrixTOR2J2000[ 3 ] = - dfSinZeta * dfCosTheta * dfCosZee -
												dfCosZeta * dfSinZee;
		m_pdfPrecessionMatrixTOR2J2000[ 6 ] = - dfSinTheta * dfCosZee;

		m_pdfPrecessionMatrixTOR2J2000[ 1 ] =   dfCosZeta * dfCosTheta * dfSinZee +
												dfSinZeta * dfCosZee;
		m_pdfPrecessionMatrixTOR2J2000[ 4 ] = - dfSinZeta * dfCosTheta * dfSinZee +
												dfCosZeta * dfCosZee;
		m_pdfPrecessionMatrixTOR2J2000[ 7 ] = - dfSinTheta * dfSinZee;

		m_pdfPrecessionMatrixTOR2J2000[ 2 ] =   dfCosZeta * dfSinTheta;
		m_pdfPrecessionMatrixTOR2J2000[ 5 ] = - dfSinZeta * dfSinTheta;
		m_pdfPrecessionMatrixTOR2J2000[ 8 ] =   dfCosTheta;

		return;
	}
	catch( ... )
	{
		return;
	}
}

/***************************************************************************************

  Compute precession matrix from J2000 to TOR

  m_pdfPrecessionMatrixJ20002TOR is vectorized by line

***************************************************************************************/

void cPrecession::ComputePrecessionMatrixJ20002TOR()
{
	try
	{
		ComputePrecessionMatrixTOR2J2000();
		g_Matrix.Vec_Assign( m_pdfPrecessionMatrixTOR2J2000, m_pdfPrecessionMatrixJ20002TOR, 9 );
		g_Matrix.Mat_Transpose( m_pdfPrecessionMatrixJ20002TOR, 3, 3 );

		return;
	}
	catch( ... )
	{
	}
}