/***************************************************************************************

 about nutation and precession parameters


***************************************************************************************/

#include <cmath>

#include "NPElement.h"
#include "Constant.h"
#include "Matrix.h"

cNPElement g_NPElement;

const int nMaxTimeSpanDay = 500;
const double dfComputeStepHour = 1.0;

cNPElement::cNPElement()
{
	m_pstNPElement.clear();
	m_nNumberEpochs = 0;
}

cNPElement::~cNPElement()
{
	m_pstNPElement.clear();
}

/***************************************************************************************

  Compute the nutation and precession parameters for the given time period at the time
  of 1 hours

***************************************************************************************/

bool cNPElement::ComputeNPElement( double dfJDStart, double dfJDEnd )
{
	try
	{
		if( dfJDEnd - dfJDStart > (double) nMaxTimeSpanDay )
		{
			return false;
		}

		m_nNumberEpochs = 0;
		m_pstNPElement.clear();

		m_dfStartJD = dfJDStart;
		m_dfEndJD = dfJDEnd;

		double dfJDStep = dfComputeStepHour / 24.0;

		for( double dfJD = dfJDStart; dfJD <= dfJDEnd; dfJD += dfJDStep )
		{
			if( !ComputeAndSetNPElement( dfJD ) )
			{
				return false;
			}
		}

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Compute the nutation and precession parameters at given time

  Fortran Subroutine: compute_nutation_precession_equinox_elements

***************************************************************************************/

bool cNPElement::ComputeNPElement( double dfJD, stNPElement &tNPElement )
{
	try
	{
		if( !m_Precession.ComputeParameters( dfJD,
											 tNPElement.dfZeta,
											 tNPElement.dfZee,
											 tNPElement.dfTheta ) )
		{
			return false;
		}

		if( !m_Nutation.ComputeParameters( dfJD,
										   tNPElement.dfPSI,
										   tNPElement.dfEPS,
										   tNPElement.dfOblm ) )
		{
			return false;
		}

		tNPElement.dfJD = dfJD;
		tNPElement.dfMJD = dfJD - g_dfMJDRef;

		tNPElement.dfEqn = tNPElement.dfPSI * cos( tNPElement.dfOblm );	// in radian

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Compute the nutation and precession matrix at given time

***************************************************************************************/

bool cNPElement::ComputeNPMatrix( double dfJD, double *matrix, bool bJ2K2TOD )
{
	try
	{
		stNPElement tNPElement;
		double nut[ 9 ], pre[ 9 ];

		if( !m_Precession.ComputeParameters( dfJD,
											 tNPElement.dfZeta,
											 tNPElement.dfZee,
											 tNPElement.dfTheta ) )
		{
			return false;
		}

		if( !m_Nutation.ComputeParameters( dfJD,
										   tNPElement.dfPSI,
										   tNPElement.dfEPS,
										   tNPElement.dfOblm ) )
		{
			return false;
		}


		if( bJ2K2TOD )
		{
			if( !m_Precession.GetPrecessionMatrixJ20002TOR( tNPElement.dfZeta, tNPElement.dfZee, tNPElement.dfTheta, pre ) )
			{
				return false;
			}

			if( !m_Nutation.GetNutationMatrixMean2True( tNPElement.dfPSI, tNPElement.dfEPS, tNPElement.dfOblm, nut ) )
			{
				return false;
			}

			g_Matrix.Mat_Multiply_Mat( pre, nut, matrix, 3, 3, 3 );
		}
		else
		{
			if( !m_Precession.GetPrecessionMatrixTOR2J2000( tNPElement.dfZeta, tNPElement.dfZee, tNPElement.dfTheta, pre ) )
			{
				return false;
			}

			if( !m_Nutation.GetNutationMatrixTrue2Mean( tNPElement.dfPSI, tNPElement.dfEPS, tNPElement.dfOblm, nut ) )
			{
				return false;
			}

			g_Matrix.Mat_Multiply_Mat( nut, pre, matrix, 3, 3, 3 );
		}

		tNPElement.dfJD = dfJD;
		tNPElement.dfMJD = dfJD - g_dfMJDRef;

		tNPElement.dfEqn = tNPElement.dfPSI * cos( tNPElement.dfOblm );	// in radian

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

  Compute, and then set, the nutation and precession parameters at given time

***************************************************************************************/

bool cNPElement::ComputeAndSetNPElement( double dfJD )
{
	try
	{
		stNPElement tNPElement;

		if( !ComputeNPElement( dfJD, tNPElement ) )
		{
			return false;
		}

		m_pstNPElement.push_back( tNPElement );
		m_nNumberEpochs++;

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the nutation and precession parameters at given time

***************************************************************************************/

bool cNPElement::GetNPElement( double dfJD, stNPElement &tNPElement )
{
	try
	{
		int i, nIndex = 0, nNum = 4;
		double pdfJD[ 10 ], pdfX[ 10 ];

		for( i = 0; i < m_nNumberEpochs - 1; i++ )
		{
			if( dfJD >= m_pstNPElement[ i ].dfJD && dfJD < m_pstNPElement[ i + 1 ].dfJD )
			{
				nIndex = i;
				break;
			}
		}

		if( nIndex <= 1 ) nIndex = 1;
		else if( nIndex >= m_nNumberEpochs - 2 ) nIndex = m_nNumberEpochs - 3;

		int nIndex0 = nIndex - 1;

		// precession: Zeta
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfJD[ i ] = m_pstNPElement[ nIndex ].dfJD;
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfZeta;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate
			( nNum, pdfJD, pdfX, dfJD, tNPElement.dfZeta ) )
		{
			return false;
		}

		// precession: Zee
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfZee;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate
			( nNum, pdfJD, pdfX, dfJD, tNPElement.dfZee ) )
		{
			return false;
		}


		// precession: Theta
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfTheta;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate
			( nNum, pdfJD, pdfX, dfJD, tNPElement.dfTheta ) )
		{
			return false;
		}

		// nutation: PSI
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfPSI;
			nIndex++;
		}

		if( m_LagInterpolation.Interpolate
			( nNum, pdfJD, pdfX, dfJD, tNPElement.dfPSI ) )
		{
			return false;
		}

		// nutation: EPS
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfEPS;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate
			( nNum, pdfJD, pdfX, dfJD, tNPElement.dfEPS ) )
		{
			return false;
		}

		// nutation: Mean Obliquity
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfOblm;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate
			( nNum, pdfJD, pdfX, dfJD, tNPElement.dfOblm ) )
		{
			return false;
		}

		// Equinox equation
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfEqn;
			nIndex++;
		}

		if( m_LagInterpolation.Interpolate
			( nNum, pdfJD, pdfX, dfJD, tNPElement.dfEqn ) )
		{
			return false;
		}

		tNPElement.dfJD = dfJD;
		tNPElement.dfMJD = dfJD - g_dfMJDRef;


		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the precession parameters at the given time

  The elements are in radians

***************************************************************************************/

bool cNPElement::GetPrecessionElement( double dfJD, double &dfZeta,
									   double &dfZee, double &dfTheta )
{
	try
	{
		int i, nIndex = -10, nNum = 4;
		double pdfJD[ 10 ], pdfX[ 10 ];

		for( i = 0; i < m_nNumberEpochs - 1; i++ )
		{
			if( dfJD >= m_pstNPElement[ i ].dfJD && dfJD < m_pstNPElement[ i + 1 ].dfJD )
			{
				nIndex = i;
				break;
			}
		}

		if( nIndex == -10 )
		{
			return m_Precession.ComputeParameters( dfJD, dfZeta, dfZee, dfTheta );
		}

		if( nIndex <= 1 ) nIndex = 1;
		else if( nIndex >= m_nNumberEpochs - 2 ) nIndex = m_nNumberEpochs - 3;

		int nIndex0 = nIndex - 1;

		// precession: Zeta
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfJD[ i ] = m_pstNPElement[ nIndex ].dfJD;
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfZeta;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate( nNum, pdfJD, pdfX, dfJD, dfZeta ) )
		{
			return false;
		}

		// precession: Zee
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfZee;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate( nNum, pdfJD, pdfX, dfJD, dfZee ) )
		{
			return false;
		}

		// precession: Theta
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfTheta;
			nIndex++;
		}

		if(!m_LagInterpolation.Interpolate( nNum, pdfJD, pdfX, dfJD, dfTheta ) )
		{
			return false;
		}

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the precession matrix at the given time dfJD,

  bTOR2J2000: true if the matrix from TOR to J2000 is required, otherwise false

  pdfMatrix is the required transformation matrix, vectorized by line

***************************************************************************************/

bool cNPElement::GetPrecessionMatrix( double dfJD, double *pdfMatrix, bool bTOR2J2000 )
{
	try
	{
		double dfZeta, dfZee, dfTheta;

		if( !GetPrecessionElement( dfJD, dfZeta, dfZee, dfTheta ) )
		{
			return false;
		}

		if( bTOR2J2000 )
			m_Precession.GetPrecessionMatrixTOR2J2000( dfZeta, dfZee, dfTheta, pdfMatrix );
		else
			m_Precession.GetPrecessionMatrixJ20002TOR( dfZeta, dfZee, dfTheta, pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the nutation parameters at the given time

  The elements are in radians

***************************************************************************************/

bool cNPElement::GetNutationElement( double dfJD, double &dfPSI, double &dfEPS, double &dfOblm )
{
	try
	{
		int i, nIndex = -10, nNum = 4;
		double pdfJD[ 10 ], pdfX[ 10 ];

		for( i = 0; i < m_nNumberEpochs - 1; i++ )
		{
			if( dfJD >= m_pstNPElement[ i ].dfJD && dfJD < m_pstNPElement[ i + 1 ].dfJD )
			{
				nIndex = i;
				break;
			}
		}

		if( nIndex == -10 )
		{
			return m_Nutation.ComputeParameters( dfJD, dfPSI, dfEPS, dfOblm );
		}

		if( nIndex <= 1 ) nIndex = 1;
		else if( nIndex >= m_nNumberEpochs - 2 ) nIndex = m_nNumberEpochs - 3;

		int nIndex0 = nIndex - 1;

		// nutation: PSI
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfJD[ i ] = m_pstNPElement[ nIndex ].dfJD;
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfPSI;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate( nNum, pdfJD, pdfX, dfJD, dfPSI ) )
		{
			return false;
		}

		// nutation: EPS
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfEPS;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate( nNum, pdfJD, pdfX, dfJD, dfEPS ) )
		{
			return false;
		}

		// nutation: Mean Obliquity
		nIndex = nIndex0;
		for( i = 0; i < nNum; i++ )
		{
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfOblm;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate( nNum, pdfJD, pdfX, dfJD, dfOblm ) )
		{
			return false;
		}

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the nutation matrix at the given time

  pdfMatrix is vectorized by line

***************************************************************************************/

bool cNPElement::GetNutationMatrix( double dfJD, double *pdfMatrix, bool bTrue2Mean )
{
	try
	{
		double dfPSI, dfEPS, dfOblm;

		if( !GetNutationElement( dfJD, dfPSI, dfEPS, dfOblm ) )
		{
			return false;
		}

		if( bTrue2Mean )
			m_Nutation.GetNutationMatrixTrue2Mean( dfPSI, dfEPS, dfOblm, pdfMatrix );
		else
			m_Nutation.GetNutationMatrixMean2True( dfPSI, dfEPS, dfOblm, pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the equinox equation at the given time

  dfEqn in radian

***************************************************************************************/

bool cNPElement::GetEquinoxEquation( double dfJD, double &dfEqn )
{
	try
	{
		if( m_nNumberEpochs < 1 )
		{
			stNPElement tNPElement;

			if( !ComputeNPElement( dfJD, tNPElement ) )
			{
				return false;
			}

			dfEqn = tNPElement.dfEqn;

			return true;
		}


		int i, nIndex = 0, nNum = 4;
		double pdfJD[ 10 ], pdfX[ 10 ];

		for( i = 0; i < m_nNumberEpochs - 1; i++ )
		{
			if( dfJD >= m_pstNPElement[ i ].dfJD && dfJD < m_pstNPElement[ i + 1 ].dfJD )
			{
				nIndex = i;
				break;
			}
		}

		if( nIndex <= 1 ) nIndex = 1;
		else if( nIndex >= m_nNumberEpochs - 2 ) nIndex = m_nNumberEpochs - 3;

		nIndex = nIndex - 1;
		for( i = 0; i < nNum; i++ )
		{
			pdfJD[ i ] = m_pstNPElement[ nIndex ].dfJD;
			pdfX[ i ] = m_pstNPElement[ nIndex ].dfEqn;
			nIndex++;
		}

		if( !m_LagInterpolation.Interpolate( nNum, pdfJD, pdfX, dfJD, dfEqn ) )
		{
			return false;
		}

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

bool cNPElement::GetMatrixTOD2J2000( double dfJD, double *pdfMatrix )
{
	double nutation[ 9 ], precession[ 9 ];

	if( GetNutationMatrix( dfJD, nutation, true ) &&
		GetPrecessionMatrix( dfJD, precession, true ) )
	{
		g_Matrix.Mat_Multiply_Mat( precession, nutation, pdfMatrix, 3, 3, 3 );
		return true;
	}

	return false;
}

bool cNPElement::GetMatrixJ20002TOD( double dfJD, double *pdfMatrix )
{
	if( !GetMatrixTOD2J2000( dfJD, pdfMatrix ) ) return false;

	g_Matrix.Mat_Transpose( pdfMatrix, 3, 3 );

	return true;
}