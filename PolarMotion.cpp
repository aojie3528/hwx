/***************************************************************************************

 polar motion parameters and rotation matrix

  Jizhang Sang, EOS

 (C) COPYRIGHT -Electro Optic Systems Pty. Ltd. Australia (EOS) 2002,
 All rights reserved. No part of this program may be photocopied, 
 reproduced or translated to another program language without the prior
 written consent of EOS

 Reference: IERS TECH NOTE 21

***************************************************************************************/


#include "PolarMotion.h"


class cPolarMotion g_PolarMotion;

cPolarMotion::cPolarMotion()
{
	m_nNumber = 0;
	m_pstPM.clear();
}

cPolarMotion::~cPolarMotion()
{
	m_pstPM.clear();
}

/***************************************************************************************

  Set the polar motion parameter at given time

***************************************************************************************/

bool cPolarMotion::SetPM()
{
	try
	{	
		double dfJD, dfX, dfY;

		m_nNumber = 0;
		m_pstPM.clear();
		stPolarMotion tPM;

		for( int i = 0; i < g_EOP.GetNumberEOP(); i++ )
		{
			if( !g_EOP.GetPolarMotion( i, dfJD, dfX, dfY, false ) ) 
			{
				return false;
			}

			tPM.dfJD = dfJD;
			tPM.dfX = dfX;
			tPM.dfY = dfY;

			m_pstPM.push_back( tPM );

			m_nNumber++;
		}

		return true;
	}
	catch( ... )
	{
		return false;
	}
}



/*
20201014
*/
bool cPolarMotion::SetPMLOCAL(cEOP & myeop)
{
	try
	{
		double dfJD, dfX, dfY;

		m_nNumber = 0;
		m_pstPM.clear();
		stPolarMotion tPM;

		for (int i = 0; i < g_EOP.GetNumberEOP(); i++)
		{
			if (!myeop.GetPolarMotion(i, dfJD, dfX, dfY, false))
			{
				return false;
			}

			tPM.dfJD = dfJD;
			tPM.dfX = dfX;
			tPM.dfY = dfY;

			m_pstPM.push_back(tPM);

			m_nNumber++;
		}

		return true;
	}
	catch (...)
	{
		return false;
	}
}






bool cPolarMotion::SetPM( double dfJD, double dfX, double dfY )
{
	try
	{
		stPolarMotion tPM;

		tPM.dfJD = dfJD;
		tPM.dfX = dfX;
		tPM.dfY = dfY;

		m_pstPM.push_back( tPM );

		m_nNumber++;

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the polar motion parameters at the given time

***************************************************************************************/

bool cPolarMotion::GetPM( double dfJD, double &dfX, double &dfY, bool bMean2True )
{
	try
	{
		int i, nIndex = 0, nNum = 4;
		double pdfJD[ 10 ], pdfX[ 10 ], pdfY[ 10 ];

		for( i = 0; i < m_nNumber - 1; i++ )
		{
			if( dfJD >= m_pstPM[ i ].dfJD && dfJD < m_pstPM[ i + 1 ].dfJD ) 
			{
				nIndex = i;
				break;
			}
		}

		if( nIndex <= 1 ) nIndex = 1;
		else if( nIndex >= m_nNumber - 2 ) nIndex = m_nNumber - 3;

		nIndex = nIndex - 1;

		for( i = 0; i < nNum; i++ )
		{
			pdfJD[ i ] = m_pstPM[ nIndex ].dfJD;
			pdfX[ i ] = m_pstPM[ nIndex ].dfX;
			pdfY[ i ] = m_pstPM[ nIndex ].dfY;
			nIndex++;
		}

		m_LagInterpolation.Interpolate( nNum, pdfJD, pdfX, dfJD, dfX );
		m_LagInterpolation.Interpolate( nNum, pdfJD, pdfY, dfJD, dfY );

		m_dfJD = dfJD;
		m_dfX = dfX;
		m_dfY = dfY;

		// correct for the diurnal and semi-diurnal tidal variations
		double eop[ 3 ];
		g_EOP.ComputeTidalVariation( dfJD, eop[ 0 ], eop[ 1 ], eop[ 2 ] );

		m_dfX += eop[ 0 ] / 206264.8;
		m_dfY += eop[ 1 ] / 206264.8;

		if( bMean2True ) ComputePMMatrixMean2True();
		else ComputePMMatrixTrue2Mean();

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the polar motion matrix from EXEF mean to true at given time

***************************************************************************************/

bool cPolarMotion::GetPMMatrixMean2True( double *pdfPMMatrix )
{
	try
	{
		for( int i = 0; i < 9; i++ )
			pdfPMMatrix[ i ] = m_pdfPMMatrixMean2True[ i ];

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the polar motion matrix from EXEF mean to true

***************************************************************************************/

bool cPolarMotion::GetPMMatrixMean2True( double dfJD, double *pdfPMMatrix )
{
	try
	{
		double dfX, dfY;

		GetPM( dfJD, dfX, dfY, true );

		GetPMMatrixMean2True( pdfPMMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the polar motion matrix from EXEF true to mean at teh given time

***************************************************************************************/

bool cPolarMotion::GetPMMatrixTrue2Mean( double dfJD, double *pdfPMMatrix )
{
	try
	{
		double dfX, dfY;

		GetPM( dfJD, dfX, dfY, false );

		GetPMMatrixTrue2Mean( pdfPMMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get the polar motion matrix from EXEF true to mean

***************************************************************************************/

bool cPolarMotion::GetPMMatrixTrue2Mean( double *pdfPMMatrix )
{
	try
	{
		for( int i = 0; i < 9; i++ )
			pdfPMMatrix[ i ] = m_pdfPMMatrixTrue2Mean[ i ];

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

// Rotation matrix from mean pole to true pole
// R = R1(Y) R2(X), see p. 21, IERS TECH NOTE 21
//       1    0   -x
//   = x*y    1    y
// 	     x   -y    1

/***************************************************************************************

  Compute polar motion matrix from ECEF mean to ECEF true

  Fortran subroutine: compute_polar_motion_matrix

***************************************************************************************/

bool cPolarMotion::ComputePMMatrixMean2True()
{
	m_pdfPMMatrixMean2True[ 0 ] = 1.0;
	m_pdfPMMatrixMean2True[ 1 ] = 0.0;
	m_pdfPMMatrixMean2True[ 2 ] = -m_dfX;
	m_pdfPMMatrixMean2True[ 3 ] = m_dfX * m_dfY;
	m_pdfPMMatrixMean2True[ 4 ] = 1.0;
	m_pdfPMMatrixMean2True[ 5 ] = m_dfY;
	m_pdfPMMatrixMean2True[ 6 ] = m_dfX;
	m_pdfPMMatrixMean2True[ 7 ] = -m_dfY;
	m_pdfPMMatrixMean2True[ 8 ] = 1.0;

/*
    double sinx = sin( m_dfX );
    double cosx = cos( m_dfX );
    double siny = sin( m_dfY );
    double cosy = cos( m_dfY );

    m_pdfPMMatrixMean2True[ 0 ] = cosx;
    m_pdfPMMatrixMean2True[ 1 ] = sinx * siny;
    m_pdfPMMatrixMean2True[ 2 ] = sinx * cosy;
    m_pdfPMMatrixMean2True[ 3 ] = 0.0;
    m_pdfPMMatrixMean2True[ 4 ] = cosy;
    m_pdfPMMatrixMean2True[ 5 ] = -siny;
    m_pdfPMMatrixMean2True[ 6 ] = -sinx;
    m_pdfPMMatrixMean2True[ 7 ] = siny * cosx;
    m_pdfPMMatrixMean2True[ 8 ] = cosy * cosx;
*/

	return true;
}

/***************************************************************************************

  Compute polar motion matrix from ECEF mean to ECEF true

***************************************************************************************/

bool cPolarMotion::ComputePMMatrixMean2True( double dfX, double dfY, double *pdfMatrix )
{
	try
	{
		m_dfX = dfX;
		m_dfY = dfY;

		ComputePMMatrixMean2True();

		GetPMMatrixMean2True( pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Compute polar motion matrix from ECEF true to ECEF mean

***************************************************************************************/

bool cPolarMotion::ComputePMMatrixTrue2Mean()
{
	m_pdfPMMatrixTrue2Mean[ 0 ] = 1.0;
	m_pdfPMMatrixTrue2Mean[ 3 ] = 0.0;
	m_pdfPMMatrixTrue2Mean[ 6 ] = -m_dfX;
	m_pdfPMMatrixTrue2Mean[ 1 ] = m_dfX * m_dfY;
	m_pdfPMMatrixTrue2Mean[ 4 ] = 1.0;
	m_pdfPMMatrixTrue2Mean[ 7 ] = m_dfY;
	m_pdfPMMatrixTrue2Mean[ 2 ] = m_dfX;
	m_pdfPMMatrixTrue2Mean[ 5 ] = -m_dfY;
	m_pdfPMMatrixTrue2Mean[ 8 ] = 1.0;

/*
    double sinx = sin( m_dfX );
    double cosx = cos( m_dfX );
    double siny = sin( m_dfY );
    double cosy = cos( m_dfY );

    m_pdfPMMatrixTrue2Mean[ 0 ] = cosx;
    m_pdfPMMatrixTrue2Mean[ 3 ] = sinx * siny;
    m_pdfPMMatrixTrue2Mean[ 6 ] = sinx * cosy;
    m_pdfPMMatrixTrue2Mean[ 1 ] = 0.0;
    m_pdfPMMatrixTrue2Mean[ 4 ] = cosy;
    m_pdfPMMatrixTrue2Mean[ 7 ] = -siny;
    m_pdfPMMatrixTrue2Mean[ 2 ] = -sinx;
    m_pdfPMMatrixTrue2Mean[ 5 ] = siny * cosx;
    m_pdfPMMatrixTrue2Mean[ 8 ] = cosy * cosx;
*/

	return true;
}

/***************************************************************************************

  Compute polar motion matrix from ECEF true to ECEF mean

***************************************************************************************/

bool cPolarMotion::ComputePMMatrixTrue2Mean( double dfX, double dfY, double *pdfMatrix )
{
	try
	{
		m_dfX = dfX;
		m_dfY = dfY;

		ComputePMMatrixTrue2Mean();

		GetPMMatrixTrue2Mean( pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}
