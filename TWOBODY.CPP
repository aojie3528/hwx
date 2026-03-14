/***************************************************************************

 Orbit computation using two body theory


 Reference: IAU SOFA EPV00 for Model 1
  
***************************************************************************/

#include <stdio.h>
#include <math.h>

#include "twobody.h"
#include "coortran.h"
#include "GreenwichSiderealTime.h"
#include "constant.h"

cTwoBody::cTwoBody()
{
	m_dfTransCoe = ( g_dfEarthSemiMajor * dfEarthGMIOD ) / 60.0;
}


cTwoBody::~cTwoBody()
{
}


void cTwoBody::SetElement( TwoBodyOrbitElement &stSatellite )
{
    m_stSatellite.semi_major = stSatellite.semi_major;
    m_stSatellite.eccentricity = stSatellite.eccentricity;
    m_stSatellite.inclination = stSatellite.inclination;
    m_stSatellite.perigee = stSatellite.perigee;
    m_stSatellite.Longi_RisingNode = stSatellite.Longi_RisingNode;
    m_stSatellite.MeanAnomaly = stSatellite.MeanAnomaly;
    m_stSatellite.referencetime_jd = stSatellite.referencetime_jd;

	m_stSatellite.referencetime_jd_int = stSatellite.referencetime_jd_int;
	m_stSatellite.referencetime_jd_fra = stSatellite.referencetime_jd_fra;

	double t = m_stSatellite.referencetime_jd_int + m_stSatellite.referencetime_jd_fra;

	if( fabs( m_stSatellite.referencetime_jd - g_dfJ2000 ) > 500.0 * 365.0 )
	{
		m_stSatellite.referencetime_jd = t;
	}
	else
	{
		if( fabs( m_stSatellite.referencetime_jd - t ) > 1.0e-5 )
		{
			m_stSatellite.referencetime_jd_int = 
				(double)( (int) m_stSatellite.referencetime_jd );
			m_stSatellite.referencetime_jd_fra = 
				m_stSatellite.referencetime_jd - 
				m_stSatellite.referencetime_jd_int;
		}
	}

	m_stSatellite.GM = stSatellite.GM;		
}


void cTwoBody::ComputeInertialXYZ( double dfTimeJD )
{
	try
	{
		double dfMM, dfDeltaT;

		if( m_stSatellite.GM == 1 ) 
			dfMM = sqrt( g_dfEarthGM / pow( m_stSatellite.semi_major, 3.0 ) );
		else dfMM = sqrt( dfMuIOD / m_stSatellite.semi_major ) / 
			        m_stSatellite.semi_major;	 		

		dfDeltaT = dfTimeJD - m_stSatellite.referencetime_jd;

		if( m_stSatellite.GM == 1 ) dfDeltaT *= 86400.0;
		else dfDeltaT *= dfEarthGMIOD * 1440.0;

		double mean_anomaly0 = m_stSatellite.MeanAnomaly + dfMM * dfDeltaT;
		double eo = mean_anomaly0;

		double ecc_ano;
		for( int j = 0; j < 50; j++ )
		{ 
			ecc_ano = mean_anomaly0 + m_stSatellite.eccentricity * sin( eo );  
			eo = ecc_ano; 
		}

		double cE = cos( ecc_ano );
		double sE = sin( ecc_ano );
		double qee = sqrt( 1.0 - m_stSatellite.eccentricity * m_stSatellite.eccentricity );
		double xw = m_stSatellite.semi_major * ( cE - m_stSatellite.eccentricity );	// (3.69)
		double yw = m_stSatellite.semi_major * qee * sE;							// (3.69)
		double Edot = dfMM / ( 1.0 - m_stSatellite.eccentricity * cE );		
		double xwdot = - m_stSatellite.semi_major * Edot * sE;						// (3.70)
		double ywdot =   m_stSatellite.semi_major * Edot * qee * cE;				// (3.70)

		double sO = sin( m_stSatellite.Longi_RisingNode );
		double cO = cos( m_stSatellite.Longi_RisingNode );
		double sw = sin( m_stSatellite.perigee );
		double cw = cos( m_stSatellite.perigee );
		double si = sin( m_stSatellite.inclination );
		double ci = cos( m_stSatellite.inclination );

		double Px = cO * cw - sO * sw * ci;
		double Py = sO * cw + cO * sw * ci;										// (3.43)
		double Pz = sw * si;

		double Qx = - cO * sw - sO * cw * ci;
		double Qy = - sO * sw + cO * cw * ci;									// (3.44)
		double Qz =   cw * si;

		m_dfX = xw * Px + yw * Qx;
		m_dfY = xw * Py + yw * Qy;												// (3.57)
		m_dfZ = xw * Pz + yw * Qz;

		m_dfXdot = xwdot * Px + ywdot * Qx;
		m_dfYdot = xwdot * Py + ywdot * Qy;
		m_dfZdot = xwdot * Pz + ywdot * Qz;

		return;
	}
	catch( ... )
	{
		return;
	}
}

void cTwoBody::ComputeInertialXYZ( double dfIntJD, double dfFraJD )
{
	try
	{
		double dfMM, dfDeltaT;

		if( m_stSatellite.GM == 1 ) 
			dfMM = sqrt( g_dfEarthGM / pow( m_stSatellite.semi_major, 3.0 ) );
		else dfMM = sqrt( dfMuIOD / m_stSatellite.semi_major ) / 
			        m_stSatellite.semi_major;	 		

		dfDeltaT = dfIntJD - m_stSatellite.referencetime_jd_int + 
			       dfFraJD - m_stSatellite.referencetime_jd_fra;
		
		if( m_stSatellite.GM == 1 ) dfDeltaT *= 86400.0;
		else dfDeltaT *= dfEarthGMIOD * 1440.0;

		double mean_anomaly0 = m_stSatellite.MeanAnomaly + dfMM * dfDeltaT;
		double eo = mean_anomaly0;

		double ecc_ano;
		for( int j = 0; j < 50; j++ )
		{ 
			ecc_ano = mean_anomaly0 + m_stSatellite.eccentricity * sin( eo );  
			eo = ecc_ano; 
		}

		double cE = cos( ecc_ano );
		double sE = sin( ecc_ano );
		double qee = sqrt( 1.0 - m_stSatellite.eccentricity * m_stSatellite.eccentricity );
		double xw = m_stSatellite.semi_major * ( cE - m_stSatellite.eccentricity );	// (3.69)
		double yw = m_stSatellite.semi_major * qee * sE;							// (3.69)
		double Edot = dfMM / ( 1.0 - m_stSatellite.eccentricity * cE );		
		double xwdot = - m_stSatellite.semi_major * Edot * sE;						// (3.70)
		double ywdot =   m_stSatellite.semi_major * Edot * qee * cE;				// (3.70)

		double sO = sin( m_stSatellite.Longi_RisingNode );
		double cO = cos( m_stSatellite.Longi_RisingNode );
		double sw = sin( m_stSatellite.perigee );
		double cw = cos( m_stSatellite.perigee );
		double si = sin( m_stSatellite.inclination );
		double ci = cos( m_stSatellite.inclination );

		double Px = cO * cw - sO * sw * ci;
		double Py = sO * cw + cO * sw * ci;										// (3.43)
		double Pz = sw * si;

		double Qx = - cO * sw - sO * cw * ci;
		double Qy = - sO * sw + cO * cw * ci;									// (3.44)
		double Qz =   cw * si;

		m_dfX = xw * Px + yw * Qx;
		m_dfY = xw * Py + yw * Qy;												// (3.57)
		m_dfZ = xw * Pz + yw * Qz;

		m_dfXdot = xwdot * Px + ywdot * Qx;
		m_dfYdot = xwdot * Py + ywdot * Qy;
		m_dfZdot = xwdot * Pz + ywdot * Qz;

		return;
	}
	catch( ... )
	{
		return;
	}
}

void cTwoBody::ComputeECEFXYZ( double dfTimeJD )
{
	ComputeInertialXYZ( dfTimeJD );
	TransformInertialToECEF( dfTimeJD );
}


void cTwoBody::TransformInertialToECEF( double dfTimeJD )
{
	try
	{
		if( m_stSatellite.GM != 1 )
		{
			m_stSatellite.semi_major = m_stSatellite.semi_major * g_dfEarthSemiMajor;

			m_dfX *= g_dfEarthSemiMajor;
			m_dfY *= g_dfEarthSemiMajor;
			m_dfZ *= g_dfEarthSemiMajor;

			m_dfXdot *= m_dfTransCoe;
			m_dfYdot *= m_dfTransCoe;
			m_dfZdot *= m_dfTransCoe;
		}

		double ra = cGreenwichST::ComputeGST( dfTimeJD );

		cCoorTrans CoorTrans;

		CoorTrans.RotatingAroundZAxis( ra, m_dfX, m_dfY, m_dfZ );	// X, Y, Z now in ECEF 

		CoorTrans.RotatingAroundZAxis( ra, m_dfXdot, m_dfYdot, m_dfZdot );	

		m_dfXdot = m_dfXdot + g_dfEarthAngVelocity * m_dfY;
		m_dfYdot = m_dfYdot - g_dfEarthAngVelocity * m_dfX;			// XYZdot now in ECEF

		if( m_stSatellite.GM != 1 )
		{
			m_dfX /= g_dfEarthSemiMajor;
			m_dfY /= g_dfEarthSemiMajor;
			m_dfZ /= g_dfEarthSemiMajor;

			m_dfXdot /= m_dfTransCoe; 
			m_dfYdot /= m_dfTransCoe;
			m_dfZdot /= m_dfTransCoe;
		}

		return;
	}
	catch( ... )
	{
		return;
	}
}


void cTwoBody::GetXYZ( double &dfX, double &dfY, double &dfZ )
{
	dfX = m_dfX;
	dfY = m_dfY;   
	dfZ = m_dfZ;

	if( m_stSatellite.GM != 1 )
	{
		dfX = m_dfX * g_dfEarthSemiMajor;   
		dfY = m_dfY * g_dfEarthSemiMajor;   
		dfZ = m_dfZ * g_dfEarthSemiMajor;
	}
}

void cTwoBody::GetXYZ( double *pdfPos )
{
	pdfPos[ 0 ] = m_dfX;   
	pdfPos[ 1 ] = m_dfY;   
	pdfPos[ 2 ] = m_dfZ;
	
	if( m_stSatellite.GM != 1 )
	{
		pdfPos[ 0 ] = m_dfX * g_dfEarthSemiMajor;   
		pdfPos[ 1 ] = m_dfY * g_dfEarthSemiMajor;   
		pdfPos[ 2 ] = m_dfZ * g_dfEarthSemiMajor;
	}
}


void cTwoBody::GetXYZdot( double &dfXdot, double &dfYdot, double &dfZdot )
{
	dfXdot = m_dfXdot;   
	dfYdot = m_dfYdot;   
	dfZdot = m_dfZdot;

	if( m_stSatellite.GM != 1 )
	{
		dfXdot = m_dfXdot * m_dfTransCoe;   
		dfYdot = m_dfYdot * m_dfTransCoe;   
		dfZdot = m_dfZdot * m_dfTransCoe;
	}
}


void cTwoBody::GetXYZdot( double *pdfVel )
{
	pdfVel[ 0 ] = m_dfXdot;   
	pdfVel[ 1 ] = m_dfYdot;   
	pdfVel[ 2 ] = m_dfZdot;
	if( m_stSatellite.GM != 1 )
	{
		pdfVel[ 0 ] = m_dfXdot / m_dfTransCoe;   
		pdfVel[ 1 ] = m_dfYdot / m_dfTransCoe;   
		pdfVel[ 2 ] = m_dfZdot / m_dfTransCoe;
	}
}


double cTwoBody::GetSemiMajor()		// in meter
{
	if( m_stSatellite.GM == 1 )	return m_stSatellite.semi_major;			// in meter	
	else 
		return ( m_stSatellite.semi_major * g_dfEarthSemiMajor );			// in meters
}


double cTwoBody::GetInclination()		// in degree
{
	return m_stSatellite.inclination * g_dfRAD2DEG;	// in degrees
}


void cTwoBody::SetPosVel( double dfX,    double dfY,    double dfZ,
		                  double dfXdot, double dfYdot, double dfZdot, 
						  double dfTimeJD )
{
	m_dfX = dfX;
	m_dfY = dfY;
	m_dfZ = dfZ;
	m_dfXdot = dfXdot;
	m_dfYdot = dfYdot;
	m_dfZdot = dfZdot;

	m_stSatellite.referencetime_jd = dfTimeJD;

	m_stSatellite.referencetime_jd_int = 
		(double)( (int) m_stSatellite.referencetime_jd );
	m_stSatellite.referencetime_jd_fra = 
		m_stSatellite.referencetime_jd - 
		m_stSatellite.referencetime_jd_int;

	if( sqrt( m_dfX * m_dfX + m_dfY * m_dfY + m_dfZ * m_dfZ ) < 100.0 ) m_stSatellite.GM = 0;
	else m_stSatellite.GM = 1;

	ComputeElementFromPosVel();
}


void cTwoBody::SetPosVel( double *pdfPos, double *pdfVel, double dfTimeJD )
{
	SetPosVel( pdfPos[ 0 ], pdfPos[ 1 ], pdfPos[ 2 ], 
		       pdfVel[ 0 ], pdfVel[ 1 ], pdfVel[ 2 ], dfTimeJD );	
}

void cTwoBody::SetPosVel( double dfX,    double dfY,    double dfZ,
	                      double dfXdot, double dfYdot, double dfZdot, 
				          double dfTimeIntJD, double dfTimeFraJD )
{
	m_dfX = dfX;
	m_dfY = dfY;
	m_dfZ = dfZ;
	m_dfXdot = dfXdot;
	m_dfYdot = dfYdot;
	m_dfZdot = dfZdot;

	m_stSatellite.referencetime_jd = dfTimeIntJD + dfTimeFraJD;

	m_stSatellite.referencetime_jd_int = dfTimeIntJD;
	m_stSatellite.referencetime_jd_fra = dfTimeFraJD;

	if( sqrt( m_dfX * m_dfX + m_dfY * m_dfY + m_dfZ * m_dfZ ) < 100.0 ) m_stSatellite.GM = 0;
	else m_stSatellite.GM = 1;

	ComputeElementFromPosVel();
}

void cTwoBody::SetPosVel( double *pdfPos, double *pdfVel, 
		                  double dfTimeIntJD, double dfTimeFraJD )
{
	SetPosVel( pdfPos[ 0 ], pdfPos[ 1 ], pdfPos[ 2 ], 
		       pdfVel[ 0 ], pdfVel[ 1 ], pdfVel[ 2 ], 
			   dfTimeIntJD, dfTimeFraJD );	
}

void cTwoBody::ComputeElementFromPosVel()
{
	try
	{
		if( m_stSatellite.GM == 1 )
		{
			m_dfX /= g_dfEarthSemiMajor;
			m_dfY /= g_dfEarthSemiMajor;
			m_dfZ /= g_dfEarthSemiMajor;

			m_dfXdot /= m_dfTransCoe;
			m_dfYdot /= m_dfTransCoe;
			m_dfZdot /= m_dfTransCoe;
		}

	//	now X,Y,Z in earth radius

		double r_value = sqrt( m_dfX * m_dfX + m_dfY * m_dfY + m_dfZ * m_dfZ );
		double vsquare = m_dfXdot * m_dfXdot + m_dfYdot * m_dfYdot + m_dfZdot * m_dfZdot;	// V^2
		double rxrdot = m_dfX * m_dfXdot + m_dfY * m_dfYdot + m_dfZ * m_dfZdot;				// r * rdot
		double rdot = rxrdot / r_value;											// rdot
			
		m_stSatellite.semi_major = 2.0 / r_value - vsquare / dfMuIOD;			// (3.153)
		m_stSatellite.semi_major = 1.0 / m_stSatellite.semi_major;		

		double ec = 1.0 - r_value / m_stSatellite.semi_major;					// (3.154)
		double es = rxrdot / sqrt( dfMuIOD * m_stSatellite.semi_major );		// (3.160)

		m_stSatellite.eccentricity = sqrt( ec * ec + es * es );					// (3.161)

		double E0 = atan2( es, ec );
		if( E0 < 0.0 ) E0 = E0 + g_dfTWOPI;

		m_stSatellite.MeanAnomaly = E0 - es;									// Kepler equ

		double U0[ 3 ] = { m_dfX / r_value, m_dfY / r_value, m_dfZ / r_value  };// (3.162)

		double rvdot = sqrt( vsquare - rdot * rdot );							// (3.132)
		double V0[ 3 ];

		V0[ 0 ] = ( r_value * m_dfXdot - rdot * m_dfX ) / ( r_value * rvdot );
		V0[ 1 ] = ( r_value * m_dfYdot - rdot * m_dfY ) / ( r_value * rvdot );
		V0[ 2 ] = ( r_value * m_dfZdot - rdot * m_dfZ ) / ( r_value * rvdot );	// (3.165)

		double si = sqrt( U0[ 2 ] * U0[ 2 ] + V0[ 2 ] * V0[ 2 ] );			// (3.195)
		double ci = sqrt( ( U0[ 0 ] + V0[ 1 ] ) * ( U0[ 0 ] + V0[ 1 ] ) +
						( U0[ 1 ] - V0[ 0 ] ) * ( U0[ 1 ] - V0[ 0 ] )   ) -
					1.0;													// (3.196)
		m_stSatellite.inclination = atan2( si, ci );

		double cl = ( U0[ 0 ] + V0[ 1 ] ) / ( 1.0 + ci );					// (3.197)
		double sl = ( U0[ 1 ] - V0[ 0 ] ) / ( 1.0 + ci );					// (3.198)

		double su = U0[ 2 ] / si;											// (3.200)
		double cu = V0[ 2 ] / si;											// (3.200)

		m_stSatellite.Longi_RisingNode = atan2( sl, cl ) - atan2( su, cu );	// (3.199)
		if( m_stSatellite.Longi_RisingNode < 0.0 ) m_stSatellite.Longi_RisingNode += g_dfTWOPI;

		double p = m_stSatellite.semi_major * 
			       ( 1.0 - m_stSatellite.eccentricity * m_stSatellite.eccentricity );	//(3.166)
		double cv = ( p / r_value - 1.0 ) / m_stSatellite.eccentricity;					// (3.204)
		double sv = rdot * sqrt( p / dfMuIOD ) / m_stSatellite.eccentricity;			// (3.205)
		
		m_stSatellite.perigee = atan2( su, cu ) - atan2( sv, cv );						// (3.206)
		if( m_stSatellite.perigee < 0.0 ) m_stSatellite.perigee += g_dfTWOPI;

		if( m_stSatellite.GM == 1 )
		{
			m_stSatellite.semi_major *= g_dfEarthSemiMajor;

			m_dfX *= g_dfEarthSemiMajor;
			m_dfY *= g_dfEarthSemiMajor;
			m_dfZ *= g_dfEarthSemiMajor;

			m_dfXdot *= m_dfTransCoe;
			m_dfYdot *= m_dfTransCoe;
			m_dfZdot *= m_dfTransCoe;
		}

		return;
	}
	catch( ... )
	{
		return;
	}
}

void cTwoBody::GetElements( double &dfSM, double &dfEcc, double &dfInc,
		                    double &dfRAAN, double &dfPerigee, double &dfMA )
{
    dfSM  = m_stSatellite.semi_major;
    dfEcc  = m_stSatellite.eccentricity;
    dfInc  = m_stSatellite.inclination;
    dfPerigee = m_stSatellite.perigee;
    dfRAAN  = m_stSatellite.Longi_RisingNode;
    dfMA = m_stSatellite.MeanAnomaly;

	if( m_stSatellite.GM == 0 )	dfSM *= g_dfEarthSemiMajor;
}


/************************************************************************************************

  Fortran subroutine: compute_kepler_elements_from_position_velocity

  Compute the orbital elements from given position, velocity and GM.

  pdfPos in m
  pdfVel in m/s
  dfGM in m^3/s^2

  pdfElements in m, radian, in the order of (a, e, i, o, w, M )
  
  bPartial is true if the partial derivatives are required

  pdfPartial are the partial derivatived of the orbital elements WRT to the pos/vel
    PKPX(I,J) = D K(I) / D X(J), pdfPartial is the vectorized PKPX by row

************************************************************************************************/

BOOL cTwoBody::ComputeElementsFromPosVel( double *pdfPos, double *pdfVel, double dfGM,
		                                  double *pdfElements, 
									 	  BOOL bPartial, double *pdfPartial )
{
	try
	{
		double RRDOT = g_Matrix.Vec_Multiply_Vec( pdfPos, pdfVel, 3 );
		double R = sqrt( g_Matrix.Vec_Multiply_Vec( pdfPos, pdfPos, 3 ) );
		double VSQ = g_Matrix.Vec_Multiply_Vec( pdfVel, pdfVel, 3 );
		pdfElements[ 0 ] = 2.0 / R - VSQ / dfGM;

		pdfElements[ 0 ] = 1.0 / pdfElements[ 0 ]; // semi-major axis

		double HX = pdfPos[ 1 ] * pdfVel[ 2 ] - pdfPos[ 2 ] * pdfVel[ 1 ];
		double HY = pdfPos[ 2 ] * pdfVel[ 0 ] - pdfPos[ 0 ] * pdfVel[ 2 ];
		double HZ = pdfPos[ 0 ] * pdfVel[ 1 ] - pdfPos[ 1 ] * pdfVel[ 0 ];

		double HSINI2 = HX * HX + HY * HY;
		double HSQ = HSINI2 + HZ * HZ;
		double H = sqrt( HSQ );
		double OME2 = HSQ / pdfElements[ 0 ] / dfGM;
		double ESQ = 1.0 - OME2;
		pdfElements[ 1 ] = sqrt( ESQ );	// ecentricity
		double TEMPE = pdfElements[ 1 ];
	 
		double COSI = HZ / H;
		double SINI2 = 1.0 - COSI * COSI;
		double SINI = sqrt( SINI2 );
		double HSINI = H * SINI;
		pdfElements[ 2 ] = atan2( SINI, COSI );	// inclination
		if( pdfElements[ 2 ] < 0.0 ) pdfElements[ 2 ] += g_dfTWOPI;

		pdfElements[ 3 ] = atan2( HX, -HY );
		if( pdfElements[ 3 ] < 0.0 ) pdfElements[ 3 ] += g_dfTWOPI;	//RA of ascending node

		double SUPROD = -HZ * ( pdfPos[ 0 ] * HX + pdfPos[ 1 ] * HY ) + 
			            pdfPos[ 2 ] * HSINI2;
		double CUPROD = H * ( -pdfPos[ 0 ] * HY + pdfPos[ 1 ] * HX );
		double RESINF = pdfElements[ 0 ] * OME2 * RRDOT / H;
		double RECOSF = pdfElements[ 0 ] * OME2 - R;
		double dfTRUE = atan2( RESINF, RECOSF );	// true anomaly
		if( dfTRUE < 0.0 ) dfTRUE += g_dfTWOPI;
	      
		double SURECF = SUPROD * RECOSF;
		double CURECF = CUPROD * RECOSF;
		double CURESF = CUPROD * RESINF;
		double SURESF = SUPROD * RESINF;
		double SWPROD = ( SURECF - CURESF );
		double CWPROD = ( CURECF + SURESF );
		pdfElements[ 4 ] = atan2( SWPROD, CWPROD );	// argument of perigee
		if( pdfElements[ 4 ] < 0.0 ) pdfElements[ 4 ] += g_dfTWOPI;
	      
		double RTOME2 = sqrt( fabs( OME2 ) );
		double AESINE = RESINF / RTOME2;
		double AECOSE = pdfElements[ 0 ] - R;

		double ECC = atan2( AESINE, AECOSE );		// eccentric anomaly
		if( ECC < 0.0 ) ECC += g_dfTWOPI;

		pdfElements[ 5 ] = ECC - AESINE / pdfElements[ 0 ];	// mean anomaly
		if( pdfElements[ 5 ] < 0.0 ) pdfElements[ 5 ] += g_dfTWOPI;
		if( pdfElements[ 5 ] >= g_dfTWOPI ) 
			pdfElements[ 5 ] = fmod( pdfElements[ 3 ], g_dfTWOPI );

		if( !bPartial ) return TRUE;

		for( int i = 0; i < 36; i++ ) pdfPartial[ i ] = 0;

		computePartials( pdfPos, pdfVel, dfGM, pdfElements, true, pdfPartial );

		// COMPUTE PARTIALS OF KEPLER ELEMENTS W/R TO CARTESIAN ELEMENTS
/*
		double CUXYZP[ 4 ], PH[ 4 ][ 7 ], PKPX[ 7 ][ 7 ];

		double C1 = 2.0 *  pdfElements[ 0 ] * pdfElements[ 0 ];
		double CA1 = C1 / R / R / R;
		double CA2 = C1 / dfGM;
		double CH1 = HX / H;
		double CH2 = HY / H;
		double CH3 = HZ / H;
		double EPROD = OME2 / TEMPE;
		double CE1 = EPROD / ( pdfElements[ 0 ] + pdfElements[ 0 ] );
		double CE2 = -EPROD / H;
		double CN1 = -HY / HSINI2;
		double CN2 = HX / HSINI2;
		double CU5 = R * HSQ * SINI;
		double COSU = CUPROD / CU5;
		double HCOSU = H * COSU;
		double HZCOSU = HZ * COSU;
		double SINU = SUPROD / CU5;
		double HSINU = H * SINU;
		double CU1 = -pdfPos[ 1 ] * HSINU - pdfPos[ 0 ] * HZCOSU;
		double CU2 = pdfPos[ 0 ] * HSINU - pdfPos[ 1 ] * HZCOSU;
		double CU3 = 
			-( pdfPos[ 0 ] * HX + pdfPos[ 1 ] * HY + pdfPos[ 2 ] * HZ ) * COSU - 
			pdfPos[ 2 ] * HZCOSU;
		double CU4 = 2.0 * pdfPos[ 2 ] * HCOSU + 
			         ( pdfPos[ 0 ] * HY - pdfPos[ 1 ] * HX ) * SINU;
		CUXYZP[ 1 ] = HY * HSINU - HX * HZCOSU;
		CUXYZP[ 2 ] = -HX * HSINU - HY * HZCOSU;
		CUXYZP[ 3 ] = H * HCOSU - HZ * HZCOSU;
		double ADR = pdfElements[ 0 ] / R;
		double RE = R * TEMPE;
		double COSF = RECOSF / RE;
		double SINF = RESINF / RE;
		C1 = COSF * RRDOT / H;
		double C2 = C1 - SINF;
		double CF1 = OME2 * C2 / RE;
		double CF2 = -( C2 + C2 ) * ADR;
		C2 = -CE2 * ADR;
		double CF3 = -C2 * C1;
		double CF4 = C2 * COSF;
		double CF5 = SINF / ( RE * R );
		double RDA = R / pdfElements[ 0 ];
		double CM1 = RDA * RDA / RTOME2;
		double CM2 = -( RDA / OME2 + 1.0 ) * AESINE / ( pdfElements[ 0 ] * TEMPE );
		
		PH[ 1 ][ 1 ] = 0.0;
		PH[ 2 ][ 1 ] = -pdfVel[ 2 ];
		PH[ 3 ][ 1 ] = pdfVel[ 1 ];
		PH[ 1 ][ 2 ] = pdfVel[ 2 ];
		PH[ 2 ][ 2 ] = 0.0;
		PH[ 3 ][ 2 ] = -pdfVel[ 0 ];
		PH[ 1 ][ 3 ] = -pdfVel[ 1 ];
		PH[ 2 ][ 3 ] = pdfVel[ 0 ];
		PH[ 3 ][ 3 ] = 0.0;
		PH[ 1 ][ 4 ] = 0.0;
		PH[ 2 ][ 4 ] = pdfPos[ 2 ];
		PH[ 3 ][ 4 ] = -pdfPos[ 1 ];
		PH[ 1 ][ 5 ] = -pdfPos[ 2 ];
		PH[ 2 ][ 5 ] = 0.0;
		PH[ 3 ][ 5 ] = pdfPos[ 0 ];
		PH[ 1 ][ 6 ] = pdfPos[ 1 ];
		PH[ 2 ][ 6 ] = -pdfPos[ 0 ];
		PH[ 3 ][ 6 ] = 0.0;
	    
		for( int I = 1; I <= 3; I++ )
		{
			PKPX[ 1 ][ I ] = CA1 * pdfPos[ I - 1 ];
			PKPX[ 1 ][ I + 3 ]= CA2 * pdfVel[ I - 1 ];
			double PPFS = CF5 * pdfPos[ I - 1 ];
			double PPUS = CUXYZP[ I ];
					
			for( int J = I; J <= 6; J = J + 3 )
			{
				double PHPS = 
					CH1 * PH[ 1 ][ J ] + CH2 * PH[ 2 ][ J ] + CH3 * PH[ 3 ][ J ];
				PKPX[ 2 ][ J ] = CE1 * PKPX[ 1 ][ J ] + CE2 * PHPS;
				PKPX[ 3 ][ J ] = ( ( COSI * PHPS - PH[ 3 ][ J ] ) / HSINI );
				PKPX[ 4 ][ J ] = ( CN1 * PH[ 1 ][ J ] + CN2 * PH[ 2 ][ J ] );
				double PUPS = ( CU1 * PH[ 1 ][ J ] + CU2 * PH[ 2 ][ J ] + 
								CU3 * PH[ 3 ][ J ] + CU4 * PHPS + PPUS ) / CU5;
				double PFPS = 
					CF1 * PKPX[ 1 ][ J ] + CF2 * PKPX[ 2 ][ J ] + CF3 * PHPS + PPFS;
				double temp = J <= 3 ? pdfPos[ I - 1 ] : pdfVel[ I - 1 ];
				PFPS += CF4 * temp;
				PKPX[ 5 ][ J ] = ( PUPS - PFPS );
				PKPX[ 6 ][ J ] = ( CM1 * PFPS + CM2 * PKPX[ 2 ][ J ] );
				PPFS = 0.0;
				PPUS = 0.0;	
			}
		}

		int k = 0;
		for( int I = 1; I <= 6; I++ )
		{
			for( int J = 1; J <= 6; J++ )
			{
				pdfPartial[ k ] = PKPX[ J ][ I ];
				k++;
			}
		}

		k = 0;
		for( int i = 0; i < 6; i++ )
		{
			for( int j = 0; j < 6; j++ )
			{
				printf( "%20.10f ", pdfPartial[ k ] );
				k++;
			}	
			printf( "\n" );
		}
*/
/*
		double JD = 1.0, element[ 6 ], pos[ 3 ], vel[ 3 ], partial[ 36 ];

		for( int i = 0; i < 6; i++ ) element[ i ] = pdfElements[ i ];

		element[ 2 ] += 0.00000001;
		ComputePosVelFromElements( JD, JD, element, g_dfEarthGM, pos, vel );
		int k = 12;
		for( int i = 0; i < 3; i++ )
		{
			pdfPartial[ k + i ]     = 0.00000001 / ( pos[ i ] - pdfPos[ i ] );
			pdfPartial[ k + i + 3 ] = 0.00000001 / ( vel[ i ] - pdfVel[ i ] );
		}
		element[ 2 ] -= 0.00000001;

		element[ 3 ] += 0.00000001;
		ComputePosVelFromElements( JD, JD, element, g_dfEarthGM, pos, vel );
		k = 18;
		for( int i = 0; i < 3; i++ )
		{
			pdfPartial[ k + i ]     = 0.00000001 / ( pos[ i ] - pdfPos[ i ] );
			pdfPartial[ k + i + 3 ] = 0.00000001 / ( vel[ i ] - pdfVel[ i ] );
		}
		element[ 3 ] -= 0.00000001;
*/
		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}

}

/************************************************************************************************

  Fortran subroutine: compute_state_from_orbit_element

  Compute the position, velocity from the given orbital elements and GM.

  dfJDRef:   reference epoch of the orbital elements
  dfJD:      epoch of interest at which the pos/vel is required
  pdfElements in m, radian, in the order of (a, e, i, o, w, M )
  dfGM in m^3/s^2

  pdfPos in m
  pdfVel in m/s  
  
************************************************************************************************/

BOOL cTwoBody::ComputePosVelFromElements( double dfJDRef, double dfJD,
										  double *pdfElements, double dfGM, 
		                                  double *pdfPos, double *pdfVel )
{
	try
	{
		double dfMM, dfDeltaT;

		dfMM = sqrt( dfGM / pow( pdfElements[ 0 ], 3.0 ) );
		
		dfDeltaT = ( dfJD - dfJDRef ) * 86400.0;

		double mean_anomaly0 = pdfElements[ 5 ] + dfMM * dfDeltaT;
		double eo = mean_anomaly0;

		double ecc_ano;
		for( int j = 0; j < 50; j++ )
		{ 
			ecc_ano = mean_anomaly0 + pdfElements[ 1 ] * sin( eo );  
			eo = ecc_ano; 
		}

		double cE = cos( ecc_ano );
		double sE = sin( ecc_ano );
		double qee = sqrt( 1.0 - pdfElements[ 1 ] * pdfElements[ 1 ] );
		double xw = pdfElements[ 0 ] * ( cE - pdfElements[ 1 ] );			// (3.69)
		double yw = pdfElements[ 0 ] * qee * sE;							// (3.69)
		double Edot = dfMM / ( 1.0 - pdfElements[ 1 ] * cE );		
		double xwdot = - pdfElements[ 0 ] * Edot * sE;						// (3.70)
		double ywdot =   pdfElements[ 0 ] * Edot * qee * cE;				// (3.70)

		double sO = sin( pdfElements[ 3 ] );
		double cO = cos( pdfElements[ 3 ] );
		double sw = sin( pdfElements[ 4 ] );
		double cw = cos( pdfElements[ 4 ] );
		double si = sin( pdfElements[ 2 ] );
		double ci = cos( pdfElements[ 2 ] );

		double Px = cO * cw - sO * sw * ci;
		double Py = sO * cw + cO * sw * ci;									// (3.43)
		double Pz = sw * si;

		double Qx = - cO * sw - sO * cw * ci;
		double Qy = - sO * sw + cO * cw * ci;								// (3.44)
		double Qz =   cw * si;

		pdfPos[ 0 ] = xw * Px + yw * Qx;
		pdfPos[ 1 ] = xw * Py + yw * Qy;									// (3.57)
		pdfPos[ 2 ] = xw * Pz + yw * Qz;

		pdfVel[ 0 ] = xwdot * Px + ywdot * Qx;
		pdfVel[ 1 ] = xwdot * Py + ywdot * Qy;
		pdfVel[ 2 ] = xwdot * Pz + ywdot * Qz;

		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}
}

/************************************************************************************************

  Conversion between the true and mean anomaly

  a is the semi major axis in meter
  e is the eccentricity

  M is the mean anomaly in radian
  f is the true anomaly in radian
  E is the eccentric anomaly in radian
   
************************************************************************************************/

BOOL cTwoBody::MeanAnomaly2TrueAnomaly( double a, double e, double &M, 
									    double &f, double &E )
{
	try
	{
		E = M;  // eccentric anomaly

		for( int k = 0; k < 50; k++ )
		{
			E = M + e * sin( E );
		}
	     
		// compute f (true anomaly) from E

		double rcosf = a * ( cos( E ) - e );					// Eq 2.26 Liu Lin
		double rsinf = a * sqrt( 1.0 - e * e ) * sin( E );		// Eq 2.27

		f = atan2( rsinf, rcosf );
		if( f < 0.0 ) f += g_dfTWOPI;

		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}
}

/******************************************************************************

  Compute the mean motion given the semi major axis and GM

******************************************************************************/

double cTwoBody::ComputeMeanMotion( double a, double dfGM )
{
	double n;

	n = sqrt( dfGM / a / a / a );

	return n;
}

/******************************************************************************

  the delta V is expressed as delta xdot and delta ydot in the orbital plane

******************************************************************************/
BOOL cTwoBody::deltaVEffect( double dxdot, double dydot, char *fileName )
{
	double dfMM, dfDeltaT;

	if( m_stSatellite.GM == 1 ) 
		dfMM = sqrt( g_dfEarthGM / pow( m_stSatellite.semi_major, 3.0 ) );
	else dfMM = sqrt( dfMuIOD / m_stSatellite.semi_major ) / 
			    m_stSatellite.semi_major;	 		

	dfDeltaT = 0;

	double mean_anomaly0 = m_stSatellite.MeanAnomaly + dfMM * dfDeltaT;
	double eo = mean_anomaly0;

	double ecc_ano;
	for( int j = 0; j < 50; j++ )
	{ 
		ecc_ano = mean_anomaly0 + m_stSatellite.eccentricity * sin( eo );  
		eo = ecc_ano; 
	}

	double cE = cos( ecc_ano );
	double sE = sin( ecc_ano );
	double qee = sqrt( 1.0 - m_stSatellite.eccentricity * m_stSatellite.eccentricity );
	double xw = m_stSatellite.semi_major * ( cE - m_stSatellite.eccentricity );	// (3.69)
	double yw = m_stSatellite.semi_major * qee * sE;							// (3.69)
	double Edot = dfMM / ( 1.0 - m_stSatellite.eccentricity * cE );		
	double xwdot = - m_stSatellite.semi_major * Edot * sE;						// (3.70)
	double ywdot =   m_stSatellite.semi_major * Edot * qee * cE;				// (3.70)

	FILE *file = NULL;
	
	fopen_s( &file, fileName, "w" );

	fprintf( file, "before\n" );
	fprintf( file, "a %20.10f\n", m_stSatellite.semi_major );
	fprintf( file, "e %20.10f\n", m_stSatellite.eccentricity );
	fprintf( file, "i %20.10f\n", m_stSatellite.inclination * g_dfRAD2DEG );
	fprintf( file, "O %20.10f\n", m_stSatellite.Longi_RisingNode * g_dfRAD2DEG );
	fprintf( file, "w %20.10f\n", m_stSatellite.perigee * g_dfRAD2DEG );
	fprintf( file, "M %20.10f\n", m_stSatellite.MeanAnomaly * g_dfRAD2DEG );
	
	fprintf( file, "x  %19.9f y  %19.9f r %19.9f\n", xw, yw, sqrt( xw * xw + yw * yw ) );
	fprintf( file, "xv %19.9f yv %19.9f v %19.9f\n", xwdot, ywdot, sqrt( xwdot * xwdot + ywdot * ywdot ) );

	fprintf( file, "dxv %18.8f dyv %18.8f\n", dxdot, dydot );

	xwdot += dxdot;
	ywdot += dydot;

	double sO = sin( m_stSatellite.Longi_RisingNode );
	double cO = cos( m_stSatellite.Longi_RisingNode );
	double sw = sin( m_stSatellite.perigee );
	double cw = cos( m_stSatellite.perigee );
	double si = sin( m_stSatellite.inclination );
	double ci = cos( m_stSatellite.inclination );

	double Px = cO * cw - sO * sw * ci;
	double Py = sO * cw + cO * sw * ci;										// (3.43)
	double Pz = sw * si;

	double Qx = - cO * sw - sO * cw * ci;
	double Qy = - sO * sw + cO * cw * ci;									// (3.44)
	double Qz =   cw * si;

	m_dfX = xw * Px + yw * Qx;
	m_dfY = xw * Py + yw * Qy;												// (3.57)
	m_dfZ = xw * Pz + yw * Qz;

	m_dfXdot = xwdot * Px + ywdot * Qx;
	m_dfYdot = xwdot * Py + ywdot * Qy;
	m_dfZdot = xwdot * Pz + ywdot * Qz;

	ComputeElementFromPosVel();

	fprintf( file, "after\n" );
	fprintf( file, "a %20.10f\n", m_stSatellite.semi_major );
	fprintf( file, "e %20.10f\n", m_stSatellite.eccentricity );
	fprintf( file, "i %20.10f\n", m_stSatellite.inclination * g_dfRAD2DEG );
	fprintf( file, "O %20.10f\n", m_stSatellite.Longi_RisingNode * g_dfRAD2DEG );
	fprintf( file, "w %20.10f\n", m_stSatellite.perigee * g_dfRAD2DEG );
	fprintf( file, "M %20.10f\n", m_stSatellite.MeanAnomaly * g_dfRAD2DEG );
	
	fprintf( file, "x  %19.9f y  %19.9f r %19.9f\n", xw, yw, sqrt( xw * xw + yw * yw ) );
	fprintf( file, "xv %19.9f yv %19.9f v %19.9f\n", xwdot, ywdot, sqrt( xwdot * xwdot + ywdot * ywdot ) );

	return TRUE;
}


// this function follows the Eqs in Montenbruck's Satellite Orbits: Models, Methods, Applications (Sections 7.1.2 and 7.1.3)
// element: a, e, i, O, w, M
// E2X: true if the derivatives of the elements WRT the pos/vel are asked
//      false if the debrivatives of the pos/vel WRT the element are asked
// partial is a 6 x 6 matrix
//   if E2X is true, then the first row is the derivatives of a WRT pos/vel, ....
//   if E2X is false, the first row is the derivatives of X WRT the element 
bool cTwoBody::computePartials( double *pos, double *vel, double GM, double *element, bool E2X, double *partial )
{
	if( element[ 0 ] < 1.0 ) ComputeElementsFromPosVel( pos, vel, GM, element );

	if( fabs( pos[ 0 ] ) + fabs( pos[ 1 ] ) + fabs( pos[ 2 ] ) < 1.0 ) ComputePosVelFromElements( 1.0, 1.0, element, GM, pos, vel );

	double sO = sin( element[ 3 ] );
	double cO = cos( element[ 3 ] );
	double sw = sin( element[ 4 ] );
	double cw = cos( element[ 4 ] );
	double si = sin( element[ 2 ] );
	double ci = cos( element[ 2 ] );

	double P[ 3 ] = { cO * cw - sO * sw * ci,  sO * cw + cO * sw * ci, sw * si };
	double Q[ 3 ] = {-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si };

	double f, E, meanMotion = ComputeMeanMotion( element[ 0 ], GM );

	MeanAnomaly2TrueAnomaly( element[ 0 ], element[ 1 ], element[ 5 ], f, E );

	double x, y, xdot, ydot, cE = cos( E ), sE = sin( E ), Oe2 = 1.0 - element[ 1 ] * element[ 1 ], s1e2 = sqrt( Oe2 ), r = element[ 0 ] * ( 1 - element[ 1 ] * cE );
	
	x = element[ 0 ] * ( cE - element[  1 ] );
	y = element[ 0 ] * s1e2 * sE;

	xdot = -sqrt( GM * element[ 0 ] ) * sE / r;
	ydot =  sqrt( GM * element[ 0 ] ) * s1e2 * cE / r;

	// 7.24
	double xa = x / element[ 0 ];
	double xe = -element[ 0 ] - y * y / r / Oe2;
	double xM = xdot / meanMotion;

	double ya = y / element[ 0 ];
	double ye = x * y / r / Oe2;
	double yM = ydot / meanMotion;

	// 7.25
	double xdota = -xdot / 2 / element[ 0 ];
	double xdote = xdot * ( element[ 0 ] * element[ 0 ] / r / r ) * ( 2 * x / element[ 0 ] + element[ 1 ] / Oe2 * y * y / element[ 0 ] / element[ 0 ] );
	double xdotM = -meanMotion * element[ 0 ] * element[ 0 ] * element[ 0 ] * x / r / r / r;

	double ydota = -ydot / 2 / element[ 0 ];
	double ydote = meanMotion / s1e2 * element[ 0 ] * element[ 0 ] / r / r * ( x * x / r - y * y / element[ 0 ] / Oe2 );
	double ydotM = -meanMotion * element[ 0 ] * element[ 0 ] * element[ 0 ] * y / r / r / r;

	// 7.13
	partial[  0 ] = xa * P[ 0 ] + ya * Q[ 0 ];			// dX / da
	partial[  6 ] = xa * P[ 1 ] + ya * Q[ 1 ];			// dY / da
	partial[ 12 ] = xa * P[ 2 ] + ya * Q[ 2 ];			// dZ / da
	partial[ 18 ] = xdota * P[ 0 ] + ydota * Q[ 0 ];	// dXdot / da
	partial[ 24 ] = xdota * P[ 1 ] + ydota * Q[ 1 ];	// dXdot / da
	partial[ 30 ] = xdota * P[ 2 ] + ydota * Q[ 2 ];	// dXdot / da

	partial[  1 ] = xe * P[ 0 ] + ye * Q[ 0 ];			// dX / de
	partial[  7 ] = xe * P[ 1 ] + ye * Q[ 1 ];			// dY / de
	partial[ 13 ] = xe * P[ 2 ] + ye * Q[ 2 ];			// dZ / de
	partial[ 19 ] = xdote * P[ 0 ] + ydote * Q[ 0 ];	// dXdot / de
	partial[ 25 ] = xdote * P[ 1 ] + ydote * Q[ 1 ];	// dXdot / de
	partial[ 31 ] = xdote * P[ 2 ] + ydote * Q[ 2 ];	// dXdot / de

	partial[  5 ] = xM * P[ 0 ] + yM * Q[ 0 ];			// dX / dM
	partial[ 11 ] = xM * P[ 1 ] + yM * Q[ 1 ];			// dY / dM
	partial[ 17 ] = xM * P[ 2 ] + yM * Q[ 2 ];			// dZ / dM
	partial[ 23 ] = xdotM * P[ 0 ] + ydotM * Q[ 0 ];	// dXdot / dM
	partial[ 29 ] = xdotM * P[ 1 ] + ydotM * Q[ 1 ];	// dXdot / dM
	partial[ 35 ] = xdotM * P[ 2 ] + ydotM * Q[ 2 ];	// dXdot / dM

	double n[ 3 ] = { cO, sO, 0 }, ez[ 3 ] = { 0, 0, 1 }, W[ 3 ];

	VectorProduct( P, Q, W );

	double P2i[ 3 ], P2O[ 3 ], P2w[ 3 ]; 
	double Q2i[ 3 ], Q2O[ 3 ], Q2w[ 3 ];

	// 7.26
	VectorProduct( n,  P, P2i );
	VectorProduct( ez, P, P2O );
	VectorProduct( W,  P, P2w );

	VectorProduct( n,  Q, Q2i );
	VectorProduct( ez, Q, Q2O );
	VectorProduct( W,  Q, Q2w );

	// 7.14
	partial[  2 ] = x * P2i[ 0 ] + y * Q2i[ 0 ];
	partial[  8 ] = x * P2i[ 1 ] + y * Q2i[ 1 ];
	partial[ 14 ] = x * P2i[ 2 ] + y * Q2i[ 2 ];
	partial[ 20 ] = xdot * P2i[ 0 ] + ydot * Q2i[ 0 ];
	partial[ 26 ] = xdot * P2i[ 1 ] + ydot * Q2i[ 1 ];
	partial[ 32 ] = xdot * P2i[ 2 ] + ydot * Q2i[ 2 ];

	partial[  3 ] = x * P2O[ 0 ] + y * Q2O[ 0 ];
	partial[  9 ] = x * P2O[ 1 ] + y * Q2O[ 1 ];
	partial[ 15 ] = x * P2O[ 2 ] + y * Q2O[ 2 ];
	partial[ 21 ] = xdot * P2O[ 0 ] + ydot * Q2O[ 0 ];
	partial[ 27 ] = xdot * P2O[ 1 ] + ydot * Q2O[ 1 ];
	partial[ 33 ] = xdot * P2O[ 2 ] + ydot * Q2O[ 2 ];

	partial[  4 ] = x * P2w[ 0 ] + y * Q2w[ 0 ];
	partial[ 10 ] = x * P2w[ 1 ] + y * Q2w[ 1 ];
	partial[ 16 ] = x * P2w[ 2 ] + y * Q2w[ 2 ];
	partial[ 22 ] = xdot * P2w[ 0 ] + ydot * Q2w[ 0 ];
	partial[ 28 ] = xdot * P2w[ 1 ] + ydot * Q2w[ 1 ];
	partial[ 34 ] = xdot * P2w[ 2 ] + ydot * Q2w[ 2 ];

	if( E2X )
	{
		g_Matrix.Mat_Inverse( partial, 6 );
	}

	return true;
}

void cTwoBody::VectorProduct( double *pdfVec1, double *pdfVec2, double *pdfVec3 )
{
	pdfVec3[ 0 ] = pdfVec1[ 1 ] * pdfVec2[ 2 ] - pdfVec1[ 2 ] * pdfVec2[ 1 ];
	pdfVec3[ 1 ] = pdfVec1[ 2 ] * pdfVec2[ 0 ] - pdfVec1[ 0 ] * pdfVec2[ 2 ];
	pdfVec3[ 2 ] = pdfVec1[ 0 ] * pdfVec2[ 1 ] - pdfVec1[ 1 ] * pdfVec2[ 0 ];
}
