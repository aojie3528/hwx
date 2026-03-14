/*************************************************************************** 

 This class is basically a translation of the FORTRAN program TLE2IRV

 Author: Jizhang Sang, Feb-Apr, 2002

 (C) COPYRIGHT -Electro Optic Systems Australia (EOS) 2002,
 All rights reserved. No part of this program may be photocopied, 
 reproduced or translated to another program language without the prior
 written consent of EOS

 Reference: SGP4 Algorithm

***************************************************************************/

/******************************************************************************************
*******************************************************************************************
*******************************************************************************************

  EXTREME CARE SHOULD BE TAKEN TO MAKE ANY MODIFICATION IN THIS CODING OF TLE TO IRV 
  CONVERSION. PLEASE CONTACT JIZHANG SANG IF NECESSARY.

*******************************************************************************************
*******************************************************************************************
*******************************************************************************************/
// #include "stdafx.h"
#undef UNICODE

#include <math.h>
#include "TLE2PosVel.h"
#include "GreenwichSiderealTime.h"
#include "DateTimeZ.h"
#include "constant.h"


cTLE2PosVel::cTLE2PosVel()
{
	m_bInit = FALSE;

	m_dfEarthRadius = 6378135.0; 
	m_dfMeanEarthRadius = 6371.0;
	m_dfEarthFlattening = 1.0 / 298.257;
	m_dfEarthRotation = 7.29211585494e-5;

	g_DateTimeZ.DateTime2JD( 1950, 1, 0, 0, 0, 0.0, m_dfJD1950 );

	// FK5 CONSTANT VALUE
    m_dfC1 = 1.72027916940703639e-2;
    m_dfC1P2P = m_dfC1 + g_dfTWOPI;
    m_dfTHGR70 = 1.7321343856509374;
	m_dfFK5R = 5.07551419432269442e-15;
	
	// WGS-72 EARTH CONSTANTS                                              
    m_dfEarthGM = 7.43669161331734132e-2;	// earth gravitational constant
    m_dfEarthRadiusKM = 6378.135;	
	m_dfXPDOTP = 229.1831180523293;	// 1440 / (2 * PI) 
    m_dfEarthRotationPerMinute = 4.37526908801129966e-3;	//  ! earth rot vel in radian/minute
    m_dfJ2 = 1.082616e-3;                                                       
    m_dfJ3 =-2.53881e-6;                                                       
    m_dfJ4 =-1.65597e-6;                                                       
    m_dfJ3OJ2 = m_dfJ3 / m_dfJ2;                                                        

	// change velocity unit from earth_radius / minute to m/s
	m_dfVelocityChange = m_dfEarthRadiusKM * m_dfEarthGM / 60.0 * 1000.0;

	// OTHER CONSTANTS *                                                     
	m_dfMinutesPerDay = 1440.0;
    m_dfMinSatGeoDis = 78.0 / m_dfEarthRadiusKM + 1.0;	// pos min dis from earth cen to sv
    m_dfQZMS2T = pow( ( 120.0 - 78.0 ) / m_dfEarthRadiusKM, 4.0 ); // ?

	m_df2O3 = 2.0 / 3.0;

	// constant used in Deep Space Effect
	m_dfZNS = 1.19459e-05;
	m_dfZES = 0.01675;                               
	m_dfZNL = 1.5835218e-04;
	m_dfZEL = .05490;                                                                          

	m_nCounter=0;
}


cTLE2PosVel::~cTLE2PosVel()
{

}

/******************************************************************************************

  Read the TLE of the given satellite ID from the given file

  each set of TLE is given in 3 lines

*******************************************************************************************/

BOOL cTLE2PosVel::ReadTLE( int nNORADID, char *szFileName )
{
	// check whether the given file exists
	WIN32_FIND_DATA FindFileData;

	HANDLE hFile = FindFirstFile( szFileName, &FindFileData );

	if( hFile == INVALID_HANDLE_VALUE ) return FALSE;

	FindClose( hFile );

	string line;
	ifstream in( szFileName );

	stSatelliteIOE tIOE;

	while( getline( in, line ) )	// 
	{
		//if( line[ 0 ] == '1' )
		string temp = line.substr(0, 2);
		if (strcmp(temp.c_str(), "1 ") == 0)
		{
			if( !ReadTLELine1( line, tIOE ) ) return FALSE;
			getline( in, line );
			if( tIOE.nSatelliteID == nNORADID )
			{
				if( !ReadTLELine2( line, tIOE ) ) return FALSE;
				if( !SetOrbitalElements( tIOE ) ) return FALSE;
				return TRUE;
			}
		}
	}

	return FALSE;
}

/******************************************************************************************

  Set the TLE of the satellite and perform the initialisation of the algorithm

*******************************************************************************************/

BOOL cTLE2PosVel::SetOrbitalElements( stSatelliteIOE stIOE )
{
	m_bInit = FALSE;

	m_pdfTLE[ 3 ] = m_dfBSTAR = stIOE.pfElement7to18[ 2 ];
	m_dfRefJD = stIOE.GetRefJD();

	int nDoY, y, m, d, h, mm;
	double s;

	g_DateTimeZ.ComputeDayOfYearFromJD( m_dfRefJD, nDoY );
	g_DateTimeZ.JD2DateTime( m_dfRefJD, y, m, d, h, mm, s );
	if( y >= 2000 ) y -= 2000;
	else y -= 1900;

	m_pdfTLE[ 1 ] = (double)y;
	m_pdfTLE[ 2 ] = (double)nDoY + h / 24.0 + mm / 1440.0 + s / 86400.0;

//	m_dfRefJD = (double) stSatelliteIOE.nIntJD + (double) stSatelliteIOE.nFractionJD * 1.0e-9;
	m_pdfTLE[ 9 ] = m_dfMM0 = (double) stIOE.pnElement1to6[ 0 ] * 1.0e-8;
	m_pdfTLE[ 4 ] = m_dfEcc0 = (double) stIOE.pnElement1to6[ 1 ] * 1.0e-7;
	m_pdfTLE[ 5 ] = m_dfInc0 = (double) stIOE.pnElement1to6[ 2 ] * 1.0e-5 * g_dfDEG2RAD;
	m_pdfTLE[ 6 ] = m_dfRAAN0 = (double) stIOE.pnElement1to6[ 3 ] * 1.0e-5 * g_dfDEG2RAD; 
	m_pdfTLE[ 7 ] = m_dfPerigee0 = (double) stIOE.pnElement1to6[ 4 ] * 1.0e-5 * g_dfDEG2RAD;
	m_pdfTLE[ 8 ] = m_dfMA0 = (double) stIOE.pnElement1to6[ 5 ] * 1.0e-5 * g_dfDEG2RAD;

	for( int i = 5; i <= 8; i++ ) m_pdfTLE[ i ] *= g_dfRAD2DEG;

	m_dfDaySince1950 = m_dfRefJD - m_dfJD1950;
	m_dfMM0 /= m_dfXPDOTP;	// sat mean motion in radian/minute

	if( !Initialise() ) return FALSE;

	m_bInit = TRUE;

	return TRUE;
}

/******************************************************************************
  Set TLE

 pdfTLE[ 1 ]: year
 pdfTLE[ 2 ]: day of year
 pdfTLE[ 3 ]: bstar
 pdfTLE[ 4 ]: Eccentricity
 pdfTLE[ 5 ]: inclination, in deg
 pdfTLE[ 6 ]: RAAN in deg
 pdfTLE[ 7 ]: argument of perigee in deg
 pdfTLE[ 8 ]: mean anomaly in deg
 pdfTLE[ 9 ]: mean motion in orbits per day
******************************************************************************/

BOOL cTLE2PosVel::SetOrbitalElements( double *pdfTLE )
{
	m_bInit = FALSE;

	for( int i = 0; i < 10; i++ ) m_pdfTLE[ i ] = pdfTLE[ i ];

	int year = (int) ( pdfTLE[ 1 ] + 0.000001 );

	if( year > 50 ) year += 1900;
	else year += 2000;

	g_DateTimeZ.DateTime2JD( year, 1, 0, 0, 0, 0.0, m_dfRefJD );
	m_dfRefJD += pdfTLE[ 2 ];

	m_dfBSTAR = pdfTLE[ 3 ];

	m_dfEcc0 = pdfTLE[ 4 ];
	m_dfInc0 = pdfTLE[ 5 ] * g_dfDEG2RAD;
	m_dfRAAN0 = pdfTLE[ 6 ] * g_dfDEG2RAD; 
	m_dfPerigee0 = pdfTLE[ 7 ] * g_dfDEG2RAD;
	m_dfMA0 = pdfTLE[ 8 ] * g_dfDEG2RAD;
	m_dfMM0 = pdfTLE[ 9 ];

	m_dfDaySince1950 = m_dfRefJD - m_dfJD1950;
	m_dfMM0 /= m_dfXPDOTP;	// sat mean motion in radian/minute

	if( !Initialise() ) return FALSE;

	m_bInit = TRUE;

	return TRUE;
}

/******************************************************************************
  Get TLE

 pdfTLE[ 1 ]: year
 pdfTLE[ 2 ]: day of year
 pdfTLE[ 3 ]: bstar
 pdfTLE[ 4 ]: Eccentricity
 pdfTLE[ 5 ]: inclination, in deg
 pdfTLE[ 6 ]: RAAN in deg
 pdfTLE[ 7 ]: argument of perigee in deg
 pdfTLE[ 8 ]: mean anomaly in deg
 pdfTLE[ 9 ]: mean motion in orbits per day
******************************************************************************/

BOOL cTLE2PosVel::GetOrbitalElements( double *pdfTLE )
{
	for( int i = 0; i < 10; i++ ) pdfTLE[ i ] = m_pdfTLE[ i ];

	return TRUE;
}

/******************************************************************************
 
 Get TLE's ref epoch in JD

******************************************************************************/

void cTLE2PosVel::GetOrbitalElementsRefJD( double &dfRefJD )
{
	dfRefJD = m_dfRefJD;
}

/******************************************************************************************

  Compute the satellite inertial positiona nd velocity at the given time

  dfJD: given time in JD
  pdfPos: computed satellite inertial position, in meter
  pdfVel: computed satellite inertial velocity, in meter/s

*******************************************************************************************/

BOOL cTLE2PosVel::ComputeInertialPosVel( double dfJD, double *pdfPos, double *pdfVel )
{
	if( !m_bInit ) return FALSE;

	if( !SGP4( dfJD, pdfPos, pdfVel ) ) return FALSE;

	return TRUE;
}

/******************************************************************************************

  Compute the satellite inertial positiona nd velocity at the reference epoch of the TLE

  pdfPos: computed satellite inertial position, in meter
  pdfVel: computed satellite inertial velocity, in meter/s

*******************************************************************************************/

BOOL cTLE2PosVel::ComputeInertialPosVel( double &dfIntJDUTCGivenEpoch, double &dfFraJDUTCGivenEpoch, 
			                             double *pdfPos, double *pdfVel )
{
	if( !m_bInit ) return FALSE;

	dfIntJDUTCGivenEpoch = m_dfRefJD;
	dfFraJDUTCGivenEpoch = 0.0;

	g_DateTimeZ.ReConstruct( dfIntJDUTCGivenEpoch, dfFraJDUTCGivenEpoch );

	if( !SGP4( m_dfRefJD, pdfPos, pdfVel ) ) return FALSE;

	return TRUE;
}


/******************************************************************************************

  Compute the satellite EFEC positiona nd velocity at the given time

  dfJD: given time in JD
  pdfPos: computed satellite EFEC position, in meter
  pdfVel: computed satellite EFEC velocity, in meter/s

*******************************************************************************************/

BOOL cTLE2PosVel::ComputeECEFPosVel( double dfJD, double *pdfPos, double *pdfVel )
{
	if( !m_bInit ) return FALSE;

	if( !ComputeInertialPosVel( dfJD, pdfPos, pdfVel ) ) return FALSE;

	FromInertialToEFEC( dfJD, pdfPos, pdfVel );	

	return TRUE;
}


/******************************************************************************************

  Initialisation of SGP4 algorithm, including the initialisation for the deep space effect
  if required.

  FORTRAN Subroutine Initl and first part of SGP4

*******************************************************************************************/

BOOL cTLE2PosVel::Initialise()
{
                                                                            
	// CALCULATE AUXILLARY EPOCH QUANTITIES 
	if( m_dfEcc0 > 0.999999 ) 
	{
		m_nError = 1;
		return FALSE;
	}

	m_dfEcc0SQ = m_dfEcc0 * m_dfEcc0;                                                            
    m_dfOMEcc0SQ = 1.0 - m_dfEcc0SQ;
    m_dfRTEcc0SQ = sqrt( m_dfOMEcc0SQ );
    m_dfCosI0 = cos( m_dfInc0 );
    m_dfCosI02 = m_dfCosI0 * m_dfCosI0;                                                    
                                                                            
	// UN-KOZAI THE MEAN MOTION, see page 10, SGP4 Algorithm
    double dfa1 = pow( m_dfEarthGM / m_dfMM0, m_df2O3 );        // sat orbit semi-major
    double dfD1 = .75 * m_dfJ2 * ( 3.0 * m_dfCosI02 - 1.0 ) / ( m_dfRTEcc0SQ * m_dfOMEcc0SQ );
    double dfDelta1 = dfD1 / ( dfa1 * dfa1 );
    double dfa0 = dfa1 * ( 1.0 - dfDelta1 * dfDelta1 - 
		                   dfDelta1 * ( 1.0 / 3.0 + 134.0 * dfDelta1 * dfDelta1 / 81.0 ) 
						  );
    double dfDelta0 = dfD1 / ( dfa0 * dfa0 );
    
	m_dfMM0 /= ( 1.0 + dfDelta0 );	// mean mean-motion at ref epoch
//	m_dfSM0 = dfa1 / ( 1.0 - dfDelta0 );
    m_dfSM0 = pow( m_dfEarthGM / m_dfMM0, m_df2O3 ); // mean semi-major at ref epoch in ER
    m_dfSinI0 = sin( m_dfInc0 );
    m_dfP0 = m_dfSM0 * m_dfOMEcc0SQ; // p parameter 
    m_dfCon42 = 1.0 - 5.0 * m_dfCosI02;
    m_dfCon41 = -m_dfCon42 - m_dfCosI02 - m_dfCosI02;
    m_dfSMInv = 1.0 / m_dfSM0;
    m_dfEccInv = 1.0/ m_dfEcc0;
    m_dfP0SQ = m_dfP0 * m_dfP0;
    m_dfQ0 = m_dfSM0 * ( 1.0 - m_dfEcc0 );	//dis from earth centre to sat at perigee
	m_dfQ1 = m_dfSM0 * ( 1.0 + m_dfEcc0 );	//dis from earth centre to sat at apogee
     
	if( m_dfQ0 < 1.0 ) // perigee is within the Earth, impossible
	{
		m_nError = 1;
		return FALSE;		
	}                                                                                                             

	// CALCULATE NUMBER OF INTEGER DAYS SINCE 0 JAN 1970            
    m_dfGSTAtRefEpoch = cGreenwichST::ComputeGST( m_dfRefJD );

	// first part of SGP4
	m_bOrbitLowHigh = FALSE;
    if( m_dfQ0 < ( 220.0 / m_dfEarthRadiusKM + 1.0 ) ) m_bOrbitLowHigh = TRUE;
	                          
    double S4 = m_dfMinSatGeoDis;
    double QZMS24 = m_dfQZMS2T;
	double PERIGE = ( m_dfQ0 - 1.0 ) * m_dfEarthRadiusKM;
      
	if( PERIGE >= 156.0 ) goto FOURTY;                                       
    S4 = PERIGE - 78.0;
      
	if( PERIGE > 98.0 ) goto THIRTY;                                        
    S4 = 20.0;
   
THIRTY:
	QZMS24 = pow( ( 120.0 - S4 ) / m_dfEarthRadiusKM, 4.0 );
      
	S4 = S4 / m_dfEarthRadiusKM + 1.0;
   
FOURTY:
	// se page 11 of SGP4 Algorithm
	double PINVSQ = 1.0 / m_dfP0SQ;
    double dfTsi = 1.0 / ( m_dfSM0 - S4 );	// Eq (11) of SGP4 Algorithm
    m_dfEta = m_dfSM0 * m_dfEcc0 * dfTsi;
    double ETASQ = m_dfEta * m_dfEta;	
    double EETA = m_dfEcc0 * m_dfEta;
    double PSISQ = fabs( 1.0 - ETASQ );
    double COEF = QZMS24 * pow( dfTsi, 4.0 );
    double COEF1 = COEF / pow( PSISQ, 3.5 );
    m_dfCC2 = COEF1 * m_dfMM0 * 
		      ( m_dfSM0 * ( 1.0 + 1.50 * ETASQ + EETA * ( 4.0 + ETASQ ) ) +
                .3750 * m_dfJ2 * dfTsi / PSISQ * m_dfCon41 * ( 8.0 + 3.0 * ETASQ * ( 8.0 + ETASQ ) ) 
              );						// Eq (14) of SGP4 Algorithm
    m_dfCC1 = m_dfBSTAR * m_dfCC2;		// Eq (15) of SGP4 Algorithm
    m_dfCC3 = 0.0;	
    if( m_dfEcc0 > 1.e-4 ) 
		m_dfCC3 = -2.0 * COEF * dfTsi * m_dfJ3OJ2 * m_dfMM0 * m_dfSinI0 / m_dfEcc0;  // Eq (16) of SGP4 Algorithm
    m_dfX1MTH2 = 1.0 - m_dfCosI02;
    m_dfCC4 = 2.0 * m_dfMM0 * COEF1 * m_dfSM0 * m_dfOMEcc0SQ * 
             ( m_dfEta * ( 2.0 + .50 * ETASQ ) + m_dfEcc0 * ( .50 + 2.0 * ETASQ ) - 
               m_dfJ2 * dfTsi / ( m_dfSM0 * PSISQ ) * 
			   ( -3.0 * m_dfCon41 * ( 1.0 - 2.0 * EETA + ETASQ * ( 1.5 - .5 * EETA ) ) +
                 .75 * m_dfX1MTH2 * ( 2.0 * ETASQ - EETA * ( 1.0 + ETASQ ) ) * cos( 2.0 * m_dfPerigee0 ) 
			   )
             );														// Eq (17) of SGP4 Algorithm
	m_dfCC5 = 2.0 * COEF1 * m_dfSM0 * m_dfOMEcc0SQ * 
		      ( 1.0 + 2.75 * ( ETASQ + EETA ) + EETA * ETASQ );		// Eq (18) of SGP4 Algorithm
    double COSIO4 = m_dfCosI02 * m_dfCosI02;
    double TEMP1 = 1.5 * m_dfJ2 * PINVSQ * m_dfMM0;                                            
    double TEMP2 = .5 * TEMP1 * m_dfJ2 * PINVSQ;
    double TEMP3 = -.46875* m_dfJ4 * PINVSQ * PINVSQ * m_dfMM0;                                 
      
	m_dfMMDot = m_dfMM0 + .5 * TEMP1 * m_dfRTEcc0SQ * m_dfCon41 + 
                .0625 * TEMP2 * m_dfRTEcc0SQ * 
				( 13.0 - 78.0 * m_dfCosI02 + 137.0 * COSIO4 );		// Eq (22) of SGP4 Algorithm
    m_dfPerigeeDot = -.5 * TEMP1 * m_dfCon42 + 
		        .0625 * TEMP2 * ( 7.0 - 114.0 * m_dfCosI02 + 395.0 * COSIO4 ) + 
                TEMP3 * ( 3.0 - 36.0 * m_dfCosI02 + 49.0 * COSIO4 );	// Eq (23) of SGP4 Algorithm
    double XHDOT1 = -TEMP1 * m_dfCosI0;                                                   
    m_dfRAANDot = XHDOT1 + ( .50 * TEMP2 * ( 4.0 - 19.0 * m_dfCosI02 ) +
                            2.0 * TEMP3 * ( 3.0 - 7.0 * m_dfCosI02 ) 
	                      ) * m_dfCosI0;								// Eq (24) of SGP4 Algorithm
                                             
	m_dfXPIDOT = m_dfPerigeeDot + m_dfRAANDot;

    m_dfPerigeeDotDrag = m_dfBSTAR * m_dfCC3 * cos( m_dfPerigee0 );		// Eq (25) of SGP4 Algorithm
    m_dfMMDotDrag = 0.0;
    if( m_dfEcc0 > 1.e-4) m_dfMMDotDrag = -m_df2O3 * COEF * m_dfBSTAR / EETA;	// Eq (26) of SGP4 Algorithm
    m_dfRAANDot2Drag = 3.5 * m_dfOMEcc0SQ * XHDOT1 * m_dfCC1;					// Eq (29) of SGP4 Algorithm
    
	m_dfT2Coe = 1.5 * m_dfCC1;							// Coe of T^2, Eq (32) of SGP4 Algorithm
    m_dfXLCOF = -.25 * m_dfJ3OJ2 * m_dfSinI0 * ( 3.0 + 5.0 * m_dfCosI0 ) / 
		        ( 1.0 + m_dfCosI0 );					// Eq (36) of SGP4 Algorithm
    m_dfAYCOF = -.5 * m_dfJ3OJ2 * m_dfSinI0;			// Eq (37) of SGP4 Algorithm
    m_dfDeltaMA0 = pow( 1.0 + m_dfEta * cos( m_dfMA0 ), 3.0 );	// last term in the bracket of Eq (26) of SGP4 Algorithm
    m_dfSinMAO = sin( m_dfMA0 );
    m_dfX7THM1 = 7.0 * m_dfCosI02 - 1.0;

	m_bDeepSpace = FALSE;				// initially, assuming deep space effect ignored    
	if( g_dfTWOPI / m_dfMM0 >= 225.0 )	// satellite higher than 5876 km ( m_bDeepSpace == TRUE )     
	{
		m_bDeepSpace = TRUE;
		m_bOrbitLowHigh = TRUE;
	
		if( !InitialiseDeepSpace() ) return FALSE;
      
		if( m_dfEcc0 < 0.0 || m_dfEcc0 > 1.0 ) 
		{
			m_nError = 1;
			return FALSE;		
		}                
	}   
	
	if( m_bOrbitLowHigh ) return TRUE;					// orbit either high or low
      
	// perigee height > 220.0km, non-linear part of the drag term coefficients
	double CC1SQ = m_dfCC1 * m_dfCC1;
	m_dfD2 = 4.0 * m_dfSM0 * dfTsi * CC1SQ;				// Eq (19) of SGP4 Algorithm
	double dfTemp = m_dfD2 * dfTsi * m_dfCC1 / 3.0;	
	m_dfD3 = ( 17.0 * m_dfSM0 + S4 ) * dfTemp;			// Eq (20) of SGP4 Algorithm                                         
    m_dfD4 = .5 * dfTemp * m_dfSM0 * dfTsi * 
		     ( 221.0 * m_dfSM0 + 31.0 * S4 ) * m_dfCC1;	// Eq (21) of SGP4 Algorithm
	
	m_dfT3Coe = m_dfD2 + 2.0 * CC1SQ;					// Coe of T^3, Eq (32) of SGP4 Algorithm
	m_dfT4Coe = .25 * ( 3.0 * m_dfD3 + 
		                m_dfCC1 * (12.0 * m_dfD2 + 10.0 * CC1SQ ) );	// Coe of T^4, Eq (32) of SGP4 Algorithm
	m_dfT5Coe = .2 * ( 3.0 * m_dfD4 + 12.0 * m_dfCC1 * m_dfD3 + 6.0 * m_dfD2 * m_dfD2 + 
		               15.0 * CC1SQ * (2.0 * m_dfD2 + CC1SQ ) 
					  );								// Coe of T^5, Eq (32) of SGP4 Algorithm
                                                                        
	return TRUE;
}

/******************************************************************************************

  Implementation of the SGP4 algorithm

*******************************************************************************************/

BOOL cTLE2PosVel::SGP4( double dfJD, double *pdfPos, double *pdfVel )
{
	double dfT = ( dfJD - m_dfRefJD ) * 1440.0;	// in minute
    double dfT2 = dfT * dfT;                                                                
                                                  
	// secular gravity terms
	double dfMAS = m_dfMA0 + m_dfMMDot * dfT;					// Eq (22) of SGP4 Algorithm
    double dfPerifeeS = m_dfPerigee0 + m_dfPerigeeDot * dfT;	// Eq (23) of SGP4 Algorithm
    double dfRAANS = m_dfRAAN0 + m_dfRAANDot * dfT;				// Eq (24) of SGP4 Algorithm
	
    m_dfPerigeeM = dfPerifeeS;
	m_dfMAM = dfMAS;                                        
    m_dfRAANM = dfRAANS + m_dfRAANDot2Drag * dfT2;				// Eq (29) of SGP4 Algorithm
	                                    
	// secular drag terms
    double dfTempA = 1.0 - m_dfCC1 * dfT;			// first two terms in the bracket of Eq (31) of SGP4 Algorithm
	double dfTempE = m_dfBSTAR * m_dfCC4 * dfT;		// second term of Eq (30) of SGP4 Algorithm
    double dfTempL = m_dfT2Coe * dfT2;				// first term in the Braket of Eq (32) of SGP4 Algorithm

	double dfTemp;

	if( !m_bOrbitLowHigh )
	{
		// ( perigee height > 220km )                                       
		double dfT3 = dfT2 * dfT;                                                               
		double dfT4 = dfT3 * dfT;   
		double dfDeltaPerigee = m_dfPerigeeDotDrag * dfT;	// Eq (25) of SGP4 Algorithm
		double dfDeltaMA = m_dfMMDotDrag * 
			               ( pow( 1.0 + m_dfEta * cos( dfMAS ),  3.0 ) - m_dfDeltaMA0 );
		dfTemp = dfDeltaPerigee + dfDeltaMA;
		m_dfMAM = dfMAS + dfTemp;				// Eq (27) of SGP4 Algorithm
		m_dfPerigeeM = dfPerifeeS - dfTemp;		// Eq (28) of SGP4 Algorithm
		dfTempA = dfTempA - m_dfD2 * dfT2 - m_dfD3 * dfT3 - m_dfD4 * dfT4;				// Ref. Eq (31) of SGP4 Algorithm
		dfTempE = dfTempE + m_dfBSTAR * m_dfCC5 * ( sin( m_dfMAM ) - m_dfSinMAO );		// Ref. Eq (30) of SGP4 Algorithm
		dfTempL = dfTempL + m_dfT3Coe * dfT3 + dfT4 *( m_dfT4Coe + dfT * m_dfT5Coe );	// Ref. Eq (31) of SGP4 Algorithm                             
	}

	// Secular Deep Space Effect
	m_dfEccM = m_dfEcc0;
	m_dfIncM = m_dfInc0;
	m_dfMMM = m_dfMM0;
		                                                      
    if( m_bDeepSpace )
	{
		m_dfTC = dfT;
		if( !DeepSpaceSecularEffect() ) return FALSE;            
	}                                                               

    m_dfSMM = pow( m_dfEarthGM / m_dfMMM, m_df2O3 ) * pow( dfTempA, 2.0 );	// Eq (31) of SGP4 Algorithm
	m_dfMMM = m_dfEarthGM / pow( m_dfSMM, 1.5 );
	m_dfEccM = m_dfEccM - dfTempE;											// Eq (30) of SGP4 Algorithm
	if( m_dfEccM >= 1.0 || m_dfEccM < -1.0e-3 ) 	
	{
		m_nError = 1;
		return FALSE;
	}
                                            
	if ( m_dfEccM < 0.0 ) m_dfEccM = 1.e-6;	
    m_dfMAM += m_dfMM0 * dfTempL;	
	double dfXLM = m_dfMAM + m_dfPerigeeM + m_dfRAANM;						// Eq (32) of SGP4 Algorithm
	double dfEMSQ = m_dfEccM * m_dfEccM;
    dfTemp = 1.0 - dfEMSQ;                            
    double dfRTEMSQ = sqrt( dfTemp );
	m_dfRAANM = fmod( m_dfRAANM, g_dfTWOPI );
	m_dfPerigeeM = fmod( m_dfPerigeeM, g_dfTWOPI );
	dfXLM = fmod( dfXLM, g_dfTWOPI );
	m_dfMAM = fmod( dfXLM - m_dfPerigeeM - m_dfRAANM, g_dfTWOPI );
                                                                           
	// COMPUTE EXTRA MEAN QUANTITIES                                                                          
	double dfSinIM = sin( m_dfIncM );
	double dfCosIM = cos( m_dfIncM );
	double dfSinPerigeeM = sin( m_dfPerigeeM );
	double dfCosPerigeeM = cos( m_dfPerigeeM );
	
	
	// ADD LUNAR-SOLAR PERIODICS                                                                                                                       
	double dfMMP = m_dfMMM;
	double dfEccP = m_dfEccM;
	double dfIncP = m_dfIncM;
	double dfPerigeeP = m_dfPerigeeM;
	double dfRAANP = m_dfRAANM;
	double dfMAP = m_dfMAM;
	double dfSinIP = dfSinIM;
	double dfCosIP = dfCosIM;
	
	if( m_bDeepSpace ) 
	{
		if( !DeepSpacePeriodicEffect( dfEccP, dfIncP, dfRAANP, dfPerigeeP, dfMAP ) ) return FALSE;
		if( dfIncP < 0.0 ) 
		{
			dfIncP = - dfIncP;
			dfRAANP += g_dfPI;
			dfPerigeeP -= g_dfPI;
		}
	}

	if( dfEccP < 0.0 || dfEccP > 1.0 ) 
	{
		m_nError = 1;
		return FALSE;
	}
                                                                           
	// LONG PERIOD PERIODICS                                                                          
	if( m_bDeepSpace ) 
	{
		dfSinIP = sin( dfIncP );
		dfCosIP = cos( dfIncP );                                               
		m_dfAYCOF = -.5 * m_dfJ3OJ2 * dfSinIP;			// Eq (37) of SGP4 Algorithm
		m_dfXLCOF = -.25 * m_dfJ3OJ2 * dfSinIP * ( 3.0 + 5.0 * dfCosIP ) / 
			        ( 1.0 + dfCosIP );					// Eq (36) of SGP4 Algorithm
	}

	double dfAXNL = dfEccP * cos( dfPerigeeP );									// Eq (35) of SGP4 Algorithm
	dfTemp = 1.0 / ( m_dfSMM * ( 1.0 - dfEccP * dfEccP ) );
	double dfAYNL = dfEccP * sin( dfPerigeeP ) + dfTemp * m_dfAYCOF;			// Eq (39) of SGP4 Algorithm
	double dfXL = dfMAP + dfPerigeeP + dfRAANP + dfTemp * m_dfXLCOF * dfAXNL;	// Eq (38) of SGP4 Algorithm
	
	// SOLVE KEPLER'S EQUATION, Eqs (40)-(43) of SGP4 Algorthm	                                                                           
    double dfU = fmod( dfXL - dfRAANP, g_dfTWOPI );
	double dfE01 = dfU;
	double dfSinE01, dfCosE01;
	double dfTEM5;
	int nIters = 0;
	do
	{
		dfSinE01 = sin( dfE01 );
		dfCosE01 = cos( dfE01 );
		dfTEM5 = 1.0 - dfCosE01 * dfAXNL - dfSinE01 * dfAYNL;
		dfTEM5=( dfU - dfAYNL * dfCosE01 + dfAXNL * dfSinE01 - dfE01 ) / dfTEM5;
		dfE01 += dfTEM5;
		nIters++;
	} while( fabs( dfTEM5 ) > 1.0e-12 && nIters < 10 );
                                                                            
	// SHORT PERIOD PRELIMINARY QUANTITIES                                                                                                             
	double dfeCosE = dfAXNL * dfCosE01 + dfAYNL * dfSinE01;		// Eq (44) of SGP4 Algorithm
	double dfeSinE = dfAXNL * dfSinE01 - dfAYNL * dfCosE01;		// Eq (45) of SGP4 Algorithm
	double dfEL2 = dfAXNL * dfAXNL + dfAYNL * dfAYNL;			// Eq (46) of SGP4 Algorithm
	double dfPL = m_dfSMM * ( 1.0 - dfEL2 );					// Eq (47) of SGP4 Algorithm
	
	if( dfPL < 0.0 ) 
	{
		m_nError = 1;
		return FALSE;
	}
		
	double dfRL = m_dfSMM * ( 1.0 - dfeCosE );				// Eq (48) of SGP4 Algorithm
    double dfRDotL = sqrt( m_dfSMM ) * dfeSinE / dfRL;		// Eq (49) of SGP4 Algorithm
	double dfRVDotL = sqrt( dfPL ) / dfRL;					// Eq (50) of SGP4 Algorithm
	double dfBETAL = sqrt( 1.0 - dfEL2 );					// sqrt root part of Eqs (51) and (52) of SGP4 Algorithm
	dfTemp = dfeSinE / ( 1.0 + dfBETAL );					// Eqs (51) and (52) of SGP4 Algorithm
	double dfSinU = m_dfSMM / dfRL * ( dfSinE01 - dfAYNL - dfAXNL * dfTemp );	// Eq (52) of SGP4 Algorithm                                    
    double dfCosU = m_dfSMM / dfRL * ( dfCosE01 - dfAXNL + dfAYNL * dfTemp );	// Eq (51) of SGP4 Algorithm
	double dfSU = atan2( dfSinU, dfCosU );			// Eq (53) of SGP4 Algorithm
	double dfSin2U = ( dfCosU + dfCosU ) * dfSinU;
	double dfCos2U = 1.0 - 2.0 * dfSinU * dfSinU;
	dfTemp = 1.0 / dfPL;                                                          
    double dfTEMP1 = .5 * m_dfJ2 * dfTemp;			// Eq (54) of SGP4 Algorithm
	double dfTEMP2 = dfTEMP1 * dfTemp;				// Eq (55) of SGP4 Algorithm
                                                                           
	// UPDATE FOR SHORT PERIOD PERIODICS                                                                          
	if( m_bDeepSpace ) 
	{
		double m_dfCosI2 = dfCosIP * dfCosIP;
		m_dfCon41 = 3.0 * m_dfCosI2 - 1.0;
		m_dfX1MTH2 = 1.0 - m_dfCosI2;
		m_dfX7THM1 = 7.0 * m_dfCosI2 - 1.0;
	}

	double dfRT = dfRL * ( 1.0 - 1.5 * dfTEMP2 * dfBETAL * m_dfCon41 ) +		// Eq (60) of SGP4 Algorithm
		          0.5 * dfTEMP1 * m_dfX1MTH2 * dfCos2U;							// Eq (54) of SGP4 Algorithm
	dfSU = dfSU - .25 * dfTEMP2 * m_dfX7THM1 * dfSin2U;							// Eqs (55) and (61)of SGP4 Algorithm
	double dfRAANT = dfRAANP + 1.5 * dfTEMP2 * dfCosIP * dfSin2U;				// Eqs (56) and (62)of SGP4 Algorithm
	double dfIncT = dfIncP + 1.5 * dfTEMP2 * dfCosIP * dfSinIP * dfCos2U;		// Eqs (57) and (63)of SGP4 Algorithm
	double dfRDot = dfRDotL - 
		            m_dfMMM * dfTEMP1 * m_dfX1MTH2 * dfSin2U / m_dfEarthGM;		// Eqs (58) and (64)of SGP4 Algorithm
	double dfRVDot = dfRVDotL + 
		             m_dfMMM * dfTEMP1 * 
					 ( m_dfX1MTH2 * dfCos2U + 1.5 * m_dfCon41 ) / m_dfEarthGM;	// Eqs (59) and (65)of SGP4 Algorithm
                                                                          
	// ORIENTATION VECTORS	                                                             
	double dfSinSU = sin( dfSU );
	double dfCosSU = cos( dfSU );                                                         
    double dfSinNode = sin( dfRAANT );
	double dfCosNode = cos( dfRAANT );
	double dfSinI = sin( dfIncT );
	double dfCosI = cos( dfIncT );     
	
    double dfMx = - dfSinNode * dfCosI;		// Eq (68) of SGP4 Algorithm
	double dfMy =   dfCosNode * dfCosI;		// Eq (68) of SGP4 Algorithm
	
	// Eq (66) of SGP4 Algorithm
	double dfUX = dfMx *  dfSinSU + dfCosNode * dfCosSU;	
	double dfUY = dfMy *  dfSinSU + dfSinNode * dfCosSU;
	double dfUZ = dfSinI*  dfSinSU;
                                                                           
	// POSITION AND VELOCITY

	// Eq (70) of SGP4 Algorithm
    pdfPos[ 0 ] = dfRT * dfUX * m_dfEarthRadius;	// in meter
	pdfPos[ 1 ] = dfRT * dfUY * m_dfEarthRadius; 
	pdfPos[ 2 ] = dfRT * dfUZ * m_dfEarthRadius; 

	if( m_bPosOnly ) return TRUE;

	// Eq (67) of SGP4 Algorithm
    double dfVX = dfMx *  dfCosSU - dfCosNode * dfSinSU;
	double dfVY = dfMy *  dfCosSU - dfSinNode * dfSinSU;
	double dfVZ = dfSinI * dfCosSU;                                                         

	// Eq (71) of SGP4 Algorithm
	pdfVel[ 0 ] = ( dfRDot * dfUX + dfRVDot * dfVX ) * m_dfVelocityChange;	// in m/s
	pdfVel[ 1 ] = ( dfRDot * dfUY + dfRVDot * dfVY ) * m_dfVelocityChange;
	pdfVel[ 2 ] = ( dfRDot * dfUZ + dfRVDot * dfVZ ) * m_dfVelocityChange;
      
	return TRUE;
}

/******************************************************************************************

  Initialisation of the deep space effect

*******************************************************************************************/

BOOL cTLE2PosVel::InitialiseDeepSpace()
{
	// first part of Fortran Subroutine DPPER
	
	// DEEP SPACE PERIODICS INITIALIZATION
	m_dfMMM = m_dfMM0;
	m_dfEccM = m_dfEcc0;
	m_dfIncM = m_dfInc0;
	m_dfSinRAANM = sin( m_dfRAAN0 );
	m_dfCosRAANM = cos( m_dfRAAN0 );
	m_dfSinPerigeeM = sin( m_dfPerigee0 );
	m_dfCosPerigeeM = cos( m_dfPerigee0 );
	m_dfSinIncM = sin( m_dfInc0 );
	m_dfCosIncM = cos( m_dfInc0 );
                                                                          
	// INITIALIZE LUNAR SOLAR TERMS
	if( !DSCOM() ) return FALSE;

	m_dfZMOL = fmod( 4.7199672 + .22997150 * m_dfDSDay - m_dfDSGAM, g_dfTWOPI );
	m_dfZMOS = fmod( 6.2565837 + .017201977 * m_dfDSDay, g_dfTWOPI );

	// DO SOLAR TERMS	                                                                           
    m_dfDSSE2  =  2.0 * m_dfDSSS1 * m_dfDSSS6;
	m_dfDSSE3  =  2.0 * m_dfDSSS1 * m_dfDSSS7;
	m_dfDSSI2  =  2.0 * m_dfDSSS2 * m_dfDSSZ12;
	m_dfDSSI3  =  2.0 * m_dfDSSS2 * ( m_dfDSSZ13 - m_dfDSSZ11 );
    m_dfDSSL2  = -2.0 * m_dfDSSS3 * m_dfDSSZ2;
	m_dfDSSL3  = -2.0 * m_dfDSSS3 * ( m_dfDSSZ3 - m_dfDSSZ1 );
	m_dfDSSL4  = -2.0 * m_dfDSSS3 * ( -21.0 - 9.0 * m_dfEccMSQ ) * m_dfZES;
	m_dfDSSGH2 =  2.0 * m_dfDSSS4 * m_dfDSSZ32;
	m_dfDSSGH3 =  2.0 * m_dfDSSS4 * ( m_dfDSSZ33 - m_dfDSSZ31 );
	m_dfDSSGH4 =-18.0 * m_dfDSSS4 * m_dfZES;
	m_dfDSSH2  = -2.0 * m_dfDSSS2 * m_dfDSSZ22;
	m_dfDSSH3  = -2.0 * m_dfDSSS2 * ( m_dfDSSZ23 - m_dfDSSZ21 );

	// DO LUNAR TERMS	                                                                           
    m_dfDSEE2  =  2.0 * m_dfDSS1 * m_dfDSS6;
	m_dfDSE3   =  2.0 * m_dfDSS1 * m_dfDSS7;
	m_dfDSXI2  =  2.0 * m_dfDSS2 * m_dfDSZ12;
	m_dfDSXI3  =  2.0 * m_dfDSS2 *( m_dfDSZ13 - m_dfDSZ11 );
	m_dfDSXL2  = -2.0 * m_dfDSS3 * m_dfDSZ2;
	m_dfDSXL3  = -2.0 * m_dfDSS3 * ( m_dfDSZ3 - m_dfDSZ1 );
	m_dfDSXL4  = -2.0 * m_dfDSS3 * ( -21.0 - 9.0 * m_dfEccMSQ ) * m_dfZEL;
	m_dfDSXGH2 =  2.0 * m_dfDSS4 * m_dfDSZ32;
	m_dfDSXGH3 =  2.0 * m_dfDSS4 * ( m_dfDSZ33 - m_dfDSZ31 );
	m_dfDSXGH4 =-18.0 * m_dfDSS4 * m_dfZEL;
	m_dfDSXH2  = -2.0 * m_dfDSS2 * m_dfDSZ22;
	m_dfDSXH3  = -2.0 * m_dfDSS2 * ( m_dfDSZ23 - m_dfDSZ21 );

	// end of first part of Fortran Subroutine DPPER

/*************************************************************************************************/

	// first part of FORTRAN subroutine DSPACE

	m_nIREZ = 0;
	if( m_dfMMM < .0052359877 && m_dfMMM > .003490658 ) m_nIREZ = 1;
	if( m_dfMMM >= 8.26e-3 && m_dfMMM <= 9.24e-3 && m_dfEccM >= 0.5 ) m_nIREZ = 2;

	// DO SOLAR TERMS
	double SES = m_dfDSSS1 * m_dfZNS * m_dfDSSS5;
	double SIS = m_dfDSSS2 * m_dfZNS * ( m_dfDSSZ11 + m_dfDSSZ13 );
	double SLS = -m_dfZNS * m_dfDSSS3 * ( m_dfDSSZ1 + m_dfDSSZ3 - 14.0 - 6.0 * m_dfEccMSQ );
    double SGHS = m_dfDSSS4 * m_dfZNS * ( m_dfDSSZ31 + m_dfDSSZ33 - 6.0 );
	double SHS = -m_dfZNS * m_dfDSSS2 * ( m_dfDSSZ21 + m_dfDSSZ23 );
	if( m_dfIncM < 5.2359877e-2 ) SHS = 0.0;
	if( m_dfSinIncM != 0.0 ) SHS = SHS / m_dfSinIncM;
    double SGS = SGHS - m_dfCosIncM * SHS;

	// DO LUNAR TERMS
	m_dfDEDT = SES + m_dfDSS1 * m_dfZNL * m_dfDSS5;
	m_dfDIDT = SIS + m_dfDSS2 * m_dfZNL * ( m_dfDSZ11 + m_dfDSZ13 );
	m_dfDMDT = SLS - m_dfZNL * m_dfDSS3 * ( m_dfDSZ1 + m_dfDSZ3 - 14.0 - 6.0 * m_dfEccMSQ );
	double SGHL = m_dfDSS4 * m_dfZNL * ( m_dfDSZ31 + m_dfDSZ33 - 6.0 );
	double SHL = -m_dfZNL * m_dfDSS2 * ( m_dfDSZ21 + m_dfDSZ23 );
	if( m_dfIncM < 5.2359877e-2 ) SHL = 0.0;
    m_dfDOMDT = SGS + SGHL;                                        
    m_dfDNODT = SHS;
	if( m_dfSinIncM != 0.0 )
	{
		m_dfDOMDT = m_dfDOMDT - m_dfCosIncM  / m_dfSinIncM * SHL;
		m_dfDNODT = m_dfDNODT + SHL / m_dfSinIncM;
	}

	// end of first part of FORTRAN subroutine DSPACE

/*************************************************************************************************/

	// first part of FORTRAN subroutine SREZ

	double Q22 = 1.7891679e-6, Q31 = 2.1460748e-6, Q33 = 2.2123015e-7;
	double ROOT22 = 1.7891679e-6, ROOT32 = 3.7393792e-7;
	double ROOT44 = 7.3636953e-9, ROOT52 = 1.1428639e-7;
	double ROOT54 = 2.1765803e-9;
	
    double AONV = pow( m_dfMMM / m_dfEarthGM, m_df2O3 );
	
	if( m_nIREZ == 2 ) // GEOPOTENTIAL RESONANCE FOR 12 HOUR ORBITS 
	{
		double COSISQ = m_dfCosIncM * m_dfCosIncM;
		double EOC = m_dfEcc0 * m_dfEcc0SQ;
		double G201 = -.306 - ( m_dfEcc0 - .64 ) * .440;
		
		double G211, G310, G322, G410, G422, G520;

		if( m_dfEcc0 < .65 )
		{
			G211 =  3.616 - 13.247 * m_dfEcc0 + 16.290 * m_dfEcc0SQ;
			G310 = -19.302 + 117.390 * m_dfEcc0 - 228.419 * m_dfEcc0SQ + 156.591 * EOC;
			G322 = -18.9068 + 109.7927 * m_dfEcc0 - 214.6334 * m_dfEcc0SQ + 146.5816 * EOC;
			G410 = -41.122 + 242.694 * m_dfEcc0 - 471.094 * m_dfEcc0SQ + 313.953 * EOC;
			G422 = -146.407 + 841.880 * m_dfEcc0 - 1629.014 * m_dfEcc0SQ + 1083.435 * EOC;
			G520 = -532.114 + 3017.977 * m_dfEcc0 - 5740.032 * m_dfEcc0SQ + 3708.276 * EOC;
		}
		else
		{
			G211 = -72.099 + 331.819 * m_dfEcc0 - 508.738 * m_dfEcc0SQ + 266.724 * EOC;
			G310 = -346.844 + 1582.851 * m_dfEcc0 - 2415.925 * m_dfEcc0SQ + 1246.113 * EOC;
			G322 = -342.585 + 1554.908 * m_dfEcc0 - 2366.899 * m_dfEcc0SQ + 1215.972 * EOC;
			G410 = -1052.797 + 4758.686 * m_dfEcc0 - 7193.992 * m_dfEcc0SQ + 3651.957 * EOC;
			G422 = -3581.69 + 16178.11 * m_dfEcc0 - 24462.77 * m_dfEcc0SQ + 12422.52 * EOC;
			if( m_dfEcc0 < .715 ) 
				G520 = 1464.74 - 4664.75 * m_dfEcc0 + 3763.64 * m_dfEcc0SQ;                  
			else
			   G520 = -5149.66 + 29936.92 * m_dfEcc0 - 54087.36 * m_dfEcc0SQ + 31324.56 * EOC;
		}

		double G533, G521, G532;

		if( m_dfEcc0 < .7 )
		{
			G533 = -919.2277 + 4988.61 * m_dfEcc0 - 9064.77 * m_dfEcc0SQ + 5542.21 * EOC;
			G521 = -822.71072 + 4568.6173 * m_dfEcc0 - 8491.4146 * m_dfEcc0SQ + 5337.524 * EOC;
			G532 = -853.666 + 4690.25 * m_dfEcc0 - 8624.77 * m_dfEcc0SQ + 5341.4 * EOC;
		}
		else
		{
			G533 = -37995.78 + 161616.52 * m_dfEcc0 - 229838.2 * m_dfEcc0SQ + 109377.94 * EOC;
			G521 = -51752.104 + 218913.95 * m_dfEcc0 - 309468.16 * m_dfEcc0SQ + 146349.42 * EOC;
			G532 = -40023.88 + 170470.89 * m_dfEcc0 - 242699.48 * m_dfEcc0SQ + 115605.82 * EOC; 
		}

		double SINI2 = m_dfSinIncM * m_dfSinIncM;
		double F220 = .75 * (1.0 + 2.0 * m_dfCosIncM + COSISQ );
		double F221 = 1.5 * SINI2;
		double F321 = 1.875 * m_dfSinIncM * ( 1.0 - 2.0 * m_dfCosIncM - 3.0 * COSISQ );
		double F322 = -1.875 * m_dfSinIncM * ( 1.0 + 2.0 * m_dfCosIncM - 3.0 * COSISQ );
		double F441 = 35.0 * SINI2 * F220;
		double F442 = 39.375 * SINI2 * SINI2;
		double F522 = 9.84375 * m_dfSinIncM * ( SINI2 * ( 1.0 - 2.0 * m_dfCosIncM - 5.0 * COSISQ ) +
                                               .33333333 * ( -2.0 + 4.0 * m_dfCosIncM + 6.0 * COSISQ )
											  );
		double F523 = m_dfSinIncM * ( 4.92187512 * SINI2 * ( -2.0 - 4.0 * m_dfCosIncM + 10.0 * COSISQ ) +
			                    6.56250012 * ( 1.0 + 2.0 * m_dfCosIncM - 3.0 * COSISQ ) 
							  );
		double F542 = 29.53125 * m_dfSinIncM * 
			          ( 2.0 - 8.0 * m_dfCosIncM + COSISQ * ( -12.0 + 8.0 * m_dfCosIncM + 10.0 * COSISQ ) );
		double F543 = 29.53125 * m_dfSinIncM * 
			          ( -2.0 - 8.0 * m_dfCosIncM + COSISQ * ( 12.0 + 8.0 * m_dfCosIncM - 10.0 * COSISQ ) );

		double XNO2 = m_dfMMM * m_dfMMM;
		double AINV2 = AONV * AONV;

		double TEMP1 = 3.0 * XNO2 * AINV2;
		double TEMP = TEMP1 * ROOT22;                                                   
		m_dfD2201 = TEMP * F220 * G201;
		m_dfD2211 = TEMP * F221 * G211;
		TEMP1 = TEMP1 * AONV;
		TEMP = TEMP1 * ROOT32;
		m_dfD3210 = TEMP * F321 * G310;
		m_dfD3222 = TEMP * F322 * G322;
		TEMP1 = TEMP1 * AONV;
		TEMP = 2.0 * TEMP1 * ROOT44;
		m_dfD4410 = TEMP * F441 * G410;
		m_dfD4422 = TEMP * F442 * G422;
		TEMP1 = TEMP1 * AONV;
		TEMP = TEMP1 * ROOT52;
		m_dfD5220 = TEMP * F522 * G520;
		m_dfD5232 = TEMP * F523 * G532;
		TEMP = 2.0 * TEMP1 * ROOT54;
		m_dfD5421 = TEMP * F542 * G521;
		m_dfD5433 = TEMP * F543 * G533;                                                
		m_dfXLAMO = fmod( m_dfMA0 + m_dfRAAN0 * 2.0 - m_dfGSTAtRefEpoch * 2.0, g_dfTWOPI );
	}
	else	// SYNCHRONOUS RESONANCE TERMS, 24 hour orbits
	{
		double G200 = 1.0 + m_dfEcc0SQ * ( -2.5 + .8125 * m_dfEcc0SQ );
		double G310 = 1.0 + 2.0 * m_dfEcc0SQ;
		double G300 = 1.0 + m_dfEcc0SQ * (-6.0 + 6.60937 * m_dfEcc0SQ );
		double F220 = .75 * ( 1.0 + m_dfCosIncM ) * ( 1.0 + m_dfCosIncM );
		double F311 = .9375 * m_dfSinIncM * m_dfSinIncM * ( 1.0 + 3.0 * m_dfCosIncM ) -
			          .75 * ( 1.0 + m_dfCosIncM );
		double F330 = 1.0 + m_dfCosIncM;
		F330 = 1.875 * F330 * F330 * F330;
		m_dfDEL1 = 3.0 * m_dfMMM * m_dfMMM * AONV * AONV;
		m_dfDEL2 = 2.0 * m_dfDEL1 * F220 * G200 * Q22;
		m_dfDEL3 = 3.0 * m_dfDEL1 * F330 * G300 * Q33 * AONV;
		m_dfDEL1 = m_dfDEL1 * F311 * G310 * Q31 * AONV;
		m_dfFASX2 =.13130908;
		m_dfFASX4 =2.8843198;
		m_dfFASX6 =.37448087;                                                     
		m_dfXLAMO = fmod( m_dfMA0 + m_dfRAAN0 + m_dfPerigee0 - m_dfGSTAtRefEpoch, g_dfTWOPI );
	}

	// FOR SGP4, INITIALIZE THE resonance INTEGRATOR                                                                           
	m_dfXLI = m_dfXLAMO;
	m_dfXNI = m_dfMM0;
    m_dfATIME = 0.0;
	m_dfSTEPP = 720.0;
	m_dfSTEPN = -720.0;
	m_dfSTEP2 = 259200.0;
	
	if( m_nIREZ == 1 ) 
		m_dfXFACT = m_dfMMDot + m_dfXPIDOT - m_dfEarthRotationPerMinute + 
		            m_dfDMDT + m_dfDOMDT + m_dfDNODT - m_dfMM0;
	else m_dfXFACT = m_dfMMDot + m_dfDMDT + 
		             2.0 * ( m_dfRAANDot + m_dfDNODT - m_dfEarthRotationPerMinute ) - m_dfMM0;

	// end of first part of FORTRAN subroutine SREZ

/*************************************************************************************************/
	
	return TRUE;
}

/******************************************************************************************

  Compute the constants of the lunar-solar effects

  FORTRTAN Subroutine: DSCOM

*******************************************************************************************/

BOOL cTLE2PosVel::DSCOM()
{
	double dfC1SS = 2.9864797e-6, dfC1L = 4.7968065e-7, dfZCOSIS = .91744867;
	double dfZSINIS = .39785416, dfZSINGS = -.98088458, dfZCOSGS = .1945905;
	
	// DEEP SPACE INITIALIZATION                                                                          
    m_dfEccMSQ = m_dfEccM * m_dfEccM;                                                            
    double dfBetaSQ = 1.0 - m_dfEccMSQ;
	double dfRTEccMSQ = sqrt( dfBetaSQ );
	                                                                           
	// INITIALIZE LUNAR SOLAR TERMS
	m_dfDSDay = m_dfDaySince1950 + 18261.5;
    double dfXNODCE = fmod( 4.5236020 - 9.2422029e-4 * m_dfDSDay, g_dfTWOPI );
	double dfSTEM = sin( dfXNODCE );
	double dfCTEM = cos( dfXNODCE );
	double dfZCOSIL = .91375164 - .03568096 * dfCTEM;
	double dfZSINIL = sqrt( 1.0 - dfZCOSIL * dfZCOSIL );
	double dfZSINHL= .089683511 * dfSTEM / dfZSINIL;
	double dfZCOSHL = sqrt( 1.0 - dfZSINHL * dfZSINHL );
	m_dfDSGAM = 5.8351514 + .0019443680 * m_dfDSDay;
	
	double dfZX = .39785416 * dfSTEM / dfZSINIL;
	double dfZY = dfZCOSHL * dfCTEM + 0.91744867 * dfZSINHL * dfSTEM;
	dfZX = atan2( dfZX, dfZY );
	dfZX = m_dfDSGAM + dfZX - dfXNODCE;
	double dfZCOSGL = cos( dfZX );
	double dfZSINGL = sin( dfZX );                                                       

    double dfXNOI = 1.0 / m_dfMMM;

    BOOL bEnd = FALSE;

	double ZCOSG = 0.0, ZSING = 0.0, ZCOSI = 0.0, ZSINI = 0.0,
		ZCOSH = 0.0, ZSINH = 0.0, CC = 0.0;

	for( int i = 0; i < 2; i++ )
	{
		if( i == 0 )	// for the Sun
		{
			ZCOSG = dfZCOSGS;
			ZSING = dfZSINGS;
			ZCOSI = dfZCOSIS;
			ZSINI = dfZSINIS;
			ZCOSH = m_dfCosRAANM; 
			ZSINH = m_dfSinRAANM;
			CC = dfC1SS;
		}

		double A1 =  ZCOSG * ZCOSH + ZSING * ZCOSI * ZSINH;
   		double A3 = -ZSING * ZCOSH + ZCOSG * ZCOSI * ZSINH;
		double A7 = -ZCOSG * ZSINH + ZSING * ZCOSI * ZCOSH;
		double A8 =  ZSING * ZSINI;
		double A9 =  ZSING * ZSINH + ZCOSG * ZCOSI * ZCOSH;
		double A10 = ZCOSG * ZSINI;
		double A2 =  m_dfCosIncM * A7 + m_dfSinIncM * A8;
		double A4 =  m_dfCosIncM * A9 + m_dfSinIncM * A10;                                               
		double A5 = -m_dfSinIncM * A7 + m_dfCosIncM * A8;
		double A6 = -m_dfSinIncM * A9 + m_dfCosIncM * A10;
		                                                                            
		double X1 =  A1 * m_dfCosPerigeeM + A2 * m_dfSinPerigeeM;
		double X2 =  A3 * m_dfCosPerigeeM + A4 * m_dfSinPerigeeM;
		double X3 = -A1 * m_dfSinPerigeeM + A2 * m_dfCosPerigeeM;
		double X4 = -A3 * m_dfSinPerigeeM + A4 * m_dfCosPerigeeM;
		double X5 =  A5 * m_dfSinPerigeeM;
		double X6 =  A6 * m_dfSinPerigeeM;
		double X7 =  A5 * m_dfCosPerigeeM;
		double X8 =  A6 * m_dfCosPerigeeM;

		m_dfDSZ31 =  12.0 * X1 * X1 - 3.0 * X3 * X3;
		m_dfDSZ32 =  24.0 * X1 * X2 - 6.0 * X3 * X4;
		m_dfDSZ33 =  12.0 * X2 * X2 - 3.0 * X4 * X4;
		m_dfDSZ1  =  3.0 * ( A1 * A1 + A2 * A2 ) + m_dfDSZ31 * m_dfEccMSQ;
		m_dfDSZ2  =  6.0 * ( A1 * A3 + A2 * A4 ) + m_dfDSZ32 * m_dfEccMSQ;
		m_dfDSZ3  =  3.0 * ( A3 * A3 + A4 * A4 ) + m_dfDSZ33 * m_dfEccMSQ;
		m_dfDSZ11 = -6.0 * A1 * A5 + m_dfEccMSQ *( -24.0 * X1 * X7 - 6.0 * X3 * X5 );
		m_dfDSZ12 = -6.0 * ( A1 * A6 + A3 * A5 ) + 
			         m_dfEccMSQ * ( -24.0 * ( X2 * X7 + X1 * X8 ) - 6.0 * ( X3 * X6 + X4 * X5 ) );
		m_dfDSZ13 = -6.0 * A3 * A6 + m_dfEccMSQ * ( -24.0 * X2 * X8 - 6.0 * X4 * X6 );
		m_dfDSZ21 =  6.0 * A2 * A5 + m_dfEccMSQ * ( 24.0 * X1 * X5 - 6.0 * X3 * X7 );
		m_dfDSZ22 =  6.0 * ( A4 * A5 + A2 * A6 ) + 
			         m_dfEccMSQ *( 24.0 * ( X2 * X5 + X1 * X6 ) - 6.0 * ( X4 * X7 + X3 * X8 ) );
		m_dfDSZ23 =  6.0 * A4 * A6 + m_dfEccMSQ *( 24.0 * X2 * X6 - 6.0 * X4 * X8 );
		m_dfDSZ1  =  m_dfDSZ1 + m_dfDSZ1 + dfBetaSQ * m_dfDSZ31;
		m_dfDSZ2  =  m_dfDSZ2 + m_dfDSZ2 + dfBetaSQ * m_dfDSZ32;
		m_dfDSZ3  =  m_dfDSZ3 + m_dfDSZ3 + dfBetaSQ * m_dfDSZ33;
		m_dfDSS3  =  CC * dfXNOI;
		m_dfDSS2  = -.50 * m_dfDSS3 / dfRTEccMSQ;
		m_dfDSS4  =  m_dfDSS3 * dfRTEccMSQ;
		m_dfDSS1  = -15.0 * m_dfEccM * m_dfDSS4;
		m_dfDSS5  =  X1 * X3 + X2 * X4;
		m_dfDSS6  =  X2 * X3 + X1 * X4;
		m_dfDSS7  =  X2 * X4 - X1 * X3;                                                        
		
		if( bEnd ) break;

		// assign to the variables WRT the Sun
		m_dfDSSS1 = m_dfDSS1;
		m_dfDSSS2 = m_dfDSS2;
		m_dfDSSS3 = m_dfDSS3;
		m_dfDSSS4 = m_dfDSS4;
		m_dfDSSS5 = m_dfDSS5;
		m_dfDSSS6 = m_dfDSS6;
		m_dfDSSS7 = m_dfDSS7;
		m_dfDSSZ1 = m_dfDSZ1;
		m_dfDSSZ2 = m_dfDSZ2;
		m_dfDSSZ3 = m_dfDSZ3;
		m_dfDSSZ11 = m_dfDSZ11;
		m_dfDSSZ12 = m_dfDSZ12;
		m_dfDSSZ13 = m_dfDSZ13;
		m_dfDSSZ21 = m_dfDSZ21;
		m_dfDSSZ22 = m_dfDSZ22;
		m_dfDSSZ23 = m_dfDSZ23;
		m_dfDSSZ31 = m_dfDSZ31;
		m_dfDSSZ32 = m_dfDSZ32;
		m_dfDSSZ33 = m_dfDSZ33;

		// get the basic variables WRT the Moon
		ZCOSG = dfZCOSGL;
		ZSING = dfZSINGL;
		ZCOSI = dfZCOSIL;
		ZSINI = dfZSINIL;
		ZCOSH = dfZCOSHL * m_dfCosRAANM + dfZSINHL * m_dfSinRAANM;
		ZSINH = m_dfSinRAANM * dfZCOSHL - m_dfCosRAANM * dfZSINHL;
		CC = dfC1L;                                                                
		
		bEnd = TRUE;
	}

	return TRUE;
}

/******************************************************************************************

  Compute the secular variation of the deep space effect

  second part of FORTRAN subroutine DPPER

*******************************************************************************************/

BOOL cTLE2PosVel::DeepSpaceSecularEffect()
{
	m_dfDNDT = 0.0;
	double dfTHETA = fmod( m_dfGSTAtRefEpoch + m_dfTC * m_dfEarthRotationPerMinute, g_dfTWOPI );
	
	m_dfEccM += m_dfDEDT * m_dfTC;
    m_dfIncM += m_dfDIDT * m_dfTC;
	m_dfPerigeeM += m_dfDOMDT * m_dfTC;
	m_dfRAANM += m_dfDNODT * m_dfTC;
	m_dfMAM += m_dfDMDT * m_dfTC;
	
	if( m_dfIncM < 0.0 ) 
	{
		m_dfIncM = -m_dfIncM;
		m_dfPerigeeM -= g_dfPI;
		m_dfRAANM += g_dfPI;
	}

	if( m_nIREZ == 0 ) return TRUE;

	if( !DeepSapceResonance( dfTHETA ) ) return TRUE;                           
    
	m_dfMMM = m_dfMM0 + m_dfDNDT;
	
	return TRUE;
}

/******************************************************************************************

  Compute the resonant part of the deep space effect

  second part of FORTRAN subroutine SREZ

*******************************************************************************************/

BOOL cTLE2PosVel::DeepSapceResonance( double dfTHETA )
{
	double G22 = 5.7686396, G32 = 0.95240898;
	double G44 = 1.8014998, G52 = 1.0508330;
	double G54 = 4.4108898;

	double DELT=0.0, FT=0.0;
	int IRET = 0, IRETN = 0;
	double XNDT, XNDDT, XLDOT;
	double XL;

	if( m_dfATIME == 0.0 ) goto T390;
	if( m_dfTC >= 0.0 && m_dfATIME < 0.0 ) goto T390;
	if( m_dfTC < 0.0 && m_dfATIME >= 0.0 ) goto T390;

T351:
	if( fabs( m_dfTC ) < fabs( m_dfATIME ) )
	{
		if( m_dfTC >= 0.0 ) DELT = m_dfSTEPN;
		else DELT = m_dfSTEPP;

		IRET = 351;
		goto T380;
	}

	if( m_dfTC < 0.0 ) DELT = m_dfSTEPN;
	else DELT = m_dfSTEPP;

T356:
	if( fabs( m_dfTC - m_dfATIME ) >= m_dfSTEPP ) 
	{
		IRET = 0;
		goto T380;
	}
	else
	{
		FT = m_dfTC - m_dfATIME;
		IRETN = 0;
		goto T370;
	}

T361:
	m_dfMMM = m_dfXNI + XNDT * FT + XNDDT * FT * FT *.5;
	XL = m_dfXLI + XLDOT * FT + XNDT * FT * FT *.5;
	if( m_nIREZ == 2 )
	{
		m_dfMAM = XL - 2.0 * m_dfRAANM + 2.0 * dfTHETA;
		m_dfDNDT = m_dfMMM - m_dfMM0;
		return TRUE;
	}
	else
	{
		m_dfMAM = XL - m_dfRAANM - m_dfPerigeeM + dfTHETA;                                           
		m_dfDNDT = m_dfMMM - m_dfMM0;
		return TRUE;
	}

	// DOT TERMS CALCULATED

T370:
	if( m_nIREZ == 1 ) //NEAR - SYNCHRONOUS RESONANCE TERMS
	{
		XNDT = m_dfDEL1 * sin( m_dfXLI - m_dfFASX2 ) + 
		       m_dfDEL2 * sin( 2.0 * ( m_dfXLI - m_dfFASX4 ) ) +
               m_dfDEL3 * sin( 3.0 * ( m_dfXLI - m_dfFASX6 ) );                                      
		XLDOT = m_dfXNI + m_dfXFACT;
		XNDDT = m_dfDEL1 * cos( m_dfXLI - m_dfFASX2 ) +
			    2.0 * m_dfDEL2 * cos( 2.0 * ( m_dfXLI - m_dfFASX4 ) ) +
				3.0 * m_dfDEL3 * cos( 3.0 * ( m_dfXLI - m_dfFASX6 ) );
		XNDDT = XNDDT * XLDOT;
	}
	else // NEAR - HALF-DAY RESONANCE TERMS                                       
	{
		double XOMI = m_dfPerigee0 + m_dfPerigeeDot * m_dfATIME;
		double X2OMI = XOMI + XOMI;
		double X2LI = m_dfXLI + m_dfXLI;
		XNDT = m_dfD2201 * sin( X2OMI + m_dfXLI - G22 ) +
               m_dfD2211 * sin( m_dfXLI - G22 ) +
			   m_dfD3210 * sin( XOMI + m_dfXLI - G32 ) +
			   m_dfD3222 * sin( -XOMI + m_dfXLI - G32 ) +
			   m_dfD4410 * sin( X2OMI + X2LI - G44 ) +
			   m_dfD4422 * sin( X2LI - G44 ) +
			   m_dfD5220 * sin( XOMI + m_dfXLI - G52 ) +
			   m_dfD5232 * sin( -XOMI + m_dfXLI - G52 ) +
			   m_dfD5421 * sin( XOMI + X2LI - G54 ) +
			   m_dfD5433 * sin( -XOMI + X2LI - G54 );
		XLDOT = m_dfXNI + m_dfXFACT;
		XNDDT = m_dfD2201 * cos( X2OMI + m_dfXLI - G22 ) +
			    m_dfD2211 * cos( m_dfXLI - G22 ) +
				m_dfD3210 * cos( XOMI + m_dfXLI - G32 ) +
				m_dfD3222 * cos( -XOMI + m_dfXLI - G32 ) +
				m_dfD5220 * cos( XOMI + m_dfXLI - G52 ) +
				m_dfD5232 * cos( -XOMI + m_dfXLI - G52 ) +
				2.0 * ( m_dfD4410 * cos( X2OMI + X2LI - G44 ) +
				        m_dfD4422 * cos( X2LI - G44 ) +
						m_dfD5421 * cos( XOMI + X2LI - G54 ) +
						m_dfD5433 * cos( -XOMI + X2LI - G54 )   );
		XNDDT = XNDDT * XLDOT;                                                     
	}
	
	if( IRETN == 381 ) goto T381;
	goto T361;                                                             
                                                                            
//	INTEGRATOR                                                            

T380:
	IRETN = 381;                                                             
	goto T370;

T381:
	m_dfXLI = m_dfXLI + XLDOT * DELT + XNDT * m_dfSTEP2; 
	m_dfXNI = m_dfXNI + XNDT * DELT + XNDDT * m_dfSTEP2;
	m_dfATIME += DELT;
	
	if( IRET == 351 ) goto T351;
	goto T356;                                                             
                                                                            
T390:
	if( m_dfTC >= 0.0 ) DELT = m_dfSTEPP;
	else DELT = m_dfSTEPN;
	
	m_dfATIME = 0.0;
	m_dfXNI = m_dfMM0;
	m_dfXLI = m_dfXLAMO;

	goto T356;

	//return TRUE; 20170313 lei
}

/******************************************************************************************

  Compute the periodic part of the deep space effect

  second part of FORTRAN subroutine DSPACE

*******************************************************************************************/

BOOL cTLE2PosVel::DeepSpacePeriodicEffect( double &dfEccP, double &dfIncP, double &dfRAANP,
										   double &dfPerigeeP, double &dfMAP )
{
	// for the Sun
	double ZM = m_dfZMOS + m_dfZNS * m_dfTC;
	double ZF = ZM + 2.0 * m_dfZES * sin( ZM );
	double SINZF = sin( ZF );
	double F2 = .5 * SINZF * SINZF -.25;
	double F3 =-.5 * SINZF * cos( ZF );
	double SES = m_dfDSSE2 * F2 + m_dfDSSE3 * F3;                                                     
    double SIS = m_dfDSSI2 * F2 + m_dfDSSI3 * F3;
	double SLS = m_dfDSSL2 * F2 + m_dfDSSL3 * F3 + m_dfDSSL4 * SINZF;
	double SGHS = m_dfDSSGH2 * F2 + m_dfDSSGH3 * F3 + m_dfDSSGH4 * SINZF;
	double SHS = m_dfDSSH2 * F2 + m_dfDSSH3 * F3;
	
	// for the Moon
	ZM = m_dfZMOL + m_dfZNL * m_dfTC;
	ZF = ZM + 2.0 * m_dfZEL * sin( ZM );
	SINZF = sin( ZF );
	F2 = .5 * SINZF * SINZF - .25;
	F3 = -.5 * SINZF * cos( ZF );
	
	double SEL = m_dfDSEE2 * F2 + m_dfDSE3 * F3;
	double SIL = m_dfDSXI2 * F2 + m_dfDSXI3 * F3;
	double SLL = m_dfDSXL2 * F2 + m_dfDSXL3 * F3 + m_dfDSXL4 * SINZF;
	double SGHL = m_dfDSXGH2 * F2 + m_dfDSXGH3 * F3 + m_dfDSXGH4 * SINZF;
	double SHL = m_dfDSXH2 * F2 + m_dfDSXH3 * F3;
	
	double PE = SES + SEL;
	double PINC = SIS + SIL;
	double PL = SLS + SLL;
	double PGH = SGHS + SGHL;
	double PH = SHS + SHL;                                                            
	
	dfIncP += PINC;
	dfEccP += PE;
	double SINIP = sin( dfIncP );
	double COSIP = cos( dfIncP );
	
	if( dfIncP > 0.2 ) //APPLY PERIODICS DIRECTLY
	{
		PH = PH / SINIP;
		PGH = PGH - COSIP * PH;
		dfPerigeeP += PGH;
		dfRAANP += PH;
		dfMAP +=  PL;                                                          
      
		return TRUE;                                                                
	}
                                                                            
	// APPLY PERIODICS WITH LYDDANE MODIFICATION 
	double SINOP = sin( dfRAANP );
	double COSOP = cos( dfRAANP );
	double ALFDP = SINIP * SINOP;
	double BETDP = SINIP * COSOP;
	double DALF = PH * COSOP + PINC * COSIP * SINOP;
	double DBET = -PH * SINOP + PINC * COSIP * COSOP;
	ALFDP = ALFDP + DALF;
	BETDP = BETDP + DBET;                                                      
    dfRAANP = fmod( dfRAANP, g_dfTWOPI );                                                 
    double XLS = dfMAP + dfPerigeeP + COSIP * dfRAANP;
	double DLS = PL + PGH - PINC * dfRAANP * SINIP;
	XLS = XLS + DLS;
	double XNOH = dfRAANP;
	dfRAANP = atan2( ALFDP, BETDP );
	if( fabs( XNOH - dfRAANP ) > g_dfPI )
	{
		if( dfRAANP < XNOH ) dfRAANP += g_dfTWOPI;
		else dfRAANP -= g_dfTWOPI;                                            
	}
	
	dfMAP += PL;
	dfPerigeeP = XLS - dfMAP - COSIP * dfRAANP;    
	
	return TRUE;
}

/**********************************************************************************************

  Convert the satellite inertial position and velocity to EFEC system

  dfJD: given time epoch in JD
  pdfPos: as input is the satellite inertial position in meter,
          as output is the satellite ECEF position in meter
  pdfvel: as input is the satellite inertial velocity in meter/s,
          as output is the satellite ECEF velocity in meter/s

**********************************************************************************************/

BOOL cTLE2PosVel::FromInertialToEFEC( double dfJD, double *pdfPos, double *pdfVel )
{
	double dfGST = cGreenwichST::ComputeGST( dfJD );
	double dfSin = sin( dfGST ), dfCos = cos( dfGST );
	double pdfP[ 3 ], pdfV[ 3 ];

	pdfP[ 0 ] =  pdfPos[ 0 ] * dfCos + pdfPos[ 1 ] * dfSin;
	pdfP[ 1 ] = -pdfPos[ 0 ] * dfSin + pdfPos[ 1 ] * dfCos;
	pdfP[ 2 ] =  pdfPos[ 2 ];

	pdfV[ 0 ] =  pdfVel[ 0 ] * dfCos + pdfVel[ 1 ] * dfSin + pdfP[ 1 ] * g_dfEarthAngVelocity;
	pdfV[ 1 ] = -pdfVel[ 0 ] * dfSin + pdfVel[ 1 ] * dfCos - pdfP[ 0 ] * g_dfEarthAngVelocity;
	pdfV[ 2 ] =  pdfVel[ 2 ];

	for( int i = 0; i < 3; i++ )
	{
		pdfPos[ i ] = pdfP[ i ];
		pdfVel[ i ] = pdfV[ i ];
	}

	return TRUE;
}



/***********************************************************************************

  Read Line 1 of TLE element

************************************************************************************/

BOOL cTLE2PosVel::ReadTLELine1( string line, stSatelliteIOE &stIOE )
{
//	if( !CheckSum( line ) ) return FALSE;

	string sub;

	// norad number of the object
	sub = line.substr( 1, 6 );
	stIOE.nSatelliteID = atoi( sub.c_str() );

	// ref date/time
	sub = line.substr( 18, 2 );
	int nYear = atoi( sub.c_str() );
	if( nYear < 90 ) nYear += 2000;
	else nYear += 1900;

	sub = line.substr( 20, 12 );
	double dfDay = atof( sub.c_str() );

	double dfRefJD;
	g_DateTimeZ.DateTime2JD( nYear, 1, 0, 0, 0, 0.0, dfRefJD );

	dfRefJD += dfDay;

	stIOE.SetRefJD( dfRefJD );
//	stIOE.nIntJD = (int) dfRefJD;
//	stIOE.nFractionJD = (int)( ( dfRefJD - stIOE.nIntJD ) * 1.0e9 );

	// n dot -- first time derivative of the mean motion or
	// ballistic coefficient (depending on ephemeris type)
	sub = line.substr( 33, 10 );
	stIOE.pfElement7to18[ 0 ] = (float) atof( sub.c_str() );

	// n dot dot - second time derivative of the mean motion
	sub = line.substr( 44, 6 );
	double dfTemp = atoi( sub.c_str() ) * 1.0e-5;
	sub = line.substr( 50, 2 );
	int nTemp = atoi( sub.c_str() );
		
	stIOE.pfElement7to18[ 1 ] = (float) ( dfTemp * pow( 10.0, (double) nTemp ) );

	// BSTAR drag term if GP4 general perturbation theory is used
	sub = line.substr( 53, 6 );
	dfTemp = atoi( sub.c_str() ) * 1.0e-5;
	sub = line.substr( 59, 2 );
	nTemp = atoi( sub.c_str() );
		
	stIOE.pfElement7to18[ 2 ] = (float) (dfTemp * pow( 10.0, (double) nTemp ) );

	// Ephemeris type
	sub = line.substr( 62, 1 );
	
	return TRUE;
}


/***********************************************************************************

  Read Line 2 of TLE element

************************************************************************************/

BOOL cTLE2PosVel::ReadTLELine2( string line, stSatelliteIOE &stIOE )
{
//	if( !CheckSum( line ) ) return FALSE;

	string sub;

	// Inlcination 
	sub = line.substr( 8, 8 );
	stIOE.pnElement1to6[ 2 ] = (int)( atof( sub.c_str() ) * 1.0e5 );

	// RA of ascending node
	sub = line.substr( 17, 8 );
	stIOE.pnElement1to6[ 3 ] = (int)( atof( sub.c_str() ) * 1.0e5 );

	// Eccentricity 
	sub = line.substr( 26, 7 );
	stIOE.pnElement1to6[ 1 ] = atoi( sub.c_str() );

	// Perigee argument
	sub = line.substr( 34, 8 );
	stIOE.pnElement1to6[ 4 ] = (int)( atof( sub.c_str() ) * 1.0e5 );

	// Mean anomaly
	sub = line.substr( 43, 8 );
	stIOE.pnElement1to6[ 5 ] = (int)( atof( sub.c_str() ) * 1.0e5 );

	// Mean motion orbits per day
	sub = line.substr( 52, 11 );
	stIOE.pnElement1to6[ 0 ] = (int)( atof( sub.c_str() ) * 1.0e8 );

	// orbital period in days
	stIOE.pfElement7to18[ 11 ] = (float) (	1.0 / ( stIOE.pnElement1to6[ 0 ] * 1.0e-8 ) );

	stIOE.cElementType = 'T';
	
	return TRUE;
}


/******************************************************************************

  Check Sum verification of the given line 

******************************************************************************/

BOOL cTLE2PosVel::CheckSum( string line )
{
	char str[ 100 ];

	strcpy( str, line.c_str() );

	int CheckSum = 0;

	for( int nn = 0; nn < 68; nn++ )
	{
		int ascii = __toascii( str[ nn ] );

		if( ascii >= 48 && ascii <= 57 ) ascii = ascii - 48;
		else if( ascii == 32 || ascii == 46 || ascii == 43 || ascii == 45 ) ascii = 0;
		
		CheckSum += ascii;
	}

	CheckSum = CheckSum - CheckSum / 10 * 10;

	int nGivenCheckSum = __toascii( str[ 68 ] );

	nGivenCheckSum = nGivenCheckSum - 48;

	if( CheckSum != nGivenCheckSum ) return FALSE;

	return TRUE;
}

/******************************************************************************

 This function can only be called after successfully initialising TLE. The 
 heights are only approximate.

******************************************************************************/

void cTLE2PosVel::GetPerigeeApogeeHeights( double &dfPerigeeHeight, 
										   double &dfApogeeHeight )
{
	dfPerigeeHeight = ( m_dfQ0 - 1.0 ) * m_dfEarthRadius;
	dfApogeeHeight = ( m_dfQ1 - 1.0 ) * m_dfEarthRadius;
}


void cTLE2PosVel::GetTLE( double *pdfTLE )
{
	for( int i = 0; i < 10; i++ ) pdfTLE[ i ] = m_pdfTLE[ i ];
}


/******************************************************************************************
******************************************************************************************
******************************************************************************************
******************************************************************************************
******************************************************************************************
*****************************************************************************************/

/******************************************************************************************

  Added by Ricky OBrien on 1/6/04 to read all of the TLE data and return a vector containing this
  information

*******************************************************************************************/

BOOL cTLE2PosVel::ReadAllTLE( std::vector<stSatelliteIOE> &TLEData, const char *szFileName, 
							  BOOL bPerigeeTest, double perigeeLimit, double apogeeLimit )
{
	try
	{
		// check whether the given file exists
		//TLEData.clear( );
		//TLEData.reserve( 100 );

		WIN32_FIND_DATA FindFileData;

		HANDLE hFile = FindFirstFile( szFileName, &FindFileData );

		if( hFile == INVALID_HANDLE_VALUE ) 
		{
			return FALSE;
		}

		FindClose( hFile );

		string line;
		ifstream in( szFileName );
		
		stSatelliteIOE tIOE;

		while( getline( in, line ) )	// 
		{
			try
			{
				if( line[0] == '#' )
				{
					do
					{
						getline( in, line );	
					}	while( line[0] == '#' );
				}

				//if( line[0] != '1' )
				string temp = line.substr(0, 2);
				if (strcmp(temp.c_str(), "1 ") != 0)
				{
					tIOE.Line1 = line.c_str( );  // added RO
					// Line 1
					getline( in, line );
				}
				else
				{
					tIOE.Line1.clear( );
				}
								
				if( line[0] == '1' && line[1] == ' ')
				{
					tIOE.Line2 = line.c_str( );  // added RO
					if( !ReadTLELine1( line, tIOE ) ) 
					{
						// error reading the data
						return FALSE;					
					}
				}
				else
				{
					//eosLOG( "error reading TLE file ");
					// end of file
					return TRUE;
				}
			
				// Line 2
				getline( in, line );
				
				if( line[0] == '2' && line[1] == ' ')
				{
					tIOE.Line3 = line.c_str( );  // added RO
					if( !ReadTLELine2( line, tIOE ) ) 
					{
						// error reading the data
						//eosLOG( "error reading tle file 2");
						return FALSE;
					}
				}

				if( bPerigeeTest )
				{
                    double sm = tIOE.GetSemiMajor();
					double ecc = tIOE.GetEccentricity();

					if( sm * ( 1.0 - ecc ) > perigeeLimit && 
						sm * ( 1.0 + ecc ) < apogeeLimit    ) TLEData.push_back( tIOE );
				}
				else TLEData.push_back( tIOE );
			}
			catch( ... )
			{
			}
		}

		return TRUE;
	}
	catch( ... )
	{
		throw"ReadAllTLE() error ...";
		return false;
	}

}

BOOL cTLE2PosVel::ReadAllTLE( std::vector<stSatelliteIOE> &TLEData, const char *szFileName, int numberID, int*noradID )
{
	try
	{
		// check whether the given file exists
		//TLEData.clear( );
		//TLEData.reserve( 100 );

		WIN32_FIND_DATA FindFileData;

		HANDLE hFile = FindFirstFile( szFileName, &FindFileData );

		if( hFile == INVALID_HANDLE_VALUE ) 
		{
			return FALSE;
		}

		FindClose( hFile );

		string line;
		ifstream in( szFileName );
		
		stSatelliteIOE tIOE;

		while( getline( in, line ) )	// 
		{
			try
			{
				if( line[0] == '#' )
				{
					do
					{
						getline( in, line );	
					}	while( line[0] == '#' );
				}

				//if( line[0] != '1' )
				string temp = line.substr(0, 2);
				if (strcmp(temp.c_str(), "1 ") != 0)
				{
					tIOE.Line1 = line.c_str( );  // added RO
					// Line 1
					getline( in, line );
				}
				else
				{
					tIOE.Line1.clear( );
				}
								
				if( line[0] == '1' && line[1] == ' ')
				{
					tIOE.Line2 = line.c_str( );  // added RO
					if( !ReadTLELine1( line, tIOE ) ) 
					{
						// error reading the data
						return FALSE;					
					}
				}
				else
				{
					//eosLOG( "error reading TLE file ");
					// end of file
					return TRUE;
				}
			
				// Line 2
				getline( in, line );
				
				if( line[0] == '2' && line[1] == ' ')
				{
					tIOE.Line3 = line.c_str( );  // added RO
					if( !ReadTLELine2( line, tIOE ) ) 
					{
						// error reading the data
						//eosLOG( "error reading tle file 2");
						return FALSE;
					}
				}

				for( int k = 0; k < numberID; k++ )
				{
					int id = tIOE.nSatelliteID;
					if( id == noradID[ id ] )
					{
						TLEData.push_back( tIOE );
						break;
					}
				}				
			}
			catch( ... )
			{
			}
		}

		return TRUE;
	}
	catch( ... )
	{
		throw"ReadAllTLE() error ...";
		return false;
	}
}


/*
BOOL cTLE2PosVel::ReadAllTLE( std::vector<stSatelliteIOE> &TLEData, const char *szFileName )
{
	try
	{
		m_nCounter++;
		if( m_nCounter == 1 )
		{
			for( int i = 0; i < 50000; ++i )
			{
				m_bTLEExists[i] = false;
			}
		}
		// check whether the given file exists
		//TLEData.clear( );
		//TLEData.reserve( 100 );

		WIN32_FIND_DATA FindFileData;

		HANDLE hFile = FindFirstFile( szFileName, &FindFileData );

		if( hFile == INVALID_HANDLE_VALUE ) 
		{
			return FALSE;
		}

		FindClose( hFile );

		string line;
		ifstream in( szFileName );
		
		double increment = 25;
		double Largest[ 1601 ];
		double Smallest[ 1601 ];
		for( int i = 0; i < 1601; ++i )
		{
			Largest[i] = 0;
			Smallest[i] = 0;
		}

		stSatelliteIOE tIOE;

		while( getline( in, line ) )	// 
		{
			try
			{
				if( line[0] == '#' )
				{
					do
					{
						getline( in, line );	
					}	while( line[0] == '#' );
				}

				if( line[0] != '1' )
				{
					tIOE.Line1 = line.c_str( );  // added RO
					// Line 1
					getline( in, line );
				}
				else
				{
					tIOE.Line1.clear( );
				}
								
				if( line[0] == '1' )
				{
					tIOE.Line2 = line.c_str( );  // added RO
					if( !ReadTLELine1( line, tIOE ) ) 
					{
						// error reading the data
						return FALSE;					
					}
				}
				else
				{
					//eosLOG( "error reading TLE file ");
					// end of file
					return TRUE;
				}
			
				// Line 2
				getline( in, line );
				
				if( line[0] == '2' )
				{
					tIOE.Line3 = line.c_str( );  // added RO
					if( !ReadTLELine2( line, tIOE ) ) 
					{
						// error reading the data
						//eosLOG( "error reading tle file 2");
						return FALSE;
					}
				}
				else
				{
					//eosLOG( "error reading tle file 3");
					// end of file
					return TRUE;
				}

				double largestAlt = 0;
				double altitude = 0;//GetSmallestAltitude( tIOE, largestAlt );
				
				if( largestAlt >= 0 && largestAlt < 40000 )
				{
					++Largest[ (int) ( largestAlt/increment ) ];
				}
				
				if( largestAlt >= 0 && largestAlt < 40000 )
				{
					++Smallest[ (int) ( altitude/increment ) ];
				}

				//eosLOG( "rev " << tIOE.pnElement1to6[0]*1e-8 );
/*				if( tIOE.pnElement1to6[0]*1e-8 > 4 /*|| (tIOE.pnElement1to6[0]*1e-8 > 7)*/ /*********)*/
/*				{
//					stSatelliteIOE tIOEA = tIOE;
//					double largestAlt;
//					double altitude = GetSmallestAltitude( tIOEA, largestAlt );
/*					double incOriginal = tIOE.pnElement1to6[2]*1e-5;
					double ecc = (double) tIOE.pnElement1to6[ 1 ] * 1.0e-7;
					if( 1==1 || (altitude < 2000 && altitude > 400
						&& largestAlt < 2000 && largestAlt > 400 
						/*&& fabs( ecc ) < 0.2 && TLEData.size( ) < 500*/ 
/*						&& tIOE.GetRefJD( ) >= 2451940) )
					{
						
						if( m_nCounter == 1 )
						{
							m_bTLEExists[ tIOE.nSatelliteID ] = true;
						}
						if( m_nCounter == 2 )
						{
							int nStop = 100;
						}
					
/*						if( /*( nCounter > 1 && bTLEExists[ tIOE.nSatelliteID ] == true ) 
//							|| nCounter == 1*/ /*****1==1******/ /***************)
/*						{
							//if( fabs( tIOE.pnElement1to6[2]*1e-5 ) < 60 /*&& fabs( tIOE.pnElement1to6[2]*1e-5 ) > 35*/ /********** )
							//{
//								TLEData.push_back( tIOE );
							//}
//							PerturbElement( tIOEA, 100000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 150000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 200000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 250000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 300000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 350000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 400000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 450000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;



//							PerturbElement( tIOEA, 500000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 550000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 600000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 650000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//
//							PerturbElement( tIOEA, 700000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 750000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 800000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 850000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 900000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 950000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 1000000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 1050000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 1100000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 1150000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 1200000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
//							PerturbElement( tIOEA, 1250000 );
//							TLEData.push_back( tIOEA );
//							tIOEA = tIOE;
/*
						}
					}

				}
			}
			catch( ... )
			{
				eosTHROW;
			}
		}
//		FILE *stream = fopen( "D:\\Logs\\Scheduler\\Altitudes.dat", "w" );
//		for( i = 0; i < 1601; i++ )
//		{
//			fprintf( stream, "%f , %f , %f\n", i*increment, Smallest[i], Largest[i ] );
//		}
//		fflush( stream );
//		fclose( stream );
//		eosLOG( "TLE Size " << TLEData.size( ) );
		return TRUE;
	}
	catch( ... )
	{
		eosTHROW;
	}
}


*/