/*************************************************************************** 

 This class is basically a translation of the FORTRAN program TLE2IRV

 Author: Jizhang Sang, Feb-Apr, 2002

 (C) COPYRIGHT -Electro Optic Systems Australia (EOS) 2002,
 All rights reserved. No part of this program may be photocopied, 
 reproduced or translated to another program language without the prior
 written consent of EOS

 Reference: SGP4 Algorithm

***************************************************************************/

#if!defined INC_TLE2POSVEL
#define INC_TLE2POSVEL

#include <windows.h>
#include "DataStructure.h"

#include <string>
#include <fstream>
#include <vector> 

using namespace std;

class cTLE2PosVel
{
	// original TLE elements
	// pdfTLE[ 1 ]: year
	// pdfTLE[ 2 ]: day of year
	// pdfTLE[ 3 ]: bstar
	// pdfTLE[ 4 ]: Eccentricity
	// pdfTLE[ 5 ]: inclination, in deg
	// pdfTLE[ 6 ]: RAAN in deg
	// pdfTLE[ 7 ]: argument of perigee in deg
	// pdfTLE[ 8 ]: mean anomaly in deg
	// pdfTLE[ 9 ]: mean motion in orbits per day

	double m_pdfTLE[ 10 ];

	// constants
	double m_dfEarthRadius, m_dfMeanEarthRadius, m_dfEarthFlattening, m_dfEarthRotation;
	double m_dfJD1950;

	// FK5 CONSTANT VALUE
    double m_dfC1, m_dfC1P2P, m_dfTHGR70, m_dfFK5R;

	// WGS-72 constant
    double m_dfEarthGM, m_dfEarthRadiusKM, m_dfXPDOTP, m_dfEarthRotationPerMinute;
    double m_dfJ2, m_dfJ3, m_dfJ4, m_dfJ3OJ2;                               
	double m_dfVelocityChange;

	// other
	double m_dfMinutesPerDay, m_dfMinSatGeoDis, m_dfQZMS2T;

	// reference epoch of given mean orbital element
	double m_dfRefJD, m_dfDaySince1950;
	// mean orbital elements at reference from TLE elements
	double m_dfSM0, m_dfEcc0, m_dfInc0, m_dfRAAN0, m_dfPerigee0, m_dfMA0;
	double m_dfMM0;	// mean motion
	double m_dfBSTAR;

	// AUXILLARY EPOCH QUANTITIES 
	double m_dfEcc0SQ, m_dfOMEcc0SQ, m_dfRTEcc0SQ, m_dfCosI0, m_dfCosI02, m_dfSinI0;
	double m_dfP0, m_dfP0SQ, m_dfQ0, m_dfQ1;
	double m_dfCon42, m_dfCon41, m_dfSMInv, m_dfEccInv;
	double m_dfGSTAtRefEpoch;
	double m_df2O3;

	BOOL m_bDeepSpace;		// TRUE if the orbit is 5876km high or above
	BOOL m_bOrbitLowHigh;	// TRUE if orbit is lower than 220km or above 5876km

	int m_nError;

	double m_dfEta;			// Eq (13) of SGP4 Algorithm

	double m_dfX1MTH2, m_dfX7THM1;	// 1-cos(Inc)**2, 7cos(i)**2-1

	// rate due to Earth asphericity
	double m_dfMMDot, m_dfPerigeeDot, m_dfRAANDot, m_dfXPIDOT;	// m_dfXPIDOT = m_dfPerigeeDot + m_dfRAANDot

	// for short periodic terms
	double m_dfXLCOF, m_dfAYCOF, m_dfDeltaMA0, m_dfSinMAO, m_dfTC;

	// ballistic terms
	double m_dfPerigeeDotDrag, m_dfMMDotDrag, m_dfRAANDot2Drag;
	double m_dfD2, m_dfD3, m_dfD4, 
		   m_dfT2Coe, m_dfT3Coe, m_dfT4Coe, m_dfT5Coe;	// coefficients for deltaT^2, ..., terms of MA+Perigee+RAAN
	double m_dfCC1, m_dfCC2, m_dfCC3, m_dfCC4, m_dfCC5;

	// mean value at epoch
	double m_dfSMM, m_dfEccM, m_dfIncM, m_dfRAANM, m_dfPerigeeM, m_dfMAM;
	double m_dfMMM;	// mean motion

	double m_dfSinRAANM, m_dfCosRAANM, m_dfSinPerigeeM, m_dfCosPerigeeM, m_dfSinIncM, m_dfCosIncM;
	double m_dfEccMSQ;

	// Deep Space Effect
	
	// constant used
	double m_dfZNS, m_dfZES, m_dfZNL, m_dfZEL;
	
	double m_dfDSS1, m_dfDSS2, m_dfDSS3, m_dfDSS4, m_dfDSS5, m_dfDSS6, m_dfDSS7,
		   m_dfDSZ1, m_dfDSZ2, m_dfDSZ3, 
		   m_dfDSZ11, m_dfDSZ12, m_dfDSZ13,
		   m_dfDSZ21, m_dfDSZ22, m_dfDSZ23, 
		   m_dfDSZ31, m_dfDSZ32, m_dfDSZ33,
		   m_dfDSSS1, m_dfDSSS2, m_dfDSSS3, m_dfDSSS4, m_dfDSSS5, m_dfDSSS6, m_dfDSSS7,
		   m_dfDSSZ1, m_dfDSSZ2, m_dfDSSZ3, 
		   m_dfDSSZ11, m_dfDSSZ12, m_dfDSSZ13,
		   m_dfDSSZ21, m_dfDSSZ22, m_dfDSSZ23,
		   m_dfDSSZ31, m_dfDSSZ32, m_dfDSSZ33,
		   m_dfDSDay, m_dfDSGAM;

	double m_dfZMOL, m_dfZMOS;
	double m_dfDSSE2, m_dfDSSE3, m_dfDSSI2, m_dfDSSI3, m_dfDSSL2, m_dfDSSL3, m_dfDSSL4,
		   m_dfDSSGH2, m_dfDSSGH3, m_dfDSSGH4, m_dfDSSH2, m_dfDSSH3;

	double m_dfDSEE2, m_dfDSE3, m_dfDSXI2, m_dfDSXI3, m_dfDSXL2, m_dfDSXL3, m_dfDSXL4,
		   m_dfDSXGH2, m_dfDSXGH3, m_dfDSXGH4, m_dfDSXH2, m_dfDSXH3;

	int m_nIREZ;	// resonance type, 1 for 24 hour resonance, 2 for 12 hour resonance
	double m_dfDEDT, m_dfDIDT, m_dfDMDT, m_dfDOMDT, m_dfDNODT;
	double m_dfDNDT;

	double m_dfD2201, m_dfD2211, m_dfD3210, m_dfD3222, m_dfD4410, m_dfD4422,
		   m_dfD5220, m_dfD5232, m_dfD5421, m_dfD5433;
	double m_dfDEL1, m_dfDEL2, m_dfDEL3, m_dfFASX2, m_dfFASX4, m_dfFASX6, m_dfXLAMO;
	double m_dfXLI, m_dfXNI, m_dfATIME, m_dfSTEPP, m_dfSTEPN, m_dfSTEP2, m_dfXFACT;
	
	BOOL m_bInit;
	BOOL m_bPosOnly;

public:

	cTLE2PosVel();
	~cTLE2PosVel();

	BOOL ReadTLE( int nNORADID, char *szFileName );

	BOOL SetOrbitalElements( stSatelliteIOE stIOE );
	BOOL SetOrbitalElements( double *pdfTLE );
	BOOL GetOrbitalElements( double *pdfTLE );
	void GetOrbitalElementsRefJD( double &dfRefJD );
	void SetComputePositionOnly( BOOL bState ) { m_bPosOnly = bState; }

	BOOL ComputeInertialPosVel( double dfJD, double *pdfPos, double *pdfVel );
	BOOL ComputeInertialPosVel( double &dfIntJDUTCGivenEpoch, double &dfFraJDUTCGivenEpoch, 
			                    double *pdfPos, double *pdfVel );
	BOOL ComputeECEFPosVel( double dfJD, double *pdfPos, double *pdfVel );

	void GetPerigeeApogeeHeights( double &dfPerigeeHeight, double &dfApogeeHeight );
	double GetInclination() { return m_dfInc0; }

	void GetTLE( double *pdfTLE );

	BOOL ReadAllTLE( std::vector<stSatelliteIOE> &TLEData, const char *szFileName, BOOL bPerigeeTest = FALSE,
		             double perigeeLimit = 6378137.0 + 250000.0, double apogeeLimit = 6378137.0 + 5000000.0 );
	BOOL ReadAllTLE( std::vector<stSatelliteIOE> &TLEData, const char *szFileName, int numberID, int*NORADID );

private:

	BOOL Initialise();
	BOOL SGP4( double dfJD, double *pdfPos, double *pdfVel );

	BOOL InitialiseDeepSpace();

	BOOL DSCOM();

	BOOL DeepSpaceSecularEffect();
	BOOL DeepSapceResonance( double dfTHETA );
	BOOL DeepSpacePeriodicEffect( double &dfEccP, double &dfIncP, double &dfRAANP,
								  double &dfPerigeeP, double &dfMAP );
	BOOL FromInertialToEFEC( double dfJD, double *pdfPos, double *pdfVel );

	BOOL ReadTLELine1( string line, struct stSatelliteIOE &stIOE );
	BOOL ReadTLELine2( string line, struct stSatelliteIOE &stIOE );

	BOOL CheckSum( string line );

	int m_nCounter;
	bool m_bTLEExists[50000];

	double GetAltitudeKM( const double pos[3] );

	void PerturbElement( stSatelliteIOE &tIOE, int nsatID );


};

#endif