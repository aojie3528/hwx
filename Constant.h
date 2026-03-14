/*********************************************************************************
 
 mathematical and physical constants used in this orbit determination software


**********************************************************************************/

#if!defined INC_ORBITCONSTANT
#define INC_ORBITCONSTANT

// angles conversion
					
const double g_dfPI = 3.14159265358979;
const double g_dfTWOPI = g_dfPI * 2.0;
const double g_dfHALFPI = g_dfPI / 2.0;

const double g_dfDEG2RAD = g_dfPI / 180.0;
const double g_dfHOUR2RAD = 15.0 * g_dfDEG2RAD;
const double g_dfRAD2DEG = 180.0 / g_dfPI;

const double g_dfSEC2RAD = g_dfDEG2RAD / 3600.0;
const double g_dfRAD2SEC = g_dfRAD2DEG * 3600.0;

// IERS constants
const double g_dfJ2000 = 2451545.0;
const double g_dfMJDRef = 2400000.5;

const double g_dfLightSpeed = 299792458.0;	// m/s

const double g_dfAstroUnitSecond = 499.0047838061;	// astronomical unit, in second of time
const double g_dfAstroUnitMeter = 149597870691.0;	// astronomical unit, in meter

const double g_dfObliquityJ2000SEC = 84381.412;		// Obliquity of the ecliptic at J2000, in seconds
const double g_dfObliquityJ2000RAD = g_dfObliquityJ2000SEC * g_dfSEC2RAD; // Obliquity of the ecliptic at J2000, in radians	

const double g_dfEarthSemiMajor = 6378136.49;	// meter
const double g_dfEarthFlattening = 1.0 / 298.25642;
const double g_dfEarthEccentricity = 2.0 * g_dfEarthFlattening - g_dfEarthFlattening * g_dfEarthFlattening;

const double g_dfEarthGM = 3.986004418e14; 
const double g_dfSolarGM = 1.327124e20;

const double g_dfMoonEarthMassRatio = 0.0123000345;
const double g_dfEarthMoonMassRatio = 1.0 / g_dfMoonEarthMassRatio;

const double g_dfEarthAngVelocity = 7.292115e-5;	// radian/s

// time conversion
const double g_dfTTmTAIInSec = 32.184;	
const double g_dfTTmTAIInJD = g_dfTTmTAIInSec / 86400.0;
const double g_dfTAImGPSInSec = 19.0;	

const int g_nMaxStations = 10;
const int g_nMaxDays = 10;
const int g_nComputeInterval = 5;		// computation step length in second of time
const int g_nMaxDarkPeriods = 2 * g_nMaxDays;
const int g_nMaxEpochs = 86400 / g_nComputeInterval * g_nMaxDays;

// reference JD and second used in EOSDOD
const double g_dfRefJDHere = 2450000.0;
const double g_dfRefSECHere = g_dfRefJDHere * 86400.0;

// numberic integrator
const int g_nMaxInteOrder = 25;
const int g_nMaxBackValue = 25;

// number of unknowns to be solved
const int g_nMaxUnknowns = 1000;
const int g_n3MaxUnknowns = 3 * g_nMaxUnknowns;

// geopotential model maximum order (degree)
const int g_nMaxGeoDegree = 180 + 1;
const int g_nMaxGeoTerms = g_nMaxGeoDegree * ( g_nMaxGeoDegree + 1 ) / 2 + g_nMaxGeoDegree * 3 + 10;

// integration interval of Gauss-Jackson integrator
const int g_nGJIntegrationInterval = 60;	

// maximum number of satellite
const int g_nMaxSatellites = 50000;

/*****************************************************************************

  for camera calibration
******************************************************************************/

const int g_nMaxNumStars = 50000;

#endif