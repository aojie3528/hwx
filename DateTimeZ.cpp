/************************************************************************************
 methods about the date and time

 Author: Jizhang Sang, EOS

 (C) COPYRIGHT -Electro Optic Systems Pty. Ltd. Australia (EOS) 2002,
 All rights reserved. No part of this program may be photocopied, 
 reproduced or translated to another program language without the prior
 written consent of EOS

*************************************************************************************/
// #include "stdafx.h"

#include <math.h>

#include "constant.h"

#include "DateTimeZ.h"

class cDateTimeZ g_DateTimeZ;

cDateTimeZ::cDateTimeZ()
{

}

cDateTimeZ::~cDateTimeZ()
{

}

/*************************************************************************************************

  Convert the given UTC Date/Time into Julian Date.

  nYear: 4-digit year
  nMonth: month of the year, 1 for January, ..., 12 for December
  nDay: day of the month, 1 for the first day of the month, 2 the seocnd day, ...
  nHour: hour of the day, 0 for the hour after the midnight, 1 for the hour after 1 hour, ...
  nMinute: minute of the hour, 0 for the minute after the hour, 1 for the minute after 1 minute, ...
  dfSecond: seocnd of the minute.
  
  dfJD: computed Julian Date.

**************************************************************************************************/

BOOL cDateTimeZ::DateTime2JD( int nYear, int nMonth, int nDay, 
							  int nHour, int nMinute, double dfSecond,
	                          double &dfJD )
{
//	if( nYear < 1900 || nYear > 2100 ) return FALSE;
//	if( nMonth < 1 || nMonth > 12 ) return FALSE;
//	if( nDay < 1 || nDay > 366 ) return FALSE;
//	if( nHour < 0 || nHour > 720 ) return FALSE;
//	if( nMinute < 0 || nMinute > 1440 ) return FALSE;
//	if( dfSecond < 0.0 || dfSecond > 86400.0 ) return FALSE;

	int nJD =  367 * nYear - 
	           7 * ( nYear + ( nMonth + 9 ) / 12 ) / 4 +
		       275 * nMonth / 9 + nDay + 1721013;
	double dfRealDay = nHour / 24.0 + nMinute / 1440.0 + 
	                   dfSecond / 86400.0 + 0.5;

	dfJD = (double) nJD + dfRealDay;

	return TRUE;

}

/*************************************************************************************************

  Convert the given UTC Date/Time into Julian Date.

  nYear: 4-digit year
  nMonth: month of the year, 1 for January, ..., 12 for December
  nDay: day of the month, 1 for the first day of the month, 2 the seocnd day, ...
  nHour: hour of the day, 0 for the hour after the midnight, 1 for the hour after 1 hour, ...
  nMinute: minute of the hour, 0 for the minute after the hour, 1 for the minute after 1 minute, ...
  dfSecond: seocnd of the minute.

  dfIntJD: the integer part of the computed Julian Date
  dfFractionJD: the fractional part of the computed Julian Date

**************************************************************************************************/

BOOL cDateTimeZ::DateTime2JD( int nYear, int nMonth, int nDay, 
							  int nHour, int nMinute, double dfSecond,
	                          double &dfIntJD, double &dfFractionJD )
{
	double dfJD;

	if( !DateTime2JD( nYear, nMonth, nDay, nHour, nMinute, dfSecond, dfJD ) ) return FALSE;

	m_dfIntJD = dfIntJD = (double) ( (int) dfJD);

	m_dfFractionJD = dfFractionJD = dfJD - dfIntJD;

	return TRUE;

}

/*************************************************************************************************

  Convert the Julian Date to the UTC Date/Time.

  dfJD: given Julian Date.

  nYear: 4-digit year
  nMonth: month of the year, 1 for January, ..., 12 for December
  nDay: day of the month, 1 for the first day of the month, 2 the seocnd day, ...
  nHour: hour of the day, 0 for the hour after the midnight, 1 for the hour after 1 hour, ...
  nMinute: minute of the hour, 0 for the minute after the hour, 1 for the minute after 1 minute, ...
  dfSecond: seocnd of the minute.

**************************************************************************************************/

BOOL cDateTimeZ::JD2DateTime( double dfJD,
	                          int &nYear, int &nMonth, int &nDay, 
							  int &nHour, int &nMinute, double &dfSecond )
{
	if( fabs( dfJD - g_dfJ2000 ) > 36525.0 ) 
	{
		//MessageBox( NULL, "Time out of year 1900-2100", 
		//	        "JD2DateTime in DateTimeZ", MB_OK | MB_ICONHAND );
		return FALSE;
	}
	
	double t1, t2, t3, t4, date_mjd, day_real;
  
    date_mjd = dfJD - g_dfMJDRef;

	int jd = ( int )dfJD;

	t4 = fmod( date_mjd, 1.0 );
	t1 = 1. + date_mjd - t4 + 2400000.0;
	int ih =  ( int )( ( t1 - 1867216.25 ) / 36524.0 );
	t2 = t1 + 1. + ih - ih / 4;
	t3 = t2 - 1720995.0;
	int ih1 = ( int ) ( ( t3 - 122.1 ) / 365.25 );

	nYear = ih1;
	t1 = 365.25 * ih1 - fmod( 365.25 * ih1, 1.0 );
	int ih2 = ( int )( ( t3 - t1 ) / 30.6001 );
	ih1 = ( int ) ( 30.6001 * ih2 );
	day_real = t3 - t1 - ih1 + t4;
	nMonth = ih2 - 1;
	if( ih2 > 13 ) nMonth = ih2 - 13;
	if( nMonth <= 2 ) nYear += 1;

	nDay = ( int ) ( day_real );
	
    t1 = ( day_real - nDay ) * 24.0;
	nHour = ( int ) ( t1 );

	t2= ( t1 - nHour ) * 60.0;
	nMinute = ( int ) ( t2 );

	dfSecond = ( t2 - nMinute ) * 60.0;

	if( fabs( dfSecond - 59.999 ) < 0.001 )
	{
		dfSecond = 0.0;
		nMinute += 1;
	}

	if( dfSecond < 0.0 ) dfSecond = 0.0;

	if( nMinute - 60 == 0 )
	{
		nMinute = 0;
		nHour += 1;
	}

	if( nHour == 24 )
	{
		nHour = 0;
		nDay++;
	}

	m_nYear = nYear;
	m_nMonth = nMonth;
	m_nDay = nDay;
	m_nHour = nHour;
	m_nMinute = nMinute;
	m_dfSecond = dfSecond;

	return TRUE;
}

/*************************************************************************************************

  Convert the Julian Date to the UTC Date/Time.

  dfIntJD: the integer part of the computed Julian Date
  dfFractionJD: the fractional part of the computed Julian Date

  nYear: 4-digit year
  nMonth: month of the year, 1 for January, ..., 12 for December
  nDay: day of the month, 1 for the first day of the month, 2 the seocnd day, ...
  nHour: hour of the day, 0 for the hour after the midnight, 1 for the hour after 1 hour, ...
  nMinute: minute of the hour, 0 for the minute after the hour, 1 for the minute after 1 minute, ...
  dfSecond: seocnd of the minute.

**************************************************************************************************/

BOOL cDateTimeZ::JD2DateTime( double dfIntJD, double dfFractionJD,
	                          int &nYear, int &nMonth, int &nDay, 
							  int &nHour, int &nMinute, double &dfSecond )
{
	double dfJD = dfIntJD + dfFractionJD;

	if( !JD2DateTime( dfJD, nYear, nMonth, nDay, nHour, nMinute, dfSecond ) ) return FALSE;

	return TRUE;
}

/*************************************************************************************************

  Conversion between JD and MJD

  dfIntJD: the integer part of the computed Julian Date
  dfFractionJD: the fractional part of the computed Julian Date

  dfIntMJD: the integer part of the computed Modified Julian Date
  dfFractionMJD: the fractional part of the computed Modified Julian Date

  bJD2MJD: TRUE for converting from JD to MJD, FALSE from MJD to JD

**************************************************************************************************/

BOOL cDateTimeZ::ConversionJDMJD( double &dfIntJD,  double &dfFractionJD, 
		                          double &dfIntMJD, double &dfFractionMJD, BOOL bJD2MJD )
{
	if( bJD2MJD )	// from JD to MJD
	{
		m_dfIntJD = dfIntJD;
		m_dfFractionJD = dfFractionJD;

		dfIntMJD = dfIntJD - g_dfMJDRef;
		dfFractionMJD = dfFractionJD;

		ReConstruct( dfIntMJD, dfFractionMJD );

		m_dfIntMJD = dfIntMJD;
		m_dfFractionMJD = dfFractionMJD;
	}
	else
	{
		m_dfIntMJD = dfIntMJD;
		m_dfFractionMJD = dfFractionMJD;

		dfIntJD = dfIntMJD + g_dfMJDRef;
		dfFractionJD = dfFractionMJD;

		ReConstruct( dfIntJD, dfFractionJD );

		m_dfIntJD = dfIntJD;
		m_dfFractionJD = dfFractionJD;
	}

	return TRUE;
}

/*************************************************************************************************
  
  Conversion between JD and Second of Time since a predifined reference epoch (g_dfRefJDHere)

  dfIntJD: the integer part of the computed Julian Date
  dfFractionJD: the fractional part of the computed Julian Date

  dfIntSec: the integer part of the seconds since the reference epoch
  dfFractionMJD: the fractional part of the seconds since the reference epoch

  bJD2MJD: TRUE for converting from JD to Second of Time, FALSE from Second of Time to JD

**************************************************************************************************/

BOOL cDateTimeZ::ConversionJDSecond( double &dfIntJD,  double &dfFractionJD, 
		                             double &dfIntSec, double &dfFractionSec, BOOL bJD2Second )
{
	if( bJD2Second ) // from JD to Second
	{
		dfIntSec = ( dfIntJD - g_dfRefJDHere ) * 86400.0;
		dfFractionSec = dfFractionJD * 86400.0;

		ReConstruct( dfIntSec, dfFractionSec );

		m_dfIntSec = dfIntSec;
		m_dfFractionSec = dfFractionSec;
	}
	else
	{
		int nInt = (int) ( dfIntSec / 86400.0 );

		dfIntJD = g_dfRefJDHere + (double) nInt;

		dfFractionJD = ( fmod( dfIntSec, 86400.0 ) + dfFractionSec ) / 86400.0;

		ReConstruct( dfIntJD, dfFractionJD );

		m_dfIntJD = dfIntJD;
		m_dfFractionJD = dfFractionJD;
	}

	return TRUE;
}

/*************************************************************************************************
  
  Conversion between JD and Second of Time since a predifined reference epoch (g_dfRefJDHere)

  dfJD: the Julian Date
  
  dfIntSec: the integer part of the seconds since the reference epoch
  dfFractionMJD: the fractional part of the seconds since the reference epoch

  bJD2MJD: TRUE for converting from JD to Second of Time, FALSE from Second of Time to JD

**************************************************************************************************/

BOOL cDateTimeZ::ConversionJDSecond( double &dfJD, 
									 double &dfIntSec, double &dfFractionSec, BOOL bJD2Second )
{
	double dfIntJD, dfFractionJD;
	if( bJD2Second ) // from JD to Second
	{
		dfIntJD = (double) ( (int) dfJD );
		dfFractionJD = dfJD - dfIntJD;

		dfIntSec = ( dfIntJD - g_dfRefJDHere ) * 86400.0;
		dfFractionSec = dfFractionJD * 86400.0;

		ReConstruct( dfIntSec, dfFractionSec );

		m_dfIntSec = dfIntSec;
		m_dfFractionSec = dfFractionSec;
	}
	else
	{
		int nInt = (int) ( dfIntSec / 86400.0 );

		dfIntJD = g_dfRefJDHere + (double) nInt;

		dfFractionJD = ( fmod( dfIntSec, 86400.0 ) + dfFractionSec ) / 86400.0;

		ReConstruct( dfIntJD, dfFractionJD );

		m_dfIntJD = dfIntJD;
		m_dfFractionJD = dfFractionJD;

		dfJD = dfIntJD + dfFractionJD;
	}

	return TRUE;
}


/*************************************************************************************************

  Get the MJD

**************************************************************************************************/

void cDateTimeZ::GetMJD( double &dfIntMJD, double &dfFractionMJD )
{
	dfIntMJD = m_dfIntMJD;
	dfFractionMJD = m_dfFractionMJD;
}

/*************************************************************************************************

  Get the Second of Time

**************************************************************************************************/

void cDateTimeZ::GetSec( double &dfIntSec, double &dfFractionSec )
{
	dfIntSec = m_dfIntSec;
	dfFractionSec = m_dfFractionSec;
}


/*************************************************************************************************

  Redetermine the integer and fractional parts of a value

**************************************************************************************************/

void cDateTimeZ::ReConstruct( double &dfInt, double &dfFraction )
{
	while( dfFraction < 0.0 )
	{
		dfFraction += 1.0;
		dfInt -= 1.0;
	}

	int nInt = (int) dfInt;
	double dfR = dfInt - (double) nInt;
	
	dfInt = dfInt - dfR;

	dfFraction += dfR;

	nInt = (int) dfFraction;
	dfFraction = dfFraction - (double) nInt;

	dfInt += (double) nInt;
}

/*****************************************************************************

  Compute day of year from Date

*****************************************************************************/

BOOL cDateTimeZ::ComputeDayOfYearFromDate( int nYear, int nMonth, int nDay, int &nDoY )
{
	if( abs( nYear - 2000 ) > 100 ) return FALSE;

	if( nMonth < 1 || nMonth > 12 ) return FALSE;
	
	if( nDay < 1 || nDay > 366 ) return FALSE;

//	day in january
	nDoY = nDay;
	if( nMonth <= 1 ) return TRUE;

//  day in february
  
    nDoY += 31;
	if( nMonth <= 2 ) return TRUE;

//  day in march
	if( nYear / 4 * 4 == nYear ) nDoY += 29;
    else nDoY += 28;
	if( nMonth <= 3 ) return TRUE;

//  day in april
    nDoY += 31;
	if( nMonth <= 4 ) return TRUE;

//  day in may
    nDoY += 30;
	if( nMonth <= 5 ) return TRUE;

//  day in june
    nDoY += 31;
	if( nMonth <= 6 ) return TRUE;

//  day in July
    nDoY += 30;
	if( nMonth <= 7 ) return TRUE;

//  day in august
    nDoY += 31;
	if( nMonth <= 8 ) return TRUE;

//  day in september
    nDoY += 31;
	if( nMonth <= 9 ) return TRUE;

//  day in october
    nDoY += 30;
	if( nMonth <= 10 ) return TRUE;

//  day in november
    nDoY += 31;
	if( nMonth <= 11 ) return TRUE;

//  day in december
	nDoY += 30;
	
	return TRUE;
}

/*************************************************************************************************


**************************************************************************************************/

BOOL cDateTimeZ::ConversionJDMJD( BOOL bJD2MJD )
{
	if( bJD2MJD )	// from JD to MJD
	{
		m_dfIntMJD = m_dfIntJD - g_dfMJDRef;
		m_dfFractionMJD = m_dfFractionJD;

		ReConstruct( m_dfIntMJD, m_dfFractionMJD );
	}
	else
	{
		m_dfIntJD = m_dfIntMJD + g_dfMJDRef;
		m_dfFractionJD = m_dfFractionMJD;

		ReConstruct( m_dfIntJD, m_dfFractionJD );
	}

	return TRUE;
}

/*************************************************************************************************


**************************************************************************************************/

BOOL cDateTimeZ::ConversionJDSecond( BOOL bJD2Second )
{
	if( bJD2Second ) // from JD to Second
	{
		m_dfIntSec = ( m_dfIntJD - g_dfRefJDHere ) * 86400.0;
		m_dfFractionSec = m_dfFractionJD * 86400.0;

		ReConstruct( m_dfIntSec, m_dfFractionSec );
	}
	else
	{
		int nInt = (int) ( m_dfIntSec / 86400.0 );

		m_dfIntJD = g_dfRefJDHere + (double) nInt;

		m_dfFractionJD = ( fmod( m_dfIntSec, 86400.0 ) + m_dfFractionSec ) / 86400.0;

		ReConstruct( m_dfIntJD, m_dfFractionJD );
	}

	return TRUE;
}

/*************************************************************************************************

  Compute the local hour of a place given the Julian Date and the longitude of the place

  dfLongitude: the longitude of the place of interest, in radian
  dfLocalHour: the computed local hour in 24h system

  nDoY: day of year of the given dfJD

**************************************************************************************************/

BOOL cDateTimeZ::ComputeLocalHour( double dfJD, double dfLongitude, double &dfLocalHour, int &nDoY )
{
	int nYear, nMonth, nDay, nHour, nMinute;
	double dfSecond;

	if( !JD2DateTime( dfJD, nYear, nMonth, nDay, nHour, nMinute, dfSecond ) ) return FALSE;

	dfLocalHour = (double) nHour * 3600.0 + (double) nMinute * 60.0 + dfSecond;

	dfLocalHour /= 3600.0;

	dfLocalHour += dfLongitude * g_dfRAD2DEG / 15.0;

	dfLocalHour = fmod( dfLocalHour, 24.0 );

	ComputeDayOfYearFromDate( nYear, nMonth, nDay, nDoY );

	return TRUE;
}

/*************************************************************************************************

  Compute the day of year of thegiven JD

**************************************************************************************************/

BOOL cDateTimeZ::ComputeDayOfYearFromJD( double dfJD, int &nDoY )
{
	int nYear, nMonth, nDay, nHour, nMinute;
	double dfSecond;

	if( !JD2DateTime( dfJD, nYear, nMonth, nDay, nHour, nMinute, dfSecond ) ) return FALSE;

	ComputeDayOfYearFromDate( nYear, nMonth, nDay, nDoY );

	return TRUE;
}

/*************************************************************************************************

  Convert between the bessilan year and julian date

  bB2J: TRUE conversion from besselian year to Julian date, 
        FALSE conversion from Julian date to Besselian year 

**************************************************************************************************/
BOOL cDateTimeZ::ConversionBesselianYearJulianDate( double &dfBesselianYear, double &dfJD, BOOL bB2J )
{
	double dfD1900 = 1900.0;
	double dfBYR00 = 2415020.31352;	// Julian date of Besselian Yaer 1900.0
	double dfDays = 365.24219879;	// days of a besselian yaer
	double dfDDays = -8.56e-09;		// rate of length of a day

	if( bB2J )
	{
		double D = dfBesselianYear - dfD1900;

		dfJD = D * ( dfDays + dfDDays * D );

		dfJD += dfBYR00;
	}
	else	// from Julian date to Besselian Date
	{
		double D = dfJD - dfBYR00;

		dfBesselianYear = D / dfDays;

		dfBesselianYear = D / ( dfDays + dfDDays * dfBesselianYear ) + dfD1900;
	}

	return TRUE;
}

/******************************************************************************

  Compute Equation of Time according to a formular on

  http://holodeck.st.usm.edu/vrcomputing/vrc_t/tutorials/solar/eot.shtml

******************************************************************************/

double cDateTimeZ::ComputeEquationOfTime( int nDayOfYear )
{
	double b = 360.0 * ( (double)nDayOfYear - 81.0 ) / 364.0 * g_dfDEG2RAD;

	double eot = 9.87 * sin( 2.0 * b ) - 7.53 * cos( b ) - 1.5 * sin( b );	// in minutes

	eot = eot / 60.0 * g_dfHOUR2RAD;	// in radians

	return eot;
}

/******************************************************************************

  Compute Equation of Time according to a formular on

  http://holodeck.st.usm.edu/vrcomputing/vrc_t/tutorials/solar/eot.shtml

******************************************************************************/
double cDateTimeZ::ComputeEquationOfTime( int nYear, int nMonth, int nDay )
{
	int nDayOfYear;

	ComputeDayOfYearFromDate( nYear, nMonth, nDay, nDayOfYear );

	return ComputeEquationOfTime( nDayOfYear );
}

/******************************************************************************

  Compute Equation of Time according to a formular on

  http://holodeck.st.usm.edu/vrcomputing/vrc_t/tutorials/solar/eot.shtml

******************************************************************************/
double cDateTimeZ::ComputeEquationOfTime( double dfJD )
{
	int nDayOfYear;

	ComputeDayOfYearFromJD( dfJD, nDayOfYear );

	return ComputeEquationOfTime( nDayOfYear );
}

/******************************************************************************

  Compute the local solar time according to steps below

  1. given UT, then, the local mean solar time, MT, is

     LMST = UT + Local Longitude

  2. local solar time, LST, is 

     LST = LMST + Equation of Time

  Note: nYear, nMonth, nDay, nHour, nMinute, dfSecond should be UT time.

        dfLocalLongitude: the geocentric (geodetic) longitude of the location
		                  of interest
******************************************************************************/

double cDateTimeZ::ComputeLocalSolarTime( int nYear, int nMonth, int nDay, 
		                                  int nHour, int nMinute, double dfSecond,
								          double dfLocalLongitude )
{
	double UT = ( (double)nHour + (double)nMinute / 60.0 + dfSecond / 3600.0 ) * 
		        g_dfHOUR2RAD;

	double LMST = UT + dfLocalLongitude;

	double LST = LMST + ComputeEquationOfTime( nYear, nMonth, nDay );

	if( LST < 0.0 ) LST += g_dfTWOPI;

	LST = fmod( LST, g_dfTWOPI );

	return LST;
}

/******************************************************************************

  Compute the local solar time according to steps below

  1. given UT, then, the local mean solar time, MT, is

     LMST = UT + Local Longitude

  2. local solar time, LST, is 

     LST = LMST + Equation of Time

  Note: dfJD should be UT time.

        dfLocalLongitude: the geocentric (geodetic) longitude of the location
		                  of interest
******************************************************************************/

double cDateTimeZ::ComputeLocalSolarTime( double dfJD, double dfLocalLongitude )
{
	int nYear, nMonth, nDay, nHour, nMinute;
	double dfSecond;

	JD2DateTime( dfJD, nYear, nMonth, nDay, nHour, nMinute, dfSecond );

	return ComputeLocalSolarTime( nYear, nMonth, nDay, nHour, nMinute, dfSecond,
								  dfLocalLongitude );
}
