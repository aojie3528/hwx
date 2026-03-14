/********************************************************************************** 

 Data Structure Defined in This Software

 Author: Jizhang Sang, Feb-Apr, 2002

 (C) COPYRIGHT -Electro Optic Systems Australia (EOS) 2002,
 All rights reserved. No part of this program may be photocopied, 
 reproduced or translated to another program language without the prior
 written consent of EOS

 Reference: IAU SOFA EPV00 for Model 1
  
**********************************************************************************/

#ifndef INC_DATASTRUCTURE
#define INC_DATASTRUCTURE

#include <string>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <deque>
#include <map>
#include <set>
#include <cstring>

#if (_MSVC_LANG > 201402L) || (__cplusplus >= 201703L) // C++ 17
#include <filesystem>
namespace fs = std::filesystem;
#else
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem::v1;
#endif

using namespace std;

/***********************************************************************************

 information about a tracking station

***********************************************************************************/
struct stTrackStation
{
	double dfLatitude, dfLongitude, dfHeight;	// Geodetic position;
	double dfSinLat, dfCosLat;					
	double dfECEFX, dfECEFY, dfECEFZ;			// in meter

	char szName[ FILENAME_MAX   ];
	int nSiteNo;

	// 1 for Laser 
	// 2 for Optical
	// 3 for both
	int nTrackingFacility;
};

/***********************************************************************************

 information about status of a tracking station

***********************************************************************************/
struct stStationStatus
{
	int nStationID;
	bool bAvailable;
	double dfJDBegin, dfJDEnd;

	// unavailability code, for example, 
	// 1 for weather, 
	// 2 for system problem, 
	// 3 for maintance, 
	// etc
	int nCode;	

	// 1 for Laser 
	// 2 for Optical
	// 3 for both
	int nTrackingFacility;
};

/***********************************************************************************

 Parameters thats defines the visibility of a satellite

***********************************************************************************/

struct stVisibilityCondition
{
	bool bSun;				// true if the sun illumination is required

	double dfElevationMask;	// in radian
	double dfZenithMask;	// 90 degree - dfElevationMask

	double dfSunElevationMask;	// in radian

	double dfSunReflectAngle;	// in radian
};


/***********************************************************************************

 A data point of a satellite pass

***********************************************************************************/

struct stPassPoint
{
	double dfJD, dfAz, dfEl;
};

/***********************************************************************************

 Parameters that defines a satellite pass

***********************************************************************************/
struct stVisiblePass
{
	int nSatID;
	int nStationID;

	// 1 for Laser 
	// 2 for Optical
	// 3 for both
	int nTrackingFacility;

	stPassPoint tRise, tSet, tTCA;

	// sun lit data
	bool bSunLit;

	stPassPoint tSunLitRise, tSunLitSet, tSunLitMid;

	double dfBenifit;
};

/********************************************************************************

  Data that defines an epochPosVel

*********************************************************************************/
struct epochPosVel1
{
	double JDUTC = 0.0;
	double ele[6] = {0.0};
	double pos[3] = { 0.0 }, vel[3] = {0.0};
	double Az = 0.0, El = 0.0, DisTwoway = 0.0;
	double Ra = 0.0, Dec = 0.0;
	double Mag = 0.0, RCS = 0.0;

	double Weight = 1.0;

	int SiteIDs = 0;
	int SiteIDr = 0;
	double satpos[3] = {0.0}; //site pos in TOD, in meter

	bool operator==(const epochPosVel1& tPV)const
	{
		const double EPS = 1e-10;
		if (fabs(this->JDUTC - tPV.JDUTC) < EPS)  //判断是否相等
			return true;
		return false;
	}

	bool operator<(const epochPosVel1& tPV)const
	{
		if (this->JDUTC < tPV.JDUTC)
			return true;
		return false;
	}

};

/********************************************************************************

  Data that defines an epochPosVel

*********************************************************************************/
struct IODInfo
{
	fs::path m_ObsPath;
	double mag = 0.0;
	double duration = 0.0;
	int numobs = 0;
	int m_Type = 0; //Type 0:LD  1:Ground Optical 2:Space Optical
	epochPosVel1 m_IOD;
};


/********************************************************************************

  Data that defines an TCA of Site and Obj

*********************************************************************************/
struct stTCA
{
	int m_SiteIDs = 0;
	int m_SiteIDr = 0;

	int m_ObjID = 0;
	string m_ArcName;

	epochPosVel1 m_Orbit;

	double JDUTC = 0.0;
	double Dis = 0.0;
};


struct stSiteInfo
{
	int m_IDr = 0;
	double m_siteBLHr[3] = {0.0};
	double m_siteECEFr[3] = {0.0};

	int m_IDs = 0;
	double m_siteBLHs[3] = {0.0};
	double m_siteECEFs[3] = {0.0};

	int m_StationType = 1; //0 地基光学,  1 分布式 Two way, 2 相控阵 One way, 3 天基光学

	//For SBO
	double opticalAxisAzEl[2]; // [0] for az, [1] for el
	int FOVType = 1;			// 1 for circular - centered at the optical axis
								// 2 for square - centered at the optical axis
								// 3 for rectangular - centered at the optical axis
								// 4 user-defined
	double FOVSize[10];		// if FOVType == 1 or FOVType == 2, FOVSize[ 0 ] is the half angle, in rad
							// if FOVType == 3, FOVSize[ 0 ] is the first half angle, FOVSize[ 1 ] is the second half angle
							// if FOVType == 4, FOVSize[ 0 ] and FOVSize[ 1 ] are the azimuth range,
							//                  FOVSize[ 2 ] and FOVSize[ 3 ] are the elevation range
	double FOVAzElRange[4];
	double m_SunCloseAngle;
	double m_MoonCloseAngle;
	double m_phaseAngle;
	double m_VisualMagnitude;

};


/********************************************************************************

  Data that defines an Radar Arc

*********************************************************************************/
struct RadarArc
{
	int m_Type = 0; //Type 0:LD  1:Ground Optical 2:Space Optical
	double m_weight = 1.0;

	bool isNew = false;
	bool isOut = false;

	int m_ID = 0;
	char m_FileName[255];

	int m_numObs = 0;
	double m_RCS = 0.0;
	double m_Mag = 0.0;
	double m_StartJD = 0.0;
	double m_EndJD = 0.0;
	double m_Duration = 0.0;		//	in Sec
	bool isIO = false;
	bool isTwoArcIO = false;
	epochPosVel1 m_InitialOrbit, m_TwoArcOrbit, m_PreciseOrbit;

	int m_SiteIDs = 0;
	int m_SiteIDr = 0;
	stSiteInfo m_Site;

	std::vector<epochPosVel1>* m_pObsPos;
	std::vector<epochPosVel1> m_ObsPos;
	std::vector<epochPosVel1> m_IODEPHOrbit; // IOD轨道传播的结果，breakup

	bool isTLE = false;
	double m_TLEEle[10];

	std::string m_PropertyStr;

	bool m_isMergeArc = false;
	vector<string> m_OriginalArcs;

	bool operator==(const RadarArc& tArc)const
	{
		if (strcmp(this->m_FileName, tArc.m_FileName))
			return false;
		return true;
	}
	bool operator<(const RadarArc& tArc)const
	{
		if (strcmp(this->m_FileName, tArc.m_FileName) < 0)
			return true;
		return false;
	}

};

struct stDayObs
{
	int m_Date;
	vector<string> m_Obs;
};

struct stSiteEpoch
{
	int SiteID = 0;
	double JDUTC = 0.0;

	//bool operator==(const stSiteEpoch& tEpoch)const
	//{
	//	if (this->SiteID != tEpoch.SiteID)
	//		return false;
	//	if (fabs(this->JDUTC - tEpoch.JDUTC) * 86400.0 > 2.0)
	//		return false;
	//	return true;
	//}

	bool operator<(const stSiteEpoch& tEpoch)const
	{
		if (this->JDUTC < tEpoch.JDUTC)
			return true;
		//if (this->SiteID < tEpoch.SiteID)
		//	return true;
		return false;
	}

};

/********************************************************************************

  Arc Info for Group division

*********************************************************************************/
struct stArcAcross
{
	string ArcName;
	double m_Inc = 0.0;
	double m_RAAN = 0.0;

	bool operator==(const stArcAcross& tP)const
	{
		if (this->ArcName != "" && this->ArcName == tP.ArcName)
			return true;
		if (this->m_Inc == tP.m_Inc && this->m_RAAN == tP.m_RAAN)
			return true;
		return false;
	}
	bool operator<(const stArcAcross& tP)const
	{
		if (this->ArcName < tP.ArcName)
			return true;
		return false;
	}
};


/********************************************************************************

  Arc Group Information

*********************************************************************************/
struct stArcGroup
{
	double m_Inc = 0.0;
	double m_RAAN = 0.0;
	int m_GroupID = 0;
	vector<string> m_Arcs;
};


/********************************************************************************

  Data that defines an TLE or Ele

*********************************************************************************/
struct stIndependentArc
{
	std::vector<RadarArc> m_Arc;
	std::vector<RadarArc*> m_pArc;
};


/********************************************************************************

  Data that defines an TLE or Ele

*********************************************************************************/
struct stTLEEle
{
	double m_TLETime = 0.0;
	vector<double> m_TLEEle;
	bool m_isOut = false;

	//pdfTLE[2] TLETime
	bool operator==(const stTLEEle& tT)const
	{
		if (fabs(this->m_TLEEle[2] - tT.m_TLEEle[2]) < 1.0e-5)
			return true;
		return false;
	}
	bool operator<(const stTLEEle& tT)const
	{
		if (this->m_TLEEle[2] < tT.m_TLEEle[2])
			return true;
		return false;
	}
};

/***********************************************************************************

 Orbital elements that define a satellite orbit

***********************************************************************************/

struct stSatelliteIOE
{
	int nSatelliteID;
	int nSIC;
	char cElementType;			// T for TLE, E for EOS, I for IRV
	int nIntJD, nFractionJD;

	int pnElement1to6[6];
	float pfElement7to18[20];
	bool bNewElement;
	double dfTimeLastTracked;
	bool bTrackingtarget;
	bool bAcquired;
	bool bTrackable;
	bool bTrackingTested;
	bool bSecondTrack;

	bool bLEO;
	bool bNORAD;

	double Diameter;//m 20180509

	std::vector< std::string > trackbenefit;

	void Assign(const stSatelliteIOE stSource, stSatelliteIOE& stResult) const
	{
		stResult.cElementType = stSource.cElementType;
		stResult.nFractionJD = stSource.nFractionJD;
		stResult.nIntJD = stSource.nIntJD;
		stResult.nSatelliteID = stSource.nSatelliteID;
		stResult.nSIC = stSource.nSIC;
		stResult.bNewElement = stSource.bNewElement;
		stResult.dfTimeLastTracked = stSource.dfTimeLastTracked;
		stResult.bTrackingtarget = stSource.bTrackingtarget;
		stResult.bAcquired = stSource.bAcquired;
		stResult.bTrackable = stSource.bTrackable;
		stResult.bTrackingTested = stSource.bTrackingTested;
		stResult.bSecondTrack = stSource.bSecondTrack;

		for (int i = 0; i < 6; i++) stResult.pnElement1to6[i] = stSource.pnElement1to6[i];
		for (int i = 0; i < 20; i++) stResult.pfElement7to18[i] = stSource.pfElement7to18[i];
		stResult.trackbenefit.clear();
		stResult.trackbenefit.reserve(stSource.trackbenefit.size());
		for (int i2 = 0; i2 < (int)stSource.trackbenefit.size(); i2++)
		{
			stResult.trackbenefit.push_back(stSource.trackbenefit[i2]);
		}
	}

	double GetRefJD() const
	{
		return ((double)nIntJD + (double)nFractionJD * 1.0e-9);
	};

	double GetBStar()
	{
		return pfElement7to18[2];
	}

	double Get_nDot()
	{
		return pfElement7to18[0];
	}

	double GetTimeLastTracked() const
	{
		// may need to scope lock, just incase the data changes
		if (dfTimeLastTracked < 0)
		{
			return GetRefJD();
		}
		else
		{
			return dfTimeLastTracked;
		}
	};
	void SetTimeLastTracked(double dfTime)
	{
		dfTimeLastTracked = dfTime;
	}
	void SetRefJD(double dfJD)
	{
		nIntJD = (int)dfJD;
		nFractionJD = (int)((dfJD - nIntJD) * 1.0e9);
	}

	void GetPosVel(double* pdfPos, double* pdfVel) const
	{
		for (int i = 0; i < 3; i++)
		{
			pdfPos[i] = (double)pnElement1to6[i] * 1.e-1;
			pdfVel[i] = (double)pnElement1to6[i + 3] * 1.0e-4;
		}
	}

	double GetOrbitalPeriod() const
	{
		return (double)pfElement7to18[11];
	};

	double GetInclination() { return (double)pnElement1to6[2] * 1.0e-5; };
	double GetEccentricity() { return (double)pnElement1to6[1] * 1.0e-7; };
	double GetRANode() { return (double)pnElement1to6[3] * 1.0e-5; };
	double GetPerigeeAugument() { return (double)pnElement1to6[4] * 1.0e-5; };
	double GetMeanAmoaly() { return (double)pnElement1to6[5] * 1.0e-5; };

	// added by RO december 2004 to store the contents of each line of the TLE
	std::string Line1, Line2, Line3;

	stSatelliteIOE& operator=(const stSatelliteIOE& RHS)
	{
		if (this == &RHS) return *this;
		Assign(RHS);
		return *this;
	};

	void Assign(const stSatelliteIOE& RHS)
	{
		nSatelliteID = RHS.nSatelliteID;
		nSIC = RHS.nSIC;
		cElementType = RHS.cElementType;			// T for TLE, E for EOS, I for IRV
		nIntJD = RHS.nIntJD;
		nFractionJD = RHS.nFractionJD;
		bNewElement = RHS.bNewElement;
		dfTimeLastTracked = RHS.dfTimeLastTracked;
		bTrackingtarget = RHS.bTrackingtarget;
		bAcquired = RHS.bAcquired;
		bTrackable = RHS.bTrackable;
		bTrackingTested = RHS.bTrackingTested;
		bSecondTrack = RHS.bSecondTrack;

		int i;
		for (i = 0; i < 6; i++)
		{
			pnElement1to6[i] = RHS.pnElement1to6[i];
		}
		for (i = 0; i < 12; i++)
		{
			pfElement7to18[i] = RHS.pfElement7to18[i];
		}
		Line1 = RHS.Line1;
		Line2 = RHS.Line2;
		Line3 = RHS.Line3;

		trackbenefit.clear();
		trackbenefit.reserve(RHS.trackbenefit.size());
		for (int i2 = 0; i2 < (int)RHS.trackbenefit.size(); i2++)
		{
			trackbenefit.push_back(RHS.trackbenefit[i2]);
		}

	}
	stSatelliteIOE(const stSatelliteIOE& RHS)
	{
		Assign(RHS);
	}
	stSatelliteIOE()
	{
		nSatelliteID = 0;
		nSIC = 0;
		cElementType = 'T';			// T for TLE, E for EOS, I for IRV
		nIntJD = 0;
		nFractionJD = 0;
		int i;
		for (i = 0; i < 6; i++)
		{
			pnElement1to6[i] = 0;
		}
		for (i = 0; i < 12; i++)
		{
			pfElement7to18[i] = 0;
		}
		Line1 = "";
		Line2 = "";
		Line3 = "";
		bNewElement = true;
		dfTimeLastTracked = -10;
		bTrackingtarget = false;
		trackbenefit.clear();
		trackbenefit.reserve(1000);
		bAcquired = false;
		bTrackable = false;
		bTrackingTested = false;
		bSecondTrack = false;
	}

	double GetSemiMajor()
	{
		double dfn = 2 * 3.1415926 * pnElement1to6[0] * 1e-8 / 86400.0;
		return pow(3.986004418e14 / (dfn * dfn), 1.0 / 3.0);
	}

	double GetPerigeeHeight()	// in m
	{
		return (GetSemiMajor() * (1.0 - GetEccentricity()));
	}

	double GetAltitudeKM()
	{
		return (GetSemiMajor() - 6378137.0) / 1000.0;
	}

};


/********************************************************************************

  struct for Obs Matching

*********************************************************************************/
struct visPass1
{
	int noradID;
	double passStartJD, passEndJD, passTCAJD;
	double passDurationSecond;
	double passStartAz, passEndAz;
	double passStartEl, passEndEl, tcaEl;
	double tcaRange, maxRange;
	double dt;
	stSatelliteIOE IOE;

	bool isEph = false;
	map<int, vector<epochPosVel1>>::iterator m_Ephiter;

	std::vector<double> JD, az, el, Dis, TODRA, TODDec;
	int numRAD;

	bool isMatched = false;
	int SiteIDr = 0;
};

/********************************************************************************

  struct for Eph Info

*********************************************************************************/
struct stEphInfo
{
	fs::path m_ephPath;
	int coorSystem = 1; //0 TOD, 1 J2000, 2 ECEF
	int timeSystem = 1; //0 beijing, 1 UTC
};

/********************************************************************************

  struct for Obs Matching

*********************************************************************************/
struct matchedPass1
{
	char fileName[200] = "";
	int noradID = 0;
	int siteID = 0;
	int numobs;
	double duration;
	double dt = 0.0;
	double RCS = 0.0;
	double Mag = 0.0;
	double Dis = 0.0;

	double maxRADif, maxDecDif, maxDisDif;
	double medianRADif = 0.0, medianDecDif = 0.0, medianDisDif = 0.0, medPosDif;
	double RArateDif = 0.0, DecrateDif = 0.0;
	double dSec = 0.0;
	double RMS = 0.0;  //CLY

	double confidence = 0.0;  // 置信度

	double maxAZDif, maxELDif;
	double medianAZDif, medianELDif;

	double StartJD = 0.0;
	double EndJD = 0.0;

	RadarArc IODELE;

	bool operator<(const matchedPass1& stMP)const
	{
		if (strcmp(this->fileName, stMP.fileName) < 0)
			return true;
		return false;
	}

};

/********************************************************************************

  struct for Obs Matching

*********************************************************************************/
struct obsPass1
{
	char fileName[200];
	int noradID;
	int numobs;
	int equipIDs, equipIDr;
	double maxd = 0.0, mind = 40000000.0;
	std::vector<double> JD, az, el, dis, ra, dec, TODRA, TODDec;
	std::vector<double> rcs_all;

	//for 天基光学
	vector<double> satx;
	vector<double> saty;
	vector<double> satz;
	vector<double> mag_all;

	int m_Type = 2; //观测值类型 0:LD   1:Ground based Optical  2:Space based Optical
	//
	int numtemp;
	double duration;

	double RaRate = 0.0;
	double DecRate = 0.0;

	double RCS = 0.0;
	double Mag = 0.0;
};

/********************************************************************************

  Data that defines an object with Arcs and TLEs

*********************************************************************************/
struct AssoObj
{
	bool isID = false;
	int m_ID = 0;
	std::deque<RadarArc> m_AssoArc;
	std::deque<stIndependentArc> m_CatArc;
	double m_TLEEle[20] = { 0.0 };
	bool isNew = false; 

	vector<stTLEEle> m_AllTLE;

	bool isPreciseTLE = false;
	bool isTLE = false;
	bool isEph = false;
	char m_EphFile[255];
	map<int, vector<epochPosVel1>>::iterator m_Ephiter;

	bool isMultiCD = false;
	double m_CD = 2.2;
	double m_CR = 1.1;
	double m_A2M = 0.02;
	double m_Mag = 0.0;
	bool isImportant = false;
	double m_ThrustAcc[3] = { 0.0 };
	double m_RMS = 0.0;

	bool isMatchedObs = false;
	bool isMineSuccess = false;
};

struct preAssoPair
{
	int m_GroupID = 0;
	RadarArc* m_pArc1 = NULL;
	RadarArc* m_pArc2 = NULL;
	char m_FileName1[255] = "";
	char m_FileName2[255] = "";
};

struct AssoArcSet
{
	int m_Date;
	std::vector<RadarArc> m_ArcLib;
};

/********************************************************************************

  Data Base of Objects

*********************************************************************************/
struct ObjectLibrary
{
	int nIDNumber = 0;  //已编号目标数
	std::vector<AssoObj> m_Objs;

	map<int, vector<epochPosVel1>> m_EphMap;
};


struct stTLEManeuver
{
	double StartJD = 0.0;
	double EndJD = 0.0;
	double StartSMA = 0.0;
	double EndSMA = 0.0;
	stSatelliteIOE startIOE;
	stSatelliteIOE endIOE;

	double m_dv[3];
	double dSMA = 0.0; //m
	double dInc = 0.0; //deg
};

/********************************************************************************

  Data that defines an object with Arcs and TLEs

*********************************************************************************/
struct stManeuverObj
{
	int m_ID = 0;
	std::vector<RadarArc> m_ManeuverArc;

	double m_TLEEle[20] = { 0.0 };
	vector<stTLEEle> m_AllTLE;

	int m_Type = 0; // 0:Impluse thrust  1: Continuous thrust
	double m_dv[3] = { 0.0 };

	double ManeuverJD = 0.0;
	epochPosVel1 preOrbit, postOrbit;

	bool isManeuver = true;
};

struct stBreakupEvent
{
	double JDUTC;
	int NORADID;//母体
	int WHUID;  //WHU编目库ID

	int DebrisNum;
	stSatelliteIOE m_IOE;
	std::vector<matchedPass1> Arcs;//事件涉及到的观测弧段
	std::vector<RadarArc> RadarArcs;
};

struct stCovisionArc
{
	set<string> m_Arcs;
	double CoTime = 0;

	RadarArc m_Arc1;
	RadarArc m_Arc2;
	RadarArc m_CommonArc;
};

struct stObjectAttribute
{
	int ID = 0;
	double CD = 2.2;
	double CR = 1.1;
	double A2M = 0.02;
	double Mag = 0.0;
	bool isImportant = false;
};

/***********************************************************************************

 Parameters that define the dark time perods over a time period

***********************************************************************************/

struct stDarkTime
{
	int nStationID;
	int nNumberDarkPeriod;
	double *pdfBeginJD, *pdfEndJD;

	stDarkTime()
	{
		pdfBeginJD = new double[ 500 ];
		pdfEndJD = new double[ 500 ];
	}

	~stDarkTime()
	{
		delete []pdfBeginJD;
		delete []pdfEndJD;
	}
};

/***********************************************************************************

 Parameters that control the computation of satellite veisibility

***********************************************************************************/

struct stVisComControl
{
	bool bInit;

	char szStationDataFile[ FILENAME_MAX ];
	char szSatelliteDataFile[ FILENAME_MAX ];

	// start date/time and prediction period of time
	int nYear, nMonth, nDay, nHour;	// in UTC
	double dfPredictionDays;
	double dfJDBegin;	// from Year/Month/Day/Hour	in UTC
	double dfJDEnd;		// dfJDBegin + dfPredictionDas	in UTC
	
	// gravity model
	char szGravityFile[ FILENAME_MAX ];
	char szGravityModel[ FILENAME_MAX ];
	int nMaxGravityDegree;

	// Initial Orbital Elements (IOE)
	int nIOEType;	// 1 for IRV, 2 for TLE, 3 for EOSMOE
	char szIOEFile[ FILENAME_MAX ];

	char szCatalogFile[ FILENAME_MAX ];

	// visibility condition
	stVisibilityCondition stVisibility;
};

extern stVisComControl g_stControlFile;

/***********************************************************************************

 Orbital elements that define a satellite orbit

***********************************************************************************/
/*
struct stSatelliteIOE
{
	int nSatelliteID;
	int nSIC;
	char cElementType;			// T for TLE, E for EOS, I for IRV
	int nIntJD, nFractionJD;

	int pnElement1to6[ 6 ];
	float pfElement7to18[ 20 ];

	double GetRefJD() { return ( (double)nIntJD + (double)nFractionJD * 1.0e-9 ); }
	void SetRefJD( double dfJD )
	{
		nIntJD = (int) dfJD;
		nFractionJD = (int) ( ( dfJD - nIntJD ) * 1.0e9 );
	}

	void GetPosVel( double *pdfPos, double *pdfVel )
	{
		for( int i = 0; i < 3; i++ )
		{
			pdfPos[ i ] = (double)pnElement1to6[ i ] * 1.e-1;
			pdfVel[ i ] = (double)pnElement1to6[ i + 3 ] * 1.0e-4;
		}
	}

	double GetOrbitalPeriod() { return (double)pfElement7to18[ 11 ]; }
};
*/

/********************************************************************************

  It is fast to use numerical integrator to compute the position/velcoity of 
  satellite. So, in the Schedule software, it is designed that all types of
  satellite orbit elements are converted to position/velocity at a reference 
  time. With the time going on, the reference time is pushed forward.

*********************************************************************************/
struct stInternalIRV
{
	int nSatelliteID;
	double dfRefJD;
	double pdfPos[ 3 ], pdfVel[ 3 ];
};

/********************************************************************************

  Parameters that control the element conversion

*********************************************************************************/

struct stElementConversionControl
{
	char nSourceType;	// 1 for IRV, 2 for TLE, 3 for EOSMOE
	char szSourceFile[ FILENAME_MAX ];

	char nResultType;
	char szResultFile[ FILENAME_MAX ];

	// reference epoch
	int nYear, nMonth, nDay, nHour;	// in UTC
	double dfJD;	// from Year/Month/Day/Hour	in UTC

};

extern stElementConversionControl g_stElementConversionControl;

/********************************************************************************

  Star catalog data of a star

*********************************************************************************/

struct stCatalogStar
{
	int nNoCatalogue;
	char cSource;		// F for FK5, ...
	double dfMagnitude;
	double dfRA, dfDEC, dfPara, dfPMRA, dfPMDEC, dfRV;
};

/********************************************************************************

  Data that defines an object on an image

*********************************************************************************/

struct stImageObject
{
	int nNo;
	double dfX, dfY;	// origin is in the centre of the image
	double dfIntensity;
};

struct stPV
{
	double JD;
	double pos[3], vel[3];

#ifndef NDEBUG
	int year, month, day, hour, minute;
	double second;
#endif

};


#endif