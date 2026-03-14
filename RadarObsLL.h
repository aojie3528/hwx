#ifndef RADAR_OBS_LL_H
#define RADAR_OBS_LL_H

#include "DataStructure.h"
#include <map>


class cRadarObsLL
{
public:
	bool readObsFile(const char* fileName, vector<stSiteInfo>* m_Sites, RadarArc* tArc);
	bool readObsFile(const char *fileName, obsPass1* pass);
	bool readObsFile1(const char* fileName, obsPass1* pass);
	bool readObsFile2(const char* fileName, obsPass1* pass);

	bool readObsFromMap(string fileName, RadarArc* tArc,
		map<string, vector<epochPosVel1>>* ArcObsMap);
	int JudgeObsType(string ObsPath);
	bool GetArcMag(std::vector<epochPosVel1>* pObsPos, double& dMag);

	bool UseLessPoint(obsPass1& ObsPass, int nPoint);
	bool pass2Arc(std::vector<stSiteInfo>* m_Sites, obsPass1* pass, RadarArc* Arc);
	bool pass2Arc_GroundOpt(vector<stSiteInfo>* m_Sites, obsPass1* pass, RadarArc* Arc);
	bool pass2Arc_SpaceOpt(vector<stSiteInfo>* m_Sites, obsPass1* pass, RadarArc* Arc);

	bool ReadSiteInfo(const char *SiteFile, std::vector<stSiteInfo>& Sites);
	bool ReadSBOSetInfo(const char* SBOFile, std::vector<stSiteInfo>& Sites);


	bool ReadArcFile(const char* FileName, int ObsType, vector<epochPosVel1>& ArcInfo);
	bool ReadArcFile(const char* FileName, std::vector<epochPosVel1>& ArcInfo);
	bool ReadArcFile(const char *FileName, RadarArc& Arc);
	bool ReadArcFile(const char* FileName, int ObsType, RadarArc& Arc);

	bool GetArcRCSMag(RadarArc& Arc);
	bool LDOutlierDetection(int n, double *JD, double *x, bool *mark);

	//RAE文件
	bool readRAEFile(const char *fileName, obsPass1* pass);

	//odp文件
	bool readodpFile(const char *fileName, obsPass1* pass);
	bool readodpFile_GroundOpt(const char* fileName, obsPass1* pass);
	bool readodpFile_SpaceOpt(const char* fileName, obsPass1* pass);

	bool readBITFile(const char* fileName, obsPass1* pass);

	//GTW文件
	bool readGTWFile(const char* filename, obsPass1* pass);
	bool ReadGTWLine(std::string filename, std::string line, double& JD, double& X, double& Y, double& Z, double& RA, double& Dec, double& Mag);

	bool ConvertRaDec2TOD(double JDUTC, double& Ra, double& Dec);
	bool ConvertRaDec2J2000(double JDUTC, double& Ra, double& Dec);

	bool computeUnitVectorfromRADec(double ra, double dec, double* v);
	bool computeRADecFromPosVector(double* v, double& ra, double& dec);

	bool ConvertTOD2J2000(double JDUTC, double* TODPos, double* J2000Pos, bool isTOD2J2000);
	bool ECEF2TOD(double JDUTC, double* ECEFPos, double* TODPos);
	double GetMedian(vector<double>& Numbers);

	std::string rtrimtab(std::string& s);
};

#endif