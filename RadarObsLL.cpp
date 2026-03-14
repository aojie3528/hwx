#include "RadarObsLL.h"
#include <string>
#include <fstream>
#include <sstream>
#include "EOP.h"

#if (_MSVC_LANG > 201402L) || (__cplusplus >= 201703L) // C++ 17
#include <filesystem>
namespace fs = std::filesystem;
#else
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem::v1;
#endif
#include <regex>

#include "Constant.h"
#include "DateTimeZ.h"
#include "SimulateDataLL.h"
#include "GreenwichSiderealTime.h"
#include "COORTRAN.H"
#include "PolyFit.h"
#include "TransformMatrix.h"
#include "NPElement.h"
//#include "utils.hpp"

using namespace std;

bool cRadarObsLL::readRAEFile(const char *fileName, obsPass1 *pass)
{
	try
	{
		strcpy(pass->fileName, fileName);
		pass->noradID = -1;
		pass->az.clear();
		pass->el.clear();
		pass->dis.clear();
		pass->rcs_all.clear();
		pass->JD.clear();
		// 01808 1407 107652 1376 8 1 13 1 20200126 023103175624 00000073584 055349 1100 000 98 0865 103 1 0030 4803 2314 1375

		bool firstpoint = true;
		string line, sub;

		ifstream in(fileName);
		if (!in.is_open())
			return false;

		while (getline(in, line))
		{
			if (line.length() < 50)
				continue;
			int temp = 0, ID = 0, year = 0, month = 0, day = 0, hour = 0, minute = 0;
			double second = 0.0, az = 0.0, el = 0.0, dis = 0.0, tJD = 0.0;

			sub = line.substr(18, 4);
			pass->equipIDr = atoi(sub.c_str());

			sub = line.substr(23, 1);
			int TimeFlag = atoi(sub.c_str());

			sub = line.substr(27, 2);
			ID = atoi(sub.c_str());

			sub = line.substr(32, 4);
			year = atoi(sub.c_str());

			sub = line.substr(36, 2);
			month = atoi(sub.c_str());

			sub = line.substr(38, 2);
			day = atoi(sub.c_str());

			sub = line.substr(41, 2);
			hour = atoi(sub.c_str());

			sub = line.substr(43, 2);
			minute = atoi(sub.c_str());

			sub = line.substr(45, 8);
			second = atof(sub.c_str()) / 1000000.0;

			if (13 == ID)
			{
				sub = line.substr(54, 11);
				dis = atoi(sub.c_str()) * 1000;

				sub = line.substr(66, 6);
				dis += atof(sub.c_str()) / 1000.0;
				if (dis < 0)
					continue;
				// if (dis >80000000.0)continue;

				// 读入RCS
				sub = line.substr(101, 5);
				int nRCS = atoi(sub.c_str());
				int XRCS = nRCS / 10;
				int YRCS = nRCS % 10;
				pass->rcs_all.emplace_back(XRCS * 10e-4 * pow(10.0, YRCS - 1));

				// getline(in, line);
				// getline(in, line);

				getline(in, line);
				sub = line.substr(27, 2);
				ID = atoi(sub.c_str());

				if (23 == ID)
				{
					getline(in, line);
					sub = line.substr(27, 2);
					ID = atoi(sub.c_str());
				}

				if (31 == ID)
				{
					sub = line.substr(54, 3);
					az = atof(sub.c_str());
					sub = line.substr(57, 2);
					az += atof(sub.c_str()) / 60.0;
					sub = line.substr(59, 4);
					az += atof(sub.c_str()) / 100 / 3600.0;

					sub = line.substr(64, 1);
					int sym = atoi(sub.c_str());

					sub = line.substr(65, 2);
					el = atof(sub.c_str());
					sub = line.substr(67, 2);
					el += atof(sub.c_str()) / 60.0;
					sub = line.substr(69, 3);
					el += atof(sub.c_str()) / 10 / 3600.0;

					if (sym == 1)
						el = -el;

					sub = line.substr(111, 4);
					pass->equipIDs = atoi(sub.c_str());

					if (TimeFlag == 8)
						g_DateTimeZ.DateTime2JD(year, month, day, hour - 8, minute, second, tJD);
					else
						g_DateTimeZ.DateTime2JD(year, month, day, hour, minute, second, tJD);

					// double delta = dis / g_dfLightSpeed / 86400.0;
					// fraJD += delta;

					az *= g_dfDEG2RAD;
					el *= g_dfDEG2RAD;

					pass->JD.push_back(tJD);
					pass->az.push_back(az);
					pass->el.push_back(el);
					pass->dis.push_back(dis);
				}
			}

			if (11 == ID)
			{
				sub = line.substr(54, 11);
				dis = atoi(sub.c_str()) * 1000;

				sub = line.substr(66, 6);
				dis += atof(sub.c_str()) / 1000.0;

				dis *= 2.0;
				if (dis < 0)
					continue;
				// if (dis >80000000.0)continue;

				// 读入RCS
				sub = line.substr(101, 5);
				int nRCS = atoi(sub.c_str());
				int XRCS = nRCS / 10;
				int YRCS = nRCS % 10;
				pass->rcs_all.emplace_back(XRCS * 10e-4 * pow(10.0, YRCS - 1));

				getline(in, line);
				sub = line.substr(27, 2);
				ID = atoi(sub.c_str());

				if (21 == ID)
				{
					getline(in, line);
					sub = line.substr(27, 2);
					ID = atoi(sub.c_str());
				}

				if (31 == ID)
				{
					sub = line.substr(54, 3);
					az = atof(sub.c_str());
					sub = line.substr(57, 2);
					az += atof(sub.c_str()) / 60.0;
					sub = line.substr(59, 4);
					az += atof(sub.c_str()) / 100 / 3600.0;

					sub = line.substr(64, 1);
					int sym = atoi(sub.c_str());

					sub = line.substr(65, 2);
					el = atof(sub.c_str());
					sub = line.substr(67, 2);
					el += atof(sub.c_str()) / 60.0;
					sub = line.substr(69, 3);
					el += atof(sub.c_str()) / 10 / 3600.0;

					if (sym == 1)
						el = -el;

					pass->equipIDs = pass->equipIDr;

					if (TimeFlag == 8)
						g_DateTimeZ.DateTime2JD(year, month, day, hour - 8, minute, second, tJD);
					else
						g_DateTimeZ.DateTime2JD(year, month, day, hour, minute, second, tJD);

					// double delta = dis / g_dfLightSpeed / 86400.0;
					// fraJD += delta;

					az *= g_dfDEG2RAD;
					el *= g_dfDEG2RAD;

					if (firstpoint)
					{
						pass->maxd = pass->mind = dis;
						firstpoint = false;
					}
					pass->maxd = max(dis, pass->maxd);
					pass->mind = min(dis, pass->mind);
					pass->JD.push_back(tJD);
					pass->az.push_back(az);
					pass->el.push_back(el);
					pass->dis.push_back(dis);
				}
			}
		}

		pass->numobs = static_cast<int>(pass->az.size());
		if (pass->numobs == 0)
		{
			return false;
		}

		// 计算RCS中位值
		pass->RCS = GetMedian(pass->rcs_all);
		return true;
	}
	catch (...)
	{
		printf("wrong file: %s\n", fileName);
		return false;
	}
}

// read Obs to obspass
bool cRadarObsLL::readObsFile1(const char *fileName, obsPass1 *pass)
{
	try
	{
		strcpy(pass->fileName, fileName);
		pass->noradID = -1;
		pass->az.clear();
		pass->el.clear();
		pass->dis.clear();
		pass->JD.clear();
		// 01808 1407 107652 1376 8 1 13 1 20200126 023103175624 00000073584 055349 1100 000 98 0865 103 1 0030 4803 2314 1375
		bool firstpoint = true;
		string line, sub;

		// by LL RCS read
		vector<double> AllRCSs;

		ifstream in(fileName);
		if (!in.is_open())
			return false;

		string tName(fileName);
		sub = tName.substr(tName.find_last_of("_") + 1, 4);
		pass->equipIDr = atoi(sub.c_str());

		while (getline(in, line))
		{
			double JDUTC = 0.0;
			int temp = 0, ID = 0, year = 0, month = 0, day = 0, hour = 0, minute = 0;
			double second = 0.0, az = 0.0, el = 0.0, dis = 0.0, tJD = 0.0;
			double x, y, z;

			stringstream ss(line);
			ss >> JDUTC >> az >> el >> dis >> x >> y >> z;

			az *= g_dfDEG2RAD;
			el *= g_dfDEG2RAD;

			if (firstpoint)
			{
				pass->maxd = pass->mind = dis;
				firstpoint = false;
			}
			pass->maxd = max(dis, pass->maxd);
			pass->mind = min(dis, pass->mind);

			tJD = JDUTC;

			pass->JD.push_back(tJD);
			pass->az.push_back(az);
			pass->el.push_back(el);
			pass->dis.push_back(dis);
		}

		pass->numobs = static_cast<int>(pass->az.size());
		if (pass->numobs == 0)
		{
			return false;
		}

		// 计算RCS中位值
		pass->RCS = 0.0;
		pass->m_Type = 0;
		return true;
	}
	catch (...)
	{
		printf("wrong file: %s\n", fileName);
		return false;
	}
}

// read OObs to obspass
bool cRadarObsLL::readObsFile2(const char *fileName, obsPass1 *pass)
{
	try
	{
		strcpy(pass->fileName, fileName);
		pass->noradID = -1;
		pass->ra.clear();
		pass->dec.clear();
		pass->dis.clear();
		pass->JD.clear();
		bool firstpoint = true;
		string line, sub;

		ifstream in(fileName);
		if (!in.is_open())
			return false;

		string tName(fileName);
		sub = tName.substr(tName.find_last_of("_") + 1, 4);
		pass->equipIDr = atoi(sub.c_str());

		pass->m_Type = 1;
		if (pass->equipIDr / 1000 == 9)
			pass->m_Type = 2;

		while (getline(in, line))
		{
			double JDUTC = 0.0;
			int temp = 0, ID = 0, year = 0, month = 0, day = 0, hour = 0, minute = 0;
			double second = 0.0, ra = 0.0, dec = 0.0, dis = 0.0, tJD = 0.0;
			double satx, saty, satz;

			stringstream ss(line);
			ss >> JDUTC >> ra >> dec >> dis >> satx >> saty >> satz;

			ra *= g_dfDEG2RAD;
			dec *= g_dfDEG2RAD;

			if (firstpoint)
			{
				pass->maxd = pass->mind = dis;
				firstpoint = false;
			}
			pass->maxd = max(dis, pass->maxd);
			pass->mind = min(dis, pass->mind);

			tJD = JDUTC;

			pass->JD.push_back(tJD);
			pass->ra.push_back(ra);
			pass->dec.push_back(dec);
			pass->dis.push_back(dis);
			pass->satx.push_back(satx);
			pass->saty.push_back(saty);
			pass->satz.push_back(satz);
		}

		pass->numobs = static_cast<int>(pass->ra.size());
		if (pass->numobs == 0)
		{
			return false;
		}

		// 计算RCS中位值
		pass->RCS = 0.0;
		return true;
	}
	catch (...)
	{
		printf("wrong file: %s\n", fileName);
		return false;
	}
}

bool cRadarObsLL::readGTWFile(const char *filename, obsPass1 *pass)
{
	auto gtwpath = fs::path(filename);

	if (!fs::exists(gtwpath))
	{
		printf("file %s not exist!\n", filename);
		return false;
	}

	strcpy(pass->fileName, filename);
	pass->noradID = -1;
	pass->az.clear();
	pass->el.clear();
	pass->dis.clear();
	pass->ra.clear();
	pass->dec.clear();
	pass->satx.clear();
	pass->saty.clear();
	pass->satz.clear();
	pass->mag_all.clear();
	pass->JD.clear();

	std::ifstream in(gtwpath);
	std::string line;
	int num = 0;
	while (getline(in, line))
	{
		// if (line.at(0) != 'C' && line.length() > 35)
		if (line.at(0) != 'C' && line.length() > 195) // CLY
		{
			double dfJD = 0;
			double sat_pos[3]{0};
			double ra = 0;
			double dec = 0;
			double mag = 0;
			if (!ReadGTWLine(filename, line, dfJD, sat_pos[0], sat_pos[1], sat_pos[2], ra, dec, mag))
				return false;

			pass->JD.emplace_back(dfJD);
			pass->satx.emplace_back(sat_pos[0]);
			pass->saty.emplace_back(sat_pos[1]);
			pass->satz.emplace_back(sat_pos[2]);
			pass->ra.emplace_back(ra);
			pass->dec.emplace_back(dec);
			pass->dis.emplace_back(0.0);
			pass->mag_all.emplace_back(mag);

			num++;
		}
	}
	in.close();

	if (num <= 0)
		printf("文件内容为空 %s\n", filename);

	pass->numobs = static_cast<int>(pass->JD.size());
	if (pass->numobs <= 0)
		return false;

	// 计算Mag中位值
	pass->Mag = GetMedian(pass->mag_all);

	string sub = string(filename);
	sub = sub.substr(sub.find_last_of("/\\") + 1);
	sub = sub.substr(22, 4);
	pass->equipIDr = atoi(sub.c_str());

	// 判断天基地基
	double tSat_pos[3] = {pass->satx[0], pass->saty[0], pass->satz[0]};
	if (g_Matrix.Vec_norm2(tSat_pos, 3) < 6400000.0)
		pass->m_Type = 1;
	else
		pass->m_Type = 2;

	return true;
}

bool cRadarObsLL::ReadGTWLine(std::string filename, std::string line, double &JD,
							  double &X, double &Y, double &Z, double &RA, double &Dec, double &Mag)
{
	if (line.at(0) != 'C' && line.length() > 35)
	{
		int noradID = std::atoi(line.substr(11, 6).c_str());
		int equipID = std::atoi(line.substr(18, 4).c_str());
		char time_flag = line[23];
		int year = std::atoi(line.substr(32, 4).c_str());
		int month = std::atoi(line.substr(36, 2).c_str());
		int day = std::atoi(line.substr(38, 2).c_str());
		int hour = std::atoi(line.substr(41, 2).c_str());
		int minute = std::atoi(line.substr(43, 2).c_str());
		double second = std::atoll(line.substr(45, 8).c_str()) * 1.0e-6;

		double dfJD = 0;
		g_DateTimeZ.DateTime2JD(year, month, day, hour, minute, second, dfJD);
		if (time_flag == '8') // beijing time to UTC
			dfJD -= 8.0 / 24.0;

		double rad = std::atoi(line.substr(54, 3).c_str());
		double ram = std::atoi(line.substr(57, 2).c_str());
		double ras = std::atoi(line.substr(59, 4).c_str()) * 0.01;
		double ra = (rad + ram / 60.0 + ras / 3600.0) * g_dfDEG2RAD;

		int decsign = std::atoi(line.substr(64, 1).c_str()) == 0 ? 1 : -1;
		double decd = std::atoi(line.substr(65, 2).c_str());
		double decm = std::atoi(line.substr(67, 2).c_str());
		double decs = std::atoi(line.substr(69, 3).c_str()) * 0.1;
		double dec = (decd + decm / 60.0 + decs / 3600.0) * g_dfDEG2RAD * decsign;

		int magsign = std::atoi(line.substr(97, 1).c_str()) == 0 ? 1 : -1;
		double mag = std::atoi(line.substr(98, 3).c_str()) / 10.0 * magsign;

		if (dfJD < g_dfMJDRef || dfJD > g_dfJ2000 + 1000 * 365.25)
		{
			printf("GTW文件时间解析错误 %s\n", filename.c_str());
			return false;
		}

		if (ra < -0.1 || ra > g_dfTWOPI + 0.1)
		{
			printf("GTW文件RA解析错误 %s\n", filename.c_str());
			return false;
		}

		if (dec < -g_dfHALFPI - 0.1 || dec > g_dfHALFPI + 0.1)
		{
			printf("GTW文件DEC解析错误 %s\n", filename.c_str());
			return false;
		}

		if (!ConvertRaDec2TOD(dfJD, ra, dec))
			return false;

		double sat_pos[3]{0};
		double sat_pos_tod[3]{0};
		// 天基读取卫星位置
		if (rtrimtab(line).size() > 112)
		{
			int sign = std::atoi(line.substr(114, 1).c_str()) == 0 ? 1 : -1;
			sat_pos[0] = std::atoll(line.substr(115, 11).c_str()) * 0.001 * sign;
			sign = std::atoi(line.substr(127, 1).c_str()) == 0 ? 1 : -1;
			sat_pos[1] = std::atoll(line.substr(128, 11).c_str()) * 0.001 * sign;
			sign = std::atoi(line.substr(140, 1).c_str()) == 0 ? 1 : -1;
			sat_pos[2] = std::atoll(line.substr(141, 11).c_str()) * 0.001 * sign;

			if (std::sqrt(std::pow(sat_pos[0], 2) + std::pow(sat_pos[1], 2) + std::pow(sat_pos[2], 2)) < g_dfEarthSemiMajor)
			{
				printf("GTW文件卫星位置解析错误 %s\n", filename.c_str());
				return false;
			}
			if (!ConvertTOD2J2000(dfJD, sat_pos_tod, sat_pos, false))
				return false;
			// trans.TransJ20002TOD(jd_tt, sat_pos);
			magsign = std::atoi(line.substr(190, 1).c_str()) == 0 ? 1 : -1;
			mag = std::atoi(line.substr(191, 5).c_str()) / 100.0 * magsign;
		}

		JD = dfJD;
		X = sat_pos_tod[0];
		Y = sat_pos_tod[1];
		Z = sat_pos_tod[2];
		RA = ra;
		Dec = dec;
		Mag = mag;
	}

	return true;
}

std::string cRadarObsLL::rtrimtab(std::string& s)
{
	if (!s.empty() && s.back() == '\r')  // Check if the last character is '\r'
		s.pop_back();  // Remove the last character
	return s;
}

// 读入北京理工大学格式
bool cRadarObsLL::readBITFile(const char *fileName, obsPass1 *pass)
{
	ifstream fin;
	fin.open(fileName);
	string tLine;
	while (getline(fin, tLine))
	{
		if (tLine.length() < 50)
			continue;

		int nYear, nMon, nDay, nHour, nMin, nSec;
		double Dis, Az, El;
		stringstream ss(tLine);
		ss >> nYear >> nMon >> nDay >> nHour >> nMin >> nSec >> Dis >> Az >> El;

		double tJD = 0.0;
		g_DateTimeZ.DateTime2JD(nYear, nMon, nDay, nHour, nMin, (double)nSec, tJD);

		pass->JD.emplace_back(tJD);
		pass->dis.emplace_back(Dis * 1000.0 * 2.0);
		pass->az.emplace_back(Az);
		pass->el.emplace_back(El);
	}
	pass->numobs = static_cast<int>(pass->JD.size());
	sprintf(pass->fileName, "%s", fileName);

	pass->equipIDr = 0001;
	pass->equipIDs = 0001;

	return true;
}

bool cRadarObsLL::readodpFile(const char *fileName, obsPass1 *pass)
{
	strcpy(pass->fileName, fileName);
	pass->noradID = -1;
	pass->az.clear();
	pass->el.clear();
	pass->dis.clear();
	pass->JD.clear();

	bool firstpoint = true;
	string line, sub;

	ifstream in(fileName);
	if (!in.is_open())
		return false;

	bool isPair = false;
	string lastTime;
	double taz = 0.0, tel = 0.0, tDis = 0.0;

	while (getline(in, line))
	{
		if (line.length() < 50)
			continue;
		int temp = 0, ID = 0, year = 0, month = 0, day = 0, hour = 0, minute = 0;
		double second = 0.0, az = 0.0, el = 0.0, dis = 0.0, tJD = 0.0;

		sub = line.substr(7, 2);
		ID = atoi(sub.c_str());

		sub = line.substr(11, 5);
		pass->equipIDr = atoi(sub.c_str());

		sub = line.substr(16, 2);
		year = atoi(sub.c_str()) + 2000;

		sub = line.substr(18, 3);
		day = atoi(sub.c_str());

		sub = line.substr(21, 5);
		second = atof(sub.c_str());

		sub = line.substr(26, 6);
		second += atof(sub.c_str()) / 1000000.0;

		string tTime = line.substr(16, 10);

		if (tTime == lastTime)
			isPair = true;
		else
			isPair = false;

		if (70 == ID)
		{
			sub = line.substr(35, 3);
			az = atof(sub.c_str());
			sub = line.substr(38, 2);
			az += atof(sub.c_str()) / 60.0;
			sub = line.substr(40, 5);
			az += atof(sub.c_str()) / 1000 / 3600.0;

			sub = line.substr(45, 3);
			el = atof(sub.c_str());
			sub = line.substr(48, 2);
			el += atof(sub.c_str()) / 60.0;
			sub = line.substr(50, 5);
			el += atof(sub.c_str()) / 1000 / 3600.0;

			taz = az * g_dfDEG2RAD;
			tel = el * g_dfDEG2RAD;

			// if (el < 5.0)continue;
		}

		if (21 == ID)
		{
			sub = line.substr(39, 6);
			dis = atof(sub.c_str()) * 1000;
			sub = line.substr(45, 6);
			dis += atof(sub.c_str()) / 1000;

			if (dis < 0)
				continue;
			if (dis > 80000000.0)
				continue;

			sub = line.substr(96, 4);
			pass->equipIDs = atoi(sub.c_str());

			tDis = dis;
		}

		if (20 == ID)
		{
			sub = line.substr(39, 6);
			dis = atof(sub.c_str()) * 1000;
			sub = line.substr(45, 6);
			dis += atof(sub.c_str()) / 1000;

			if (dis < 0)
				continue;
			if (dis > 80000000.0)
				continue;

			pass->equipIDs = pass->equipIDr;

			tDis = dis * 2.0;
		}

		g_DateTimeZ.DateTime2JD(year, 1, day, 0, 0, second, tJD);

		if (isPair)
		{
			pass->JD.push_back(tJD);
			pass->az.push_back(taz);
			pass->el.push_back(tel);
			pass->dis.push_back(tDis);
			lastTime = "";
			// printf("%.3f  %.10f   %.3f   %.3f   %.3f\n",second, intJD+fraJD, taz * g_dfRAD2DEG, tel * g_dfRAD2DEG, tDis);
		}
		else
		{
			lastTime = tTime;
		}
	}

	pass->numobs = static_cast<int>(pass->az.size());

	if (pass->numobs <= 0)
		return false;

	return true;
}

bool cRadarObsLL::readodpFile_GroundOpt(const char *fileName, obsPass1 *pass)
{
	strcpy(pass->fileName, fileName);
	pass->noradID = -1;
	pass->az.clear();
	pass->el.clear();
	pass->dis.clear();
	pass->ra.clear();
	pass->dec.clear();
	pass->JD.clear();

	pass->m_Type = 1;

	bool firstpoint = true;
	string line, sub;

	ifstream in(fileName);
	if (!in.is_open())
		return false;

	bool isPair = false;
	string lastTime;
	double tra = 0.0, tdec = 0.0;

	while (getline(in, line))
	{
		if (line.length() < 50)
			continue;
		int temp = 0, ID = 0, year = 0, month = 0, day = 0, hour = 0, minute = 0;
		double second = 0.0, ra = 0.0, dec = 0.0, tJD = 0.0;

		sub = line.substr(7, 2);
		ID = atoi(sub.c_str());

		sub = line.substr(11, 5);
		pass->equipIDr = atoi(sub.c_str());

		sub = line.substr(16, 2);
		year = atoi(sub.c_str()) + 2000;

		sub = line.substr(18, 3);
		day = atoi(sub.c_str());

		sub = line.substr(21, 5);
		second = atof(sub.c_str());

		sub = line.substr(26, 6);
		second += atof(sub.c_str()) / 1000000.0;

		string tTime = line.substr(16, 10);

		if (10 == ID)
		{
			sub = line.substr(35, 3);
			ra = atof(sub.c_str());
			sub = line.substr(38, 2);
			ra += atof(sub.c_str()) / 60.0;
			sub = line.substr(40, 5);
			ra += atof(sub.c_str()) / 1000 / 3600.0;

			sub = line.substr(45, 3);
			dec = atof(sub.c_str());
			sub = line.substr(48, 2);
			dec += atof(sub.c_str()) / 60.0;
			sub = line.substr(50, 5);
			dec += atof(sub.c_str()) / 1000 / 3600.0;

			tra = ra * g_dfDEG2RAD;
			tdec = dec * g_dfDEG2RAD;
		}

		g_DateTimeZ.DateTime2JD(year, 1, day, 0, 0, second, tJD);

		ConvertRaDec2TOD(tJD, tra, tdec);

		pass->JD.push_back(tJD);
		pass->ra.push_back(tra);
		pass->dec.push_back(tdec);
	}

	pass->numobs = static_cast<int>(pass->ra.size());

	if (pass->numobs <= 0)
		return false;

	return true;
}

bool cRadarObsLL::readodpFile_SpaceOpt(const char *fileName, obsPass1 *pass)
{
	strcpy(pass->fileName, fileName);
	pass->noradID = -1;
	pass->az.clear();
	pass->el.clear();
	pass->dis.clear();
	pass->ra.clear();
	pass->dec.clear();
	pass->JD.clear();

	pass->m_Type = 2;

	bool firstpoint = true;
	string line, sub;

	ifstream in(fileName);
	if (!in.is_open())
		return false;

	bool isPair = false;
	string lastTime;
	double tra = 0.0, tdec = 0.0;

	while (getline(in, line))
	{
		if (line.length() < 50)
			continue;
		int temp = 0, ID = 0, year = 0, month = 0, day = 0, hour = 0, minute = 0;
		double second = 0.0, ra = 0.0, dec = 0.0, tJD = 0.0;

		sub = line.substr(7, 2);
		ID = atoi(sub.c_str());

		sub = line.substr(11, 5);
		pass->equipIDr = atoi(sub.c_str());

		sub = line.substr(16, 2);
		year = atoi(sub.c_str()) + 2000;

		sub = line.substr(18, 3);
		day = atoi(sub.c_str());

		sub = line.substr(21, 5);
		second = atof(sub.c_str());

		sub = line.substr(26, 6);
		second += atof(sub.c_str()) / 1000000.0;

		string tTime = line.substr(16, 10);

		if (10 == ID)
		{
			sub = line.substr(35, 3);
			ra = atof(sub.c_str());
			sub = line.substr(38, 2);
			ra += atof(sub.c_str()) / 60.0;
			sub = line.substr(40, 5);
			ra += atof(sub.c_str()) / 1000 / 3600.0;

			sub = line.substr(45, 3);
			dec = atof(sub.c_str());
			sub = line.substr(48, 2);
			dec += atof(sub.c_str()) / 60.0;
			sub = line.substr(50, 5);
			dec += atof(sub.c_str()) / 1000 / 3600.0;

			tra = ra * g_dfDEG2RAD;
			tdec = dec * g_dfDEG2RAD;
		}

		g_DateTimeZ.DateTime2JD(year, 1, day, 0, 0, second, tJD);

		ConvertRaDec2TOD(tJD, tra, tdec);

		pass->JD.push_back(tJD);
		pass->ra.push_back(tra);
		pass->dec.push_back(tdec);
	}

	pass->numobs = static_cast<int>(pass->ra.size());

	if (pass->numobs <= 0)
		return false;

	return true;
}

bool cRadarObsLL::ConvertRaDec2TOD(double JDUTC, double &Ra, double &Dec)
{
	double v1[3], v2[3];
	computeUnitVectorfromRADec(Ra, Dec, v1);
	// ra, dec from J2000 to TOD
	if (!ConvertTOD2J2000(JDUTC, v2, v1, false))
	{
		return false;
	}
	computeRADecFromPosVector(v2, Ra, Dec);

	return true;
}

bool cRadarObsLL::ConvertRaDec2J2000(double JDUTC, double &Ra, double &Dec)
{
	double v1[3], v2[3];
	computeUnitVectorfromRADec(Ra, Dec, v1);
	// ra, dec from J2000 to TOD
	if (!ConvertTOD2J2000(JDUTC, v1, v2, true))
	{
		return false;
	}
	computeRADecFromPosVector(v2, Ra, Dec);

	return true;
}

// 调用函数时，TODPos与J2000Pos不能是同一变量
//  by LL 20250620
bool cRadarObsLL::ConvertTOD2J2000(double JDUTC, double *TODPos, double *J2000Pos, bool isTOD2J2000)
{
	double jd_et{0};
	if (!g_EOP.ConvertETUTC(jd_et, JDUTC, false))
	{
		printf("EOP not initialized!\n");
		return false;
	}

	double matrix_convert[9]{0};
	cNPElement np;
	if (isTOD2J2000)
	{
		if (!np.GetMatrixTOD2J2000(jd_et, matrix_convert))
			return false;
	}
	else
	{
		if (!np.GetMatrixJ20002TOD(jd_et, matrix_convert))
			return false;
	}

	// 计算TOD位置
	if (isTOD2J2000)
	{
		g_Matrix.Mat_Multiply_Vec(matrix_convert, TODPos, J2000Pos, 3, 3);
	}
	else
	{
		g_Matrix.Mat_Multiply_Vec(matrix_convert, J2000Pos, TODPos, 3, 3);
	}
	return true;
}

bool cRadarObsLL::computeUnitVectorfromRADec(double ra, double dec, double *v)
{
	v[0] = cos(ra) * cos(dec);
	v[1] = sin(ra) * cos(dec);
	v[2] = sin(dec);

	return true;
}

bool cRadarObsLL::computeRADecFromPosVector(double *v, double &ra, double &dec)
{
	double dis = sqrt(g_Matrix.Vec_Multiply_Vec(v, v, 3));

	dec = asin(v[2] / dis);

	ra = atan2(v[1], v[0]);
	if (ra < 0.0)
		ra += g_dfTWOPI;

	return true;
}

bool cRadarObsLL::readObsFile(const char *fileName, obsPass1 *pass)
{
	// 判断文件类型
	string tName(fileName);
	string tType = tName.substr(tName.find_last_of(".") + 1);
	if (tType == "RAE")
	{
		return readRAEFile(fileName, pass);
	}
	else if (tType == "odp")
	{
		string tName(fileName);
		tName = tName.substr(tName.find_last_of("/\\") + 1);
		tName = tName.substr(22, 1);
		if (tName == "7")
			return readodpFile_GroundOpt(fileName, pass);
		else if (tName == "9")
			return readodpFile_SpaceOpt(fileName, pass);
		else
			return readodpFile(fileName, pass);
	}
	else if (tType == "GTW")
	{
		return readGTWFile(fileName, pass);
	}
	else if (tType == "txt")
	{
		return readBITFile(fileName, pass);
	}
	else if (tType == "Obs")
	{
		return readObsFile1(fileName, pass);
	}
	else if (tType == "OObs")
	{
		return readObsFile2(fileName, pass);
	}
	return false;
}

bool cRadarObsLL::pass2Arc(vector<stSiteInfo> *m_Sites, obsPass1 *pass, RadarArc *Arc)
{
	cSimulateDataLL SD;
	stSiteInfo mySite;
	bool isSite = false;
	for (int i = 0; i < m_Sites->size(); i++)
	{
		if ((*m_Sites)[i].m_IDr == pass->equipIDr)
		{
			mySite = (*m_Sites)[i];
			isSite = true;
			break;
		}
	}
	if (!isSite)
		return false;

	epochPosVel1 pv;
	int nobs = static_cast<int>(pass->az.size());
	double ECEFPos[3];

	// 判断测站类型，不同处理方案
	if (mySite.m_StationType == 1)
	{
		for (int j = 0; j < nobs; j++)
		{
			pv.JDUTC = pass->JD[j] /*+ pass.dis[j] / g_dfLightSpeed / 86400.0*/;

			if (!SD.NewtonRaphson_SR_RadarPositioning(mySite.m_siteECEFs, mySite.m_siteECEFr, mySite.m_siteBLHr, pass->az[j], pass->el[j], pass->dis[j], ECEFPos))
				continue;

			double dfGST = cGreenwichST::ComputeGST(pv.JDUTC);
			double dfSin = sin(dfGST), dfCos = cos(dfGST);

			pv.pos[0] = ECEFPos[0] * dfCos - ECEFPos[1] * dfSin;
			pv.pos[1] = ECEFPos[0] * dfSin + ECEFPos[1] * dfCos;
			pv.pos[2] = ECEFPos[2];

			// add by LL
			pv.Az = pass->az[j];
			pv.El = pass->el[j];
			pv.DisTwoway = pass->dis[j];
			pv.RCS = pass->rcs_all[j];

			Arc->m_ObsPos.push_back(pv);
		}
	}
	else
	{
		for (int j = 0; j < nobs; j++)
		{
			pv.JDUTC = pass->JD[j] /*+ pass.dis[j] / g_dfLightSpeed / 86400.0*/;

			if (!SD.SingleEquip_RadarPositioning(mySite.m_siteECEFr, mySite.m_siteBLHr, pass->az[j], pass->el[j], pass->dis[j], ECEFPos))
				continue;

			double dfGST = cGreenwichST::ComputeGST(pv.JDUTC);
			double dfSin = sin(dfGST), dfCos = cos(dfGST);

			pv.pos[0] = ECEFPos[0] * dfCos - ECEFPos[1] * dfSin;
			pv.pos[1] = ECEFPos[0] * dfSin + ECEFPos[1] * dfCos;
			pv.pos[2] = ECEFPos[2];

			// add by LL
			pv.Az = pass->az[j];
			pv.El = pass->el[j];
			pv.DisTwoway = pass->dis[j];
			pv.RCS = pass->rcs_all[j];
			Arc->m_ObsPos.push_back(pv);
		}
	}

	sprintf(Arc->m_FileName, "%s", pass->fileName);
	Arc->m_RCS = pass->RCS;
	Arc->m_numObs = static_cast<int>(Arc->m_ObsPos.size());
	Arc->m_Duration = (Arc->m_ObsPos[Arc->m_numObs - 1].JDUTC - Arc->m_ObsPos[0].JDUTC) * 86400.0;
	Arc->m_SiteIDs = mySite.m_IDs;
	Arc->m_SiteIDr = mySite.m_IDr;
	Arc->m_Site = mySite;
	Arc->m_pObsPos = &Arc->m_ObsPos;

	return true;
}

bool cRadarObsLL::pass2Arc_GroundOpt(vector<stSiteInfo> *m_Sites, obsPass1 *pass, RadarArc *Arc)
{
	stSiteInfo mySite;
	bool isSite = false;
	for (int i = 0; i < m_Sites->size(); i++)
	{
		if ((*m_Sites)[i].m_IDr == pass->equipIDr)
		{
			mySite = (*m_Sites)[i];
			isSite = true;
			pass->equipIDs = mySite.m_IDs;
			break;
		}
	}
	if (!isSite)
		return false;

	int nobs = static_cast<int>(pass->ra.size());
	double SiteTODPos[3];

	Arc->m_Type = pass->m_Type;
	for (int j = 0; j < nobs; j++)
	{
		epochPosVel1 pv;
		pv.JDUTC = pass->JD[j];

		// add by LL
		pv.Ra = pass->ra[j];
		pv.Dec = pass->dec[j];
		pv.Mag = pass->mag_all[j];

		// printf("%20.10f   %.3f   %.3f\n", pv.JDUTC, pv.Ra*g_dfRAD2DEG, pv.Dec * g_dfRAD2DEG);

		pv.SiteIDr = pass->equipIDr;
		pv.SiteIDs = pass->equipIDs;

		// add by Zhaogy
		ECEF2TOD(pv.JDUTC, mySite.m_siteECEFr, SiteTODPos);
		pv.satpos[0] = SiteTODPos[0];
		pv.satpos[1] = SiteTODPos[1];
		pv.satpos[2] = SiteTODPos[2];

		Arc->m_ObsPos.push_back(pv);
	}

	Arc->m_SiteIDs = pass->equipIDs;
	Arc->m_SiteIDr = pass->equipIDr;
	Arc->m_Site = mySite;

	sprintf(Arc->m_FileName, "%s", pass->fileName);
	Arc->m_RCS = pass->RCS;
	Arc->m_Mag = pass->Mag;
	Arc->m_numObs = static_cast<int>(Arc->m_ObsPos.size());
	Arc->m_StartJD = Arc->m_ObsPos[0].JDUTC;
	Arc->m_EndJD = Arc->m_ObsPos[Arc->m_numObs - 1].JDUTC;
	Arc->m_Duration = (Arc->m_ObsPos[Arc->m_numObs - 1].JDUTC - Arc->m_ObsPos[0].JDUTC) * 86400.0;
	Arc->m_pObsPos = &Arc->m_ObsPos;

	return true;
}

bool cRadarObsLL::pass2Arc_SpaceOpt(vector<stSiteInfo> *m_Sites, obsPass1 *pass, RadarArc *Arc)
{
	stSiteInfo mySite;

	epochPosVel1 pv;
	int nobs = static_cast<int>(pass->ra.size());

	Arc->m_Type = pass->m_Type;
	for (int j = 0; j < nobs; j++)
	{
		pv.JDUTC = pass->JD[j];

		// add by LL
		pv.Ra = pass->ra[j];
		pv.Dec = pass->dec[j];
		pv.Mag = pass->mag_all[j];

		pv.satpos[0] = pass->satx[j];
		pv.satpos[1] = pass->saty[j];
		pv.satpos[2] = pass->satz[j];

		// printf("%20.10f   %.3f   %.3f\n", pv.JDUTC, pv.Ra*g_dfRAD2DEG, pv.Dec * g_dfRAD2DEG);

		pv.SiteIDr = pass->equipIDr;
		pv.SiteIDs = pass->equipIDs;

		Arc->m_ObsPos.push_back(pv);
	}

	Arc->m_SiteIDs = pass->equipIDs;
	Arc->m_SiteIDr = pass->equipIDr;

	sprintf(Arc->m_FileName, "%s", pass->fileName);
	Arc->m_RCS = pass->RCS;
	Arc->m_Mag = pass->Mag;
	Arc->m_numObs = static_cast<int>(Arc->m_ObsPos.size());
	Arc->m_StartJD = Arc->m_ObsPos[0].JDUTC;
	Arc->m_EndJD = Arc->m_ObsPos[Arc->m_numObs - 1].JDUTC;
	Arc->m_Duration = (Arc->m_ObsPos[Arc->m_numObs - 1].JDUTC - Arc->m_ObsPos[0].JDUTC) * 86400.0;
	Arc->m_pObsPos = &Arc->m_ObsPos;

	return true;
}

bool cRadarObsLL::ReadSiteInfo(const char *SiteFile, std::vector<stSiteInfo> &Sites)
{
	ifstream fin;
	fin.open(SiteFile);
	stSiteInfo tSite;
	while (fin >> tSite.m_IDs >> tSite.m_siteBLHs[2] >> tSite.m_siteBLHs[1] >> tSite.m_siteBLHs[0] >> tSite.m_siteECEFs[0] >> tSite.m_siteECEFs[1] >> tSite.m_siteECEFs[2] >> tSite.m_IDr >> tSite.m_siteBLHr[2] >> tSite.m_siteBLHr[1] >> tSite.m_siteBLHr[0] >> tSite.m_siteECEFr[0] >> tSite.m_siteECEFr[1] >> tSite.m_siteECEFr[2])
	{
		// 判断测站类型
		if (tSite.m_IDr / 1000 == 9)
		{
			tSite.m_StationType = 3;
		}
		else
		{
			if (tSite.m_IDs != tSite.m_IDr)
			{
				tSite.m_StationType = 1;
			}
			else
			{
				tSite.m_StationType = 2;
			}
		}

		for (int i = 0; i < 2; i++)
		{
			tSite.m_siteBLHr[i] *= g_dfDEG2RAD;
			tSite.m_siteBLHs[i] *= g_dfDEG2RAD;
		}
		Sites.push_back(tSite);
	}
	fin.close();

	return true;
}

bool cRadarObsLL::ReadSBOSetInfo(const char *SBOFile, std::vector<stSiteInfo> &Sites)
{
	ifstream fin;
	fin.open(SBOFile);
	stSiteInfo tSite;
	string tline;
	while (getline(fin, tline))
	{
		stringstream ss(tline);
		int tID = 0;
		ss >> tID;

		bool isFind = false;
		int kSite = 0;
		for (int i = 0; i < Sites.size(); i++)
		{
			if (Sites[i].m_IDr == tID)
			{
				isFind = true;
				kSite = i;
				break;
			}
		}
		if (!isFind)
			continue;

		// 读入SBO参数
		ss >> Sites[kSite].FOVType >> Sites[kSite].opticalAxisAzEl[0] >> Sites[kSite].opticalAxisAzEl[1] >> Sites[kSite].FOVSize[0] >> Sites[kSite].m_VisualMagnitude >> Sites[kSite].m_phaseAngle >> Sites[kSite].m_SunCloseAngle >> Sites[kSite].m_MoonCloseAngle;

		Sites[kSite].opticalAxisAzEl[0] *= g_dfDEG2RAD;
		Sites[kSite].opticalAxisAzEl[1] *= g_dfDEG2RAD;
		Sites[kSite].FOVSize[0] *= g_dfDEG2RAD;
		Sites[kSite].m_phaseAngle *= g_dfDEG2RAD;
		Sites[kSite].m_SunCloseAngle *= g_dfDEG2RAD;
		Sites[kSite].m_MoonCloseAngle *= g_dfDEG2RAD;
	}
	fin.close();

	return true;
}

// 粗略转换ECEF2TOD
bool cRadarObsLL::ECEF2TOD(double JDUTC, double *ECEFPos, double *TODPos)
{
	double dfJDET, dfIntSecET, dfFraSecET, dfIntSecUTC, dfFraSecUTC;
	{
		g_DateTimeZ.ConversionJDSecond(JDUTC, dfIntSecUTC, dfFraSecUTC, true);
		if (!g_EOP.ConvertETUTC(dfIntSecET, dfFraSecET, dfIntSecUTC, dfFraSecUTC, false))
			return false;
		g_DateTimeZ.ConversionJDSecond(dfJDET, dfIntSecET, dfFraSecET, false);
	}

	double pdfMatrix[9];
	double dfIntJD = (int)dfJDET, dfFractionJD = dfJDET - (int)dfJDET;
	double dfEqn, dfJD = dfIntJD + dfFractionJD;

	cNPElement np;
	if (!np.GetEquinoxEquation(dfJD, dfEqn))
	{
		return false;
	}

	double dfIntJDUT, dfFraJDUT;
	if (!g_EOP.ConvertJDETUT(dfIntJD, dfFractionJD, dfIntJDUT, dfFraJDUT, true))
	{
		return false;
	}
	if (!cGreenwichST::ComputeGSTMatrixECEF2TOD(dfIntJDUT, dfFraJDUT, dfEqn, pdfMatrix))
	{
		return false;
	}
	g_Matrix.Mat_Multiply_Vec(pdfMatrix, ECEFPos, TODPos, 3, 3);

	return true;
}

bool cRadarObsLL::ReadArcFile(const char *FileName, vector<epochPosVel1> &ArcInfo)
{
	ifstream fin;
	fin.open(FileName);
	if (!fin.is_open())
		return false;

	string strType(FileName);
	strType = strType.substr(strType.find_last_of("/\\") + 1);
	strType = strType.substr(22, 1);

	epochPosVel1 tPV;
	if (strType == "7" || strType == "9")
	{
		// 光学.Obs
		while (fin >> tPV.JDUTC >> tPV.Ra >> tPV.Dec >> tPV.DisTwoway >> tPV.satpos[0] >> tPV.satpos[1] >> tPV.satpos[2])
		{
			tPV.Ra *= g_dfDEG2RAD;
			tPV.Dec *= g_dfDEG2RAD;
			ArcInfo.emplace_back(tPV);
		}
	}
	else
	{
		// 雷达.Obs
		while (fin >> tPV.JDUTC >> tPV.Az >> tPV.El >> tPV.DisTwoway >> tPV.pos[0] >> tPV.pos[1] >> tPV.pos[2])
		{
			tPV.Az *= g_dfDEG2RAD;
			tPV.El *= g_dfDEG2RAD;
			ArcInfo.emplace_back(tPV);
		}
	}
	fin.close();
	return true;
}

bool cRadarObsLL::ReadArcFile(const char *FileName, int ObsType, vector<epochPosVel1> &ArcInfo)
{
	ifstream fin;
	fin.open(FileName);
	if (!fin.is_open())
		return false;

	epochPosVel1 tPV;
	if (ObsType == 0)
	{
		// 雷达.Obs
		while (fin >> tPV.JDUTC >> tPV.Az >> tPV.El >> tPV.DisTwoway >> tPV.pos[0] >> tPV.pos[1] >> tPV.pos[2] >> tPV.RCS)
		{
			tPV.Az *= g_dfDEG2RAD;
			tPV.El *= g_dfDEG2RAD;
			ArcInfo.emplace_back(tPV);
		}
	}
	else
	{
		// 光学.OObs
		while (fin >> tPV.JDUTC >> tPV.Ra >> tPV.Dec >> tPV.DisTwoway >> tPV.satpos[0] >> tPV.satpos[1] >> tPV.satpos[2] >> tPV.Mag)
		{
			tPV.Ra *= g_dfDEG2RAD;
			tPV.Dec *= g_dfDEG2RAD;
			ArcInfo.emplace_back(tPV);
		}
	}
	fin.close();
	return true;
}

bool cRadarObsLL::ReadArcFile(const char *FileName, int ObsType, RadarArc &Arc)
{
	if (!ReadArcFile(FileName, ObsType, Arc.m_ObsPos))
		return false;

	string sub = string(FileName);
	sub = sub.substr(sub.find_last_of("/\\") + 1);
	sub = sub.substr(22, 4);
	Arc.m_SiteIDr = atoi(sub.c_str());
	Arc.m_numObs = static_cast<int>(Arc.m_ObsPos.size());

	if (Arc.m_numObs == 0)
		return false;

	sprintf(Arc.m_FileName, "%s", FileName);
	Arc.m_StartJD = Arc.m_ObsPos[0].JDUTC;
	Arc.m_EndJD = Arc.m_ObsPos[Arc.m_numObs - 1].JDUTC;
	Arc.m_Duration = (Arc.m_ObsPos[Arc.m_numObs - 1].JDUTC - Arc.m_ObsPos[0].JDUTC) * 86400.0;
	Arc.m_pObsPos = &Arc.m_ObsPos;

	// 获取RCS或Mag的中位数
	GetArcRCSMag(Arc);
	return true;
}

bool cRadarObsLL::GetArcRCSMag(RadarArc &Arc)
{
	if (Arc.m_pObsPos->size() == 0)
		return false;

	if (Arc.m_Type == 0)
	{
		vector<double> all_rcs;
		for (auto &eachpv : *Arc.m_pObsPos)
		{
			all_rcs.emplace_back(eachpv.RCS);
		}
		Arc.m_RCS = GetMedian(all_rcs);
	}
	else
	{
		vector<double> all_mag;
		for (auto &eachpv : *Arc.m_pObsPos)
		{
			all_mag.emplace_back(eachpv.Mag);
		}
		Arc.m_Mag = GetMedian(all_mag);
	}

	return true;
}

bool cRadarObsLL::ReadArcFile(const char *FileName, RadarArc &tArc)
{
	if (!ReadArcFile(FileName, tArc.m_ObsPos))
		return false;
	tArc.m_numObs = static_cast<int>(tArc.m_ObsPos.size());
	tArc.m_Duration = (tArc.m_ObsPos[tArc.m_numObs - 1].JDUTC - tArc.m_ObsPos[0].JDUTC) * 86400.0;
	return true;
}

bool cRadarObsLL::readObsFile(const char *fileName, vector<stSiteInfo> *m_Sites, RadarArc *tArc)
{
	// 判断文件类型
	string tName(fileName);
	string tType = tName.substr(tName.find_last_of(".") + 1);
	int Siter = atoi(tName.substr(tName.find_last_of("_") + 1, 4).c_str());
	if (tType == "Obs")
	{
		tArc->m_Type = 0;
		return ReadArcFile(fileName, tArc->m_Type, *tArc);
	}
	else if (tType == "OObs")
	{
		tArc->m_Type = 1;
		if (Siter / 1000 == 9)
			tArc->m_Type = 2;
		return ReadArcFile(fileName, tArc->m_Type, *tArc);
	}
	else
	{
		obsPass1 pass;
		if (!readObsFile(fileName, &pass))
			return false;

		//点数限制，若点数太多则均匀选点
		if (pass.numobs > 300)
		{
			UseLessPoint(pass, 300);
		}

		if (pass.m_Type == 0) // 0:LD   1:Ground based Optical  2:Space based Optical
			return pass2Arc(m_Sites, &pass, tArc);
		else if (pass.m_Type == 1)
			return pass2Arc_GroundOpt(m_Sites, &pass, tArc);
		else
			return pass2Arc_SpaceOpt(m_Sites, &pass, tArc);
	}
	return true;
}

bool cRadarObsLL::UseLessPoint(obsPass1& ObsPass, int nPoint)
{
	// 计算观测持续时间和临时点数
	ObsPass.numtemp = ObsPass.numobs;
	ObsPass.duration = (ObsPass.JD[ObsPass.numobs - 1] - ObsPass.JD[0]) * 86400.0;

	// 如果原始点数 <= 目标点数，无需缩减
	if (ObsPass.numobs <= nPoint)
		return true;

	// 保存原始数据副本
	obsPass1 pass = ObsPass;

	// 清空原始观测数据向量
	ObsPass.JD.clear();
	ObsPass.az.clear();
	ObsPass.el.clear();
	ObsPass.dis.clear();
	ObsPass.ra.clear();
	ObsPass.dec.clear();
	ObsPass.satx.clear();
	ObsPass.saty.clear();
	ObsPass.satz.clear();

	// 释放内存（C++11 后 shrink_to_fit 更合适，但保留原风格）
	std::vector<double>().swap(ObsPass.JD);
	std::vector<double>().swap(ObsPass.az);
	std::vector<double>().swap(ObsPass.el);
	std::vector<double>().swap(ObsPass.dis);
	std::vector<double>().swap(ObsPass.ra);
	std::vector<double>().swap(ObsPass.dec);
	std::vector<double>().swap(ObsPass.satx);
	std::vector<double>().swap(ObsPass.saty);
	std::vector<double>().swap(ObsPass.satz);

	// 计算采样步长（浮点型确保均匀分布）
	double step = static_cast<double>(pass.numobs - 1) / (nPoint - 1);

	// 根据观测类型抽取数据
	if (ObsPass.m_Type == 0) {  // 方位角/俯仰角类型
		for (int i = 0; i < nPoint; ++i) {
			// 计算索引（四舍五入取整）
			int idx = static_cast<int>(i * step + 0.5);
			// 边界保护
			if (idx < 0) idx = 0;
			if (idx >= pass.numobs) idx = pass.numobs - 1;

			ObsPass.JD.push_back(pass.JD[idx]);
			ObsPass.az.push_back(pass.az[idx]);
			ObsPass.el.push_back(pass.el[idx]);
			ObsPass.dis.push_back(pass.dis[idx]);
			ObsPass.satx.push_back(pass.satx[idx]);
			ObsPass.saty.push_back(pass.saty[idx]);
			ObsPass.satz.push_back(pass.satz[idx]);
		}
		ObsPass.numobs = static_cast<int>(ObsPass.dis.size());
	}
	else {  // 赤经/赤纬类型
		for (int i = 0; i < nPoint; ++i) {
			int idx = static_cast<int>(i * step + 0.5);
			if (idx < 0) idx = 0;
			if (idx >= pass.numobs) idx = pass.numobs - 1;

			ObsPass.JD.push_back(pass.JD[idx]);
			ObsPass.ra.push_back(pass.ra[idx]);
			ObsPass.dec.push_back(pass.dec[idx]);
			ObsPass.satx.push_back(pass.satx[idx]);
			ObsPass.saty.push_back(pass.saty[idx]);
			ObsPass.satz.push_back(pass.satz[idx]);
		}
		ObsPass.numobs = static_cast<int>(ObsPass.ra.size());
	}

	return true;
}

bool cRadarObsLL::readObsFromMap(string fileName, RadarArc *tArc,
								 map<string, vector<epochPosVel1>> *ArcObsMap)
{
	map<string, vector<epochPosVel1>>::iterator ObsIter;
	string tName(fileName);
	tName = tName.substr(tName.find_last_of("/\\") + 1);
	tName = tName.substr(0, tName.find_last_of("."));
	ObsIter = ArcObsMap->find(tName);
	if (ObsIter == ArcObsMap->end())
		return false;
	tArc->m_pObsPos = &ObsIter->second;

	sprintf(tArc->m_FileName, "%s", fileName.c_str());
	string sub = string(tName);
	sub = sub.substr(sub.find_last_of("_") + 1, 4);
	tArc->m_SiteIDr = atoi(sub.c_str());

	if (tArc->m_SiteIDr / 1000 == 1)
		tArc->m_Type = 0;
	else if (tArc->m_SiteIDr / 1000 == 9)
		tArc->m_Type = 2;
	else
		tArc->m_Type = 1;

	tArc->m_numObs = static_cast<int>(tArc->m_pObsPos->size());
	tArc->m_StartJD = (*tArc->m_pObsPos)[0].JDUTC;
	tArc->m_EndJD = (*tArc->m_pObsPos)[tArc->m_numObs - 1].JDUTC;
	tArc->m_Duration = ((*tArc->m_pObsPos)[tArc->m_numObs - 1].JDUTC - (*tArc->m_pObsPos)[0].JDUTC) * 86400.0;
	GetArcMag(tArc->m_pObsPos, tArc->m_Mag);
	return true;
}

bool cRadarObsLL::GetArcMag(std::vector<epochPosVel1> *pObsPos, double &dMag)
{
	std::vector<double> MagVec;
	for (auto &eachObs : *pObsPos)
	{
		MagVec.emplace_back(eachObs.Mag);
	}
	if (MagVec.size() == 0)
		return false;

	cRadarObsLL RO;
	dMag = RO.GetMedian(MagVec);
	return true;
}

int cRadarObsLL::JudgeObsType(string ObsPath)
{
	string tName = ObsPath.substr(ObsPath.find_last_of("/\\") + 1);
	tName = tName.substr(tName.find_last_of("_") + 1, 4);
	int SiteID = atoi(tName.c_str());
	if (SiteID / 1000 == 1)
		return 0;
	else if (SiteID / 1000 == 9)
		return 2;
	else
		return 1;

	return 0;
}

bool cRadarObsLL::LDOutlierDetection(int n, double *JD, double *x, bool *mark)
{
	int i;
	for (i = 0; i < n; i++)
		mark[i] = true;

	cPolyFit PF;
	double *t = new double[n];
	for (i = 0; i < n; i++)
		t[i] = (JD[i] - JD[0]) * 86400.0;

	PF.SetObsAndTime(x, t, n, 2);
	PF.Fitting(mark);

	return true;
}

double cRadarObsLL::GetMedian(vector<double> &Numbers)
{
	int n = static_cast<int>(Numbers.size());
	if (n == 0)
		return 0.0;

	sort(Numbers.begin(), Numbers.end());

	if (n % 2 != 0)
	{
		return Numbers[n / 2];
	}
	else
	{
		return (Numbers[n / 2] + Numbers[n / 2 - 1]) / 2;
	}
}