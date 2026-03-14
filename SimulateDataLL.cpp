#include "SimulateDataLL.h"
#include "Constant.h"
#include "COORTRAN.H"
#include <Eigen/Dense>

using namespace Eigen;


//  gyzhao   20200713
//参考系：接收站B的地平坐标系（右手系）；
//输入：接收站B观测的Azimuth、Elevation, 发射站A在参考系下的三维坐标(以收站为原点，计算发站坐标)、两站间收发总距离；
//输出：空间目标在接收站地平坐标系中的三维坐标
//单位：位置m,  角度radian 
bool cSimulateDataLL::NewtonRaphson_SR_RadarPositioning(double* SendStationECEF, double* ReceiveStationECEF, double* ReceiveStationBLH,
	double Azimuth, double Elevation, double Rou, double *satpos)
{
	cCoorTrans CT;
	double deltaXYZ[3], NEUs[3];

	deltaXYZ[0] = SendStationECEF[0] - ReceiveStationECEF[0];
	deltaXYZ[1] = SendStationECEF[1] - ReceiveStationECEF[1];
	deltaXYZ[2] = SendStationECEF[2] - ReceiveStationECEF[2];
	
	CT.ConvertNEU2DXYZ(deltaXYZ, NEUs, ReceiveStationBLH[0], ReceiveStationBLH[1], 0);
	NEUs[1] *= -1;

	if (Elevation == g_dfHALFPI && Elevation == -g_dfHALFPI) { return 0; }

	Vector3d X;
	Vector3d Fx;
	Vector3d staB;

	staB(0) = NEUs[0];
	staB(1) = NEUs[1];
	staB(2) = NEUs[2];

	//Jacobi矩阵：
	Matrix3d MatJacobi;

	//计算目标初值，注意方位角和右手系旋角相反
	X(0) = 0.5*Rou*cos(Elevation)*cos(-Azimuth);
	X(1) = 0.5*Rou*cos(Elevation)*sin(-Azimuth);
	X(2) = 0.5*Rou*sin(Elevation);

	bool BreakNow = 0;
	Vector3d lastX = X;
	int iterate_times = 0;
	double DELTA = 1;//精度限定

	 //Newton-Raphson迭代
	do
	{
		iterate_times++;
		if (iterate_times > 1000)
		{
			return 0;
		}

		//构造Jacobi阵
		MatJacobi(0, 0) = X(0) / sqrt(X(0)*X(0) + X(1)*X(1) + X(2)*X(2))
			+ (X(0) - staB(0)) / sqrt((X(0) - staB(0))*(X(0) - staB(0)) + (X(1) - staB(1))*(X(1) - staB(1)) + (X(2) - staB(2))*(X(2) - staB(2)));

		MatJacobi(0, 1) = X(1) / sqrt(X(0)*X(0) + X(1)*X(1) + X(2)*X(2))
			+ (X(1) - staB(1)) / sqrt((X(0) - staB(0))*(X(0) - staB(0)) + (X(1) - staB(1))*(X(1) - staB(1)) + (X(2) - staB(2))*(X(2) - staB(2)));

		MatJacobi(0, 2) = X(2) / sqrt(X(0)*X(0) + X(1)*X(1) + X(2)*X(2))
			+ (X(2) - staB(2)) / sqrt((X(0) - staB(0))*(X(0) - staB(0)) + (X(1) - staB(1))*(X(1) - staB(1)) + (X(2) - staB(2))*(X(2) - staB(2)));

		MatJacobi(1, 0) = 1 / sqrt(1 - X(2)*X(2) / (X(0)*X(0) + X(1)*X(1) + X(2)*X(2)))
			*(-1)*X(0)*X(2)*pow(X(0)*X(0) + X(1)*X(1) + X(2)*X(2), -3 / 2);

		MatJacobi(1, 1) = 1 / sqrt(1 - X(2)*X(2) / (X(0)*X(0) + X(1)*X(1) + X(2)*X(2)))
			*(-1)*X(1)*X(2)*pow(X(0)*X(0) + X(1)*X(1) + X(2)*X(2), -3 / 2);

		MatJacobi(1, 2) = 1 / sqrt(1 - X(2)*X(2) / (X(0)*X(0) + X(1)*X(1) + X(2)*X(2)))
			*(pow(X(0)*X(0) + X(1)*X(1) + X(2)*X(2), -1 / 2) - X(2)*X(2)*pow(X(0)*X(0) + X(1)*X(1) + X(2)*X(2), -3 / 2));


		MatJacobi(2, 0) = -X(1) / (X(0)*X(0) + X(1)*X(1));

		MatJacobi(2, 1) = X(0) / (X(0)*X(0) + X(1)*X(1));

		MatJacobi(2, 2) = 0;

		if (g_dfPI <= Azimuth && Azimuth < g_dfTWOPI)
		{
			MatJacobi(2, 0) *= -1;

			MatJacobi(2, 1) *= -1;
		}


		//计算Fx向量
		if (0 <= Azimuth && Azimuth < g_dfPI)
		{
			Fx(0) = sqrt(X(0)*X(0) + X(1)*X(1) + X(2)*X(2))
				+ sqrt((X(0) - staB(0))*(X(0) - staB(0)) + (X(1) - staB(1))*(X(1) - staB(1)) + (X(2) - staB(2))*(X(2) - staB(2))) - Rou;
			Fx(1) = asin(X(2) / sqrt(X(0)*X(0) + X(1)*X(1) + X(2)*X(2))) - Elevation;
			Fx(2) = acos(X(0) / sqrt(X(0)*X(0) + X(1)*X(1))) - Azimuth;

		}
		else
		{
			Fx(0) = sqrt(X(0)*X(0) + X(1)*X(1) + X(2)*X(2))
				+ sqrt((X(0) - staB(0))*(X(0) - staB(0)) + (X(1) - staB(1))*(X(1) - staB(1)) + (X(2) - staB(2))*(X(2) - staB(2))) - Rou;
			Fx(1) = asin(X(2) / sqrt(X(0)*X(0) + X(1)*X(1) + X(2)*X(2))) - Elevation;

			Fx(2) = g_dfTWOPI - Azimuth - acos(X(0) / sqrt(X(0)*X(0) + X(1)*X(1)));

		}

		//计算新的x y z 
		X = X - MatJacobi.inverse()*Fx;

		//有效数字核验 
		if (abs(X(0) - lastX(0)) < DELTA && abs(X(1) - lastX(1)) < DELTA && abs(X(2) - lastX(2)) < DELTA)
		{
			BreakNow = 1;
		}

		lastX = X;

	} while (BreakNow == 0);

	double NEUr[3], deltaXYZr[3];
	NEUr[0] = X(0);
	NEUr[1] = -X(1);
	NEUr[2] = X(2);

	CT.ConvertNEU2DXYZ(deltaXYZr, NEUr, ReceiveStationBLH[0], ReceiveStationBLH[1], 1);

	satpos[0] = deltaXYZr[0] + ReceiveStationECEF[0];
	satpos[1] = deltaXYZr[1] + ReceiveStationECEF[1];
	satpos[2] = deltaXYZr[2] + ReceiveStationECEF[2];

	return 1;
}


//  LL   20201118
//参考系：接收站B的地平坐标系（右手系）；
//输入：接收站B观测的Azimuth、Elevation, 发射站A在参考系下的三维坐标(以收站为原点，计算发站坐标)、单站双程距离；
//输出：空间目标在接收站地平坐标系中的三维坐标
//单位：位置m,  角度radian 
bool cSimulateDataLL::SingleEquip_RadarPositioning(double* StationECEF, double* StationBLH,
	double Azimuth, double Elevation, double Rou, double *satpos)
{
	cCoorTrans CT;



	if (Elevation == g_dfHALFPI && Elevation == -g_dfHALFPI) { return 0; }

	Vector3d X;
	//计算目标初值，注意方位角和右手系旋角相反
	X(0) = 0.5*Rou*cos(Elevation)*cos(-Azimuth);
	X(1) = 0.5*Rou*cos(Elevation)*sin(-Azimuth);
	X(2) = 0.5*Rou*sin(Elevation);


	double NEUr[3], deltaXYZr[3];
	NEUr[0] = X(0);
	NEUr[1] = -X(1);
	NEUr[2] = X(2);

	CT.ConvertNEU2DXYZ(deltaXYZr, NEUr, StationBLH[0], StationBLH[1], 1);

	satpos[0] = deltaXYZr[0] + StationECEF[0];
	satpos[1] = deltaXYZr[1] + StationECEF[1];
	satpos[2] = deltaXYZr[2] + StationECEF[2];

	return 1;
}

