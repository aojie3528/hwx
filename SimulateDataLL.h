#ifndef SIMULATE_DATA_LL_H
#define SIMULATE_DATA_LL_H

class cSimulateDataLL
{
public:
	bool NewtonRaphson_SR_RadarPositioning(double* SendStationECEF, double* ReceiveStationECEF, double* ReceiveStationBLH,
		double Azimuth, double Elevation, double Rou, double *satpos);

	bool SingleEquip_RadarPositioning(double* StationECEF, double* StationBLH,
		double Azimuth, double Elevation, double Rou, double *satpos);
};

#endif