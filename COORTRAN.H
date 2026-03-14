/******************************************************************************

  transformation between geodetic and cartesian coordinate systems 


******************************************************************************/

#ifndef COORTRAN_H
#define COORTRAN_H

#include <math.h>
#include "Matrix.h"


class cCoorTrans
{
public:
	cCoorTrans();
	~cCoorTrans();
	void XYZ_to_BLH( double x , double y , double z ,
					 double &b, double &l, double &h,
					 double semi_major = 6378136.5,
					 double flatening_denominator = 298.25642 );
	void BLH_to_XYZ( double b,  double l,  double h,
					 double &x, double &y, double &z,
					 double semi_major = 6378136.5,
					 double flatening_denominator = 298.25642 );
	void ConvertXYZ2BLH( double &x, double &y, double &z ,
					     double &b, double &l, double &h,
						 bool bXYZ2BLH,
					     double semi_major = 6378136.5,
					     double eccentricity = 0.006694395747586 );
	void trans_matrix_XYZ_to_BLH( double lati, double longi,
								  double *trans_matrix );
	void deltaXYZ_to_deltaBLH( double dx, double dy, double dz,
							   double &db,double &dl, double &dh,
							   double lati, double longi );
	void deltaBLH_to_deltaXYZ( double db, double dl, double dh,
							   double &dx,double &dy,double &dz,
							   double lati, double longi );
	void covBLH_to_covXYZ( double *cov_blh, double *cov_xyz,
						   double lati, double longi );
	void covXYZ_to_covBLH( double *cov_xyz, double *cov_blh,
						   double lati, double longi );
	void meridian_length_coefficient( double *coefficient,
		                              double semi_major = 6378136.5,
					                  double flatening_denominator = 298.25642 );
	void latilongi_to_northeast( double lati, double longi,
								 double &north, double &east,
								 double central_longi,
								 double semi_major = 6378136.5,
					             double flatening_denominator = 298.25642 );
	void northeast_to_loatlongi( double north, double east,
								 double &lati, double &longi,
								 double central_longi,
								 double semi_major = 6378136.5,
					             double flatening_denominator = 298.25642 );
	void geo_azimuth_distance( double latifrom, double longifrom,
							   double latoto,   double longito,
							   double &azimuth, double &distance,
							   double semi_major = 6378136.5,
					           double flatening_denominator = 298.25642 );
	void horistow( double delta_azimuth,
				   double delta_lati, double delta_longi,
				   double lati_begin, double longi_begin,
				   double &lati_end,  double &longi_end,
				   double semi_major = 6378136.5,
				   double flatening_denominator = 298.25642 );
	void DDMMSS2RAD( double &deg );
	void RAD2DDMMSS( double &rad );
	void DegreeFraction( double deg, int &Deg, int &Min, int &Sec, int &IntDec );
	void DeltaXYZ_to_EleAzi( double dX,  double dY, double dZ,
							 double Lat, double Longi,
							 double &elevation, double &azimuth, double &distance );
	void EleAziDis_to_DeltaXYZ( double azimuth, double elevation, double distance,
							    double Lat, double Longi,
							    double &dX, double &dY, double &dZ );
	void DeltaXYZ_to_EleAzi( double dX, double dY, double dZ,
							 double dfSinLat, double dfCosLat,
							 double dfSinLong, double dfCosLong,
							 double &elevation, double &azimuth, double &distance );
	void DeltaXYZ_to_RAD( double dX,  double dY, double dZ,
						  double &RA, double &D, double &distance );
	void RotatingAroundXAxis( double RotationAngle, 
							  double &X, double &Y, double &Z );
	void RotatingAroundYAxis( double RotationAngle, 
							  double &X, double &Y, double &Z );
	void RotatingAroundZAxis( double RotationAngle, 
							  double &X, double &Y, double &Z );
	void dXdYdZ_to_DisAziZenith( double Lat, double Longi,
								 double dX, double dY, double dZ,
								 double &Distance, double &Azimuth, double &Zenith );
	void ComputedXdYdZFromAzimuthElevation( double Azimuth, 
										    double Elevation, 
										    double &dx, 
										    double &dy, 
										    double &dz );
	void ComputeAzimuthElevationFromdXdYdZ( double &Azimuth, 
										    double &Elevation, 
								            double dx, 
								            double dy, 
								            double dz );

	void ComputeECEFUnitVector( double dfLat, double dfLong, 
		                        double *pdfEast, double *pdfNorth, double *pdfUp );

	void ComputeRotationMatrix( int nAxis, double dfAngle, double *pdfRotationMatrix );

	void ConvertNEU2DXYZ( double &dX, double &dY, double &dZ,
		                  double &North, double &East, double &Up,
						  double latitude, double longitude, bool bNEU2DXYZ );
	void ConvertNEU2DXYZ( double *pdfDeltaXYZ, double *pdfNEU, 
						  double latitude, double longitude, bool bNEU2DXYZ );

	void ConvertPosDif2OrbitBias( double *pos, double *vel, double *posDif, double *bias );
	void ConvertPosDif2OrbitBias( double *pos, double *vel, double *posDif, double *bias, double &flight_angle );
	void ConvertOrbitBias2PosDif( double *pos, double *vel, double *bias, double *posDis );

	void ConvertDeltaAzElDis2DeltaNEU( double az, double el, double dis, 
		                               double delta_az, double delta_el, double delta_dis, 
				 					   double &delta_north, double &delta_east, double &delta_up );

	void convertAzElBias2AlongCrossBias( double deltaAz, double deltaEl, double El, double azDot, double elDot,
		                                 double &along, double &cross );
	void convertAlongCrossBias2AzElBias( double along, double cross, double El, double azDot, double elDot,
		                                 double &deltaAz, double &deltaEl );

};

#endif