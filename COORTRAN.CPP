/******************************************************************************

  transformation between geodetic and cartesian coordinate systems 


******************************************************************************/

#include "Constant.h"
#include "CoorTran.h"
#include "matrix.h"

cCoorTrans::cCoorTrans()
{ 
 ;
}

cCoorTrans::~cCoorTrans()
{
 ;
}

/******************************************************************************

  Convert from rectangular X, Y, Z to geodetic latitude, longitude and height

******************************************************************************/

void cCoorTrans::XYZ_to_BLH( double x, double y, double z,
							 double &b,double &l,double &h,
							 double semi_major,
							 double flatening_denominator )
{
    double e = 1.0 / flatening_denominator ;   
	e = 2 * e - e * e;

	ConvertXYZ2BLH( x, y, z, b, l, h, true, semi_major, e );
}

/******************************************************************************

  Convert from geodetic latitude, longitude and height to rectangular X, Y, Z

******************************************************************************/

void cCoorTrans::BLH_to_XYZ( double b, double l, double h,
							 double &x,double &y,double &z,
							 double semi_major,
							 double flatening_denominator )
{
    double e = 1.0 / flatening_denominator ;   
	e = 2 * e - e * e;

	ConvertXYZ2BLH( x, y, z, b, l, h, false, semi_major, e );    
}

/******************************************************************************

  Compute the transformation matrix from rectangular delta XYZ to 
  delta latitude, longitude and height (north, east and up)

******************************************************************************/

void cCoorTrans::trans_matrix_XYZ_to_BLH( double lati,
										  double longi,
										  double *trans_matrix )
		  /*
		    dB(m)                      dX
		    dL(m)   =   trans_matrix * dY
		    dH(m)                      dZ
		  */
{
    trans_matrix[ 0 ] = - sin ( lati ) * cos ( longi ) ;
    trans_matrix[ 1 ] = - sin ( lati ) * sin ( longi ) ;
    trans_matrix[ 2 ] =   cos ( lati )                 ;
    trans_matrix[ 3 ] = - sin ( longi )                ;
    trans_matrix[ 4 ] =   cos ( longi )                ;
    trans_matrix[ 5 ] =   0.0                          ;
    trans_matrix[ 6 ] =   cos ( lati ) * cos ( longi ) ;
    trans_matrix[ 7 ] =   cos ( lati ) * sin ( longi ) ;
    trans_matrix[ 8 ] =   sin ( lati )                 ;
}

/******************************************************************************

  Convert from rectangular delta X,Y,Z and delta latitide, longitude and height

******************************************************************************/

void cCoorTrans::deltaXYZ_to_deltaBLH( double dx, double dy, double dz,
									   double &db,double &dl,double &dh,
									   double lati, double longi )
{
    double trans_matrix[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	
    trans_matrix_XYZ_to_BLH( lati, longi, trans_matrix );

    db =  trans_matrix[ 0 ] * dx + trans_matrix[ 1 ] * dy
		+ trans_matrix[ 2 ] * dz;
    dl =  trans_matrix[ 3 ] * dx + trans_matrix[ 4 ] * dy
		+ trans_matrix[ 5 ] * dz;
    dh =  trans_matrix[ 6 ] * dx + trans_matrix[ 7 ] * dy
		+ trans_matrix[ 8 ] * dz;
}

/******************************************************************************

  Convert from delta latitide, longitude and height to rectangular delta X,Y,Z

******************************************************************************/

void cCoorTrans::deltaBLH_to_deltaXYZ( double db, double dl, double dh,
									   double &dx,double &dy,double &dz,
									   double lati, double longi )
{
    double trans_matrix[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	
    trans_matrix_XYZ_to_BLH( lati, longi, trans_matrix );

    dx =  trans_matrix[ 0 ] * db + trans_matrix[ 3 ] * dl
		+ trans_matrix[ 6 ] * dh;
    dy =  trans_matrix[ 1 ] * db + trans_matrix[ 4 ] * dl
		+ trans_matrix[ 7 ] * dh;
    dz =  trans_matrix[ 2 ] * db + trans_matrix[ 5 ] * dl
		+ trans_matrix[ 8 ] * dh;
}

/******************************************************************************

  Convert the covariance matrix from delta latitude/longitude/height to 
  rectabgular delta x, y, z

******************************************************************************/

void cCoorTrans::covBLH_to_covXYZ( double *cov_blh,  double *cov_xyz,
								   double lati,		double longi )
{
    double trans_matrix[ 9 ],temp[ 9 ], temp1[ 9 ];
    int m = 3;

	g_Matrix.Set_Vec_Zero( temp, 9 );
	g_Matrix.Set_Vec_Zero( temp1, 9 );
	g_Matrix.Set_Vec_Zero( trans_matrix, 9 );

    trans_matrix_XYZ_to_BLH( lati, longi, trans_matrix );

    g_Matrix.TriMat_Multiply_Mat( cov_blh, trans_matrix, temp, m, m );

    g_Matrix.Mat_Transpose( trans_matrix, m, m );

    g_Matrix.Mat_Multiply_Mat( trans_matrix, temp, temp1, m, m, m );

    cov_xyz[ 0 ] = temp1[ 0 ];         
	cov_xyz[ 1 ] = temp1[ 3 ]; 
	cov_xyz[ 2 ] = temp1[ 4 ];
    cov_xyz[ 3 ] = temp1[ 6 ]; 
	cov_xyz[ 4 ] = temp1[ 7 ]; 
	cov_xyz[ 5 ] = temp1[ 8 ];
}

/******************************************************************************

  Convert the covariance matrix from rectabgular delta x, y, z to
  delta latitude/longitude/height

******************************************************************************/

void cCoorTrans::covXYZ_to_covBLH( double *cov_xyz,  double *cov_blh,
								   double lati,		double longi )
{
    double trans_matrix[ 9 ], temp[ 9 ], temp1[ 9 ];
    int m = 3;
	
	g_Matrix.Set_Vec_Zero( temp, 9 );
	g_Matrix.Set_Vec_Zero( temp1, 9 );
	g_Matrix.Set_Vec_Zero( trans_matrix, 9 );

    trans_matrix_XYZ_to_BLH( lati, longi, trans_matrix );
	
    g_Matrix.Mat_Multiply_TriMat( trans_matrix, cov_xyz, temp, m, m );

    g_Matrix.Mat_Transpose( trans_matrix, m, m );

    g_Matrix.Mat_Multiply_Mat( temp, trans_matrix, temp1, m, m, m );

    cov_blh[ 0 ] = temp1[ 0 ];         
	cov_blh[ 1 ] = temp1[ 3 ]; 
	cov_blh[ 2 ] = temp1[ 4 ];
    cov_blh[ 3 ] = temp1[ 6 ]; 
	cov_blh[ 4 ] = temp1[ 7 ]; 
	cov_blh[ 5 ] = temp1[ 8 ];
}

/******************************************************************************

  Compute the coefficients of formula of calculating meridian length

******************************************************************************/

void cCoorTrans::meridian_length_coefficient( double *coefficient, 
											  double semi_major,
											  double flatening_denominator )
{
    double e = 0, ratio = 0;

//  e_square 
    e = 1.0 / flatening_denominator; 
	e = 2.0 * e - e * e;

    ratio = semi_major * ( 1.0 - e );

    coefficient[ 0 ] = (   pow( e, 5.0 ) * 43659.0 / 65536.0
		                 + pow( e, 4.0 ) * 11025.0 / 16384.0
		                 + pow( e, 3.0 ) * 175.0 / 256.0
		                 + pow( e, 2.0 ) * 45.0 / 64.0
		                 + e * 0.75                     + 1.0 ) * ratio;
    coefficient[ 1 ] =-(   pow( e, 5.0 ) * 72765.0 / 65536.0
		                 + pow( e, 4.0 ) * 2205.0 / 2048.0
		                 + pow( e, 3.0 ) * 525.0 / 512.0
		                 + pow( e, 2.0 ) * 15.0 / 16.0
					     + e * 0.75                           ) * ratio / 2.0;
    coefficient[ 2 ] = (   pow( e, 5.0 ) * 10395.0 / 16384.0
		                 + pow( e, 4.0 ) * 2205.0 / 4096.0
		                 + pow( e, 3.0 ) * 105.0 / 256.0
		                 + pow( e, 2.0 ) * 15.0 / 64.0        ) * ratio / 4.0;
    coefficient[ 3 ] =-(   pow( e, 5.0 ) * 31185.0 / 131072.0
		                 + pow( e, 4.0 ) * 315.0 / 2048.0
		                 + pow( e, 3.0 ) * 35.0 / 512.0       ) * ratio / 6.0;
    coefficient[ 4 ] = (   pow( e, 5.0 ) * 3465.0 / 65536.0
		                 + pow( e, 4.0 ) * 315.0 / 16384.0    ) * ratio / 8.0;
    coefficient[ 5 ] = -pow( e, 5.0 ) * 693.0 / 131072.0        * ratio / 10.0;
}

/******************************************************************************

  gauss-kruge projection, latitude/longitude to north/east 

******************************************************************************/

void cCoorTrans::latilongi_to_northeast( double lati, double longi,
										 double &north, double &east,
										 double central_longi,
										 double semi_major,
										 double flatening_denominator )
{
    double e = 0, e1 = 0, coeff[ 6 ] = { 0, 0, 0, 0, 0, 0 }, dl = 0;
    double N = 0, s = 0, c = 0, t = 0;

//  e_square 
    e = 1.0 / flatening_denominator; 
	e = 2.0 * e - e * e;
    meridian_length_coefficient( coeff, semi_major,
								 flatening_denominator );

//  delta_longitude 
	dl = longi - central_longi;

	north = coeff[ 0 ] * lati + coeff[ 1 ] * sin(  2.0 * lati )
					          + coeff[ 2 ] * sin(  4.0 * lati )
						      + coeff[ 3 ] * sin(  6.0 * lati )
						      + coeff[ 4 ] * sin(  8.0 * lati )
						      + coeff[ 5 ] * sin( 10.0 * lati );

    s = sin ( lati );
    c = cos ( lati );
    t = tan ( lati );

//  radius of prime 
    N = semi_major / sqrt( 1.0 - e * s * s );

//  eta_square 
    e1 = e / ( 1.0 - e ) * c * c;

    north = north + s * c *             N * pow( dl,2.0 ) / 2.0
		          + s * pow( c, 3.0 ) * N * pow( dl,4.0 ) / 24.0
		              * ( 5.0 - t * t + 9.0 * e1 + 4.0 * e1 * e1 )
		          + s * pow( c, 5.0 ) * N * pow( dl,6.0 ) / 720.0
		              * ( 61.0 - 58.0 * t * t - pow( t, 4.0 ) );

    east = 500000.0 + N * c *                  dl
		            + N * pow( c, 3.0 ) * pow( dl, 3.0 ) / 6.0
		                * ( 1.0 - t * t + e1 )
		            + N * pow( c, 5.0 ) * pow( dl, 5.0 ) / 120.0
		                * ( 5. - 18. * t * t + pow( t, 4.0 ) + 14.0 * e1 - 58. * e1 * t * t );
}

/******************************************************************************

  gauss-kruge projection, north/east to latitude/longitude 

******************************************************************************/

void cCoorTrans::northeast_to_loatlongi( double north, double east,
										 double &lati, double &longi,
										 double central_longi,
										 double semi_major,
										 double flatening_denominator )
{
    double e = 0, e1 = 0, coeff[ 6 ] = { 0, 0, 0, 0, 0, 0 }, dl = 0;
    double N = 0, M = 0, c = 0, t = 0, east1 = 0, lat0 = 0, lat = 0;

//  e_square 
    e = 1.0 / flatening_denominator; 
	e = 2.0 * e - e * e;
    meridian_length_coefficient( coeff, semi_major,
								 flatening_denominator );

    lat0 = north / coeff[ 0 ];
RET:
    lat = - coeff[ 1 ] * sin(  2.0 * lat0 ) - coeff[ 2 ] * sin( 4.0 * lat0 )
	      - coeff[ 3 ] * sin(  6.0 * lat0 ) - coeff[ 4 ] * sin( 8.0 * lat0 )
	      - coeff[ 5 ] * sin( 10.0 * lat0 );
    lat = ( lat + north ) / coeff[ 0 ];
    if( fabs( lat - lat0 ) < 1.0e-10 ) goto OOUT;
    lat0 = lat;
    goto RET;

OOUT:
//  eta_square 
    e1 = e / ( 1.0 - e ) * cos( lat ) * cos( lat );

//  radius of prime 
    N = semi_major / sqrt( 1.0 - e * sin( lat ) * sin( lat ) );

//  radius of meridian 
    M = N * ( 1.0 - e ) / ( 1.0 - e * sin( lat ) * sin( lat ) );

	t = tan( lat );

	east1 = east - 500000.0;

    lat0 = - t / M * pow( east1, 2.0 ) / N             / 2.0
	       + t / M * pow( east1, 4.0 ) / pow( N, 3.0 ) / 24.0
	               * ( 5.0 + 3.0 * t * t + e1 - 9.0 * e1 * t * t )
	       - t / M * pow( east1, 6.0 ) / pow( N, 5.0 ) / 720.0
	               * ( 61.0 + 90.0 * t * t + 45.0 * pow( t, 4.0 ) );
    lati = lat0 + lat;

    c = 1.0 / cos( lat );
    east1 = east1 / N;
    dl =   c * east1
	     - c * pow( east1, 3.0 ) / 6.0
	         * (1.0+2.0*t*t+e1)
	     + c * pow( east1, 5.0 ) / 120.0
	         * ( 5.0 + 28. * t * t + 24. * pow( t, 4.0 ) + 6. * e1 + 8. * e1 * t * t );
    longi = central_longi + dl;
}

/******************************************************************************

  compute the geodetic azimuth and distance between two points of given 
  geodetic latitudes and longitudes on a given ellipsoid

******************************************************************************/

void cCoorTrans::geo_azimuth_distance( double latifrom, double longifrom,
									   double latito,   double longito,
									   double &azimuth, double &distance,
									   double semi_major,
									   double flatening_denominator )
{
    double e = 0, e1 = 0;
    double N = 0, c = 0, t = 0, a1 = 0, a2 = 0, db = 0, dl = 0, alfa = 0;

    t = tan( latifrom );
    c = cos( latifrom );
    e = 1.0 / flatening_denominator; 
	e = 2.0 * e - e * e;
    e1 = e / ( 1.0 - e ) * c * c;
    N = semi_major / sqrt( 1.0 - e * sin( latifrom ) * sin( latifrom ) );

    double b10 = N * ( 1.0 - e1 + e1 * e1 - pow( e1,3.0 ) );
    double b20 = N * t * ( 3.0 * e1 - 6.0 * e1 * e1 ) / 2.0;
    double b02 = N * c * c * t / 2.0;
    double b30 = N * ( e1 - t * t * e1 - 2.0 * e1 * e1 + 7.0 * t * t * e1 * e1 ) / 2.0;
    double b12 = 1.0 - 3.0 * t * t + 3.0 * t * t * e1 - 3.0 * t * t * e1 * e1;
		   b12 = N * c * c * b12 / 6.0;
    double b40 = N * t * ( 0.0 - e1 ) / 2.0;
    double b22 = N * t * c * c * ( -4.0 + 5.0 * e1 - 9.0 * t * t * e1 ) / 12.0;
    double b04 = N * t * pow( c, 4.0 ) * ( 1.0 - t * t + e1 ) / 24.0;

    double l01 = N * c;
    double l11 = N * t * c * ( -1.0 + e1 - e1 * e1 );
    double l21 = N * c * ( -2.0 + 2.0 * e1 - 9.0 * t * t * e1 - 2.0 * e1 * e1 
		                   +18.0 * t * t * e1 * e1 ) / 6.0;
    double l03 = N * pow( c, 3.0 ) * ( 0.0 - t * t ) / 6.0;
    double l31 = N * c * t * ( -7.0 * e1 + 3.0 * t * t * e1 ) / 6.0;
    double l13 = N * pow( c, 3.0 ) * t * ( -1.0 + t * t - t * t * e1 ) / 6.0;

    db = latito - latifrom;
    dl = longito - longifrom;

    a1 =  b10 * db                  + b20 * pow( db, 2.0 )
        + b02 *      pow( dl, 2.0)  + b30 * pow( db, 3.0 )
        + b40 *      pow( db, 4.0)  + b22 * pow( db, 2.0 ) * pow(dl, 2.0)
        + b12 * db * pow( dl, 2.0)
        + b04 *      pow( dl, 4.0);

    a2 =  l01 * dl                      + l11 * db * dl
        + l21 *      pow( db, 2.0 ) * dl
        + l03 *      pow( dl, 3.0 )
        + l31 *      pow( db, 3.0 ) * dl
        + l13 * db * pow( dl, 3.0 );

    alfa = atan2( a2, a1 );

	if( alfa < 0.0 ) alfa += g_dfTWOPI;

    azimuth = alfa;
    distance = sqrt( a1 * a1 + a2 * a2 );
}

/******************************************************************************

  Horistow differential formula to compute latitude/longitude corrections 
  due to azimuth correction 

******************************************************************************/

void cCoorTrans::horistow( double delta_azimuth,
						   double delta_lati, double delta_longi,
						   double lati_begin, double longi_begin,
						   double &lati_end,  double &longi_end,
						   double semi_major,
						   double flatening_denominator )
{
    double t = tan( lati_begin );
    double c = cos( lati_begin );
    double e1 = 1.0 / flatening_denominator; 
	e1 = 2.0 * e1 - e1 * e1;
    e1 = e1 / ( 1.0 - e1 ) * c * c;

    double db = lati_end- lati_begin;
    double dl = longi_end - longi_begin;

    double p1 = - c * c * ( 1.0 + t * t + e1 ) * dl * dl / 2.0 
		        - 1.5 * e1 * ( 1.0 - t * t ) * db * db
				- 3.0 * t * e1 * ( 1.0 - e1 ) * db + 1.0;
//  double p3=db-1.5*t*t*e1*db*db-0.5*c*c*t*(1.0+e1)*dl*dl-c*c*db*dl*dl/3.0;
    double p4 = - c * ( 1.0 + e1 ) * dl
		        + 3.0 * c * t * e1 * db * dl 
				+ c * c * c * ( 1.0 + t * t ) * dl * dl * dl / 6.0;

    double q1 =   t * ( 1.0 - e1 + e1 * e1 ) * dl
		        + ( 1.0 + t * t - e1 - 2.0 * t * t * e1 ) * db * dl
				+ t * ( 1.0 + t * t ) * db * db * dl
				- t * dl * dl * dl / 6.0;
//  double q3= dl+t*(1.0-e1)*db*dl+(2.0+3.0*t*t)*db*db*dl/3.0-c*c*t*t*dl*dl*dl/6.0;
    double q4 =   ( 1.0 - e1 + e1 * e1 ) * db / c 
		        + t * ( 1.0 - 0.5 * e1 ) * db * db / c
				- 0.5 * c * t * dl * dl
				+ ( 1.0 + 3.0 * t * t ) * pow( db, 3.0 ) / c / 3.0
				- 0.5 * c * ( 1.0 + t * t ) * db * dl * dl;

    double b2 = p1 * delta_lati + p4 * delta_azimuth;
    double l2 = q1 * delta_lati + q4 * delta_azimuth;

    lati_end  = lati_end + b2;
    longi_end = longi_end + l2 + delta_longi;
}

/******************************************************************************

  convert from DDD.MMSSSSSS to radians

******************************************************************************/

void cCoorTrans::DDMMSS2RAD( double &deg )
{
   int sign = 1;
   if( deg < 0) sign = -1;
   deg *= sign;

   int degree = ( int ) ( deg );
   double temp1 = deg - 1.0 * degree;
   int minute = ( int ) ( temp1 * 100.0 );
   temp1 = deg - 1.0 * degree - minute / 100.0;
   double seconds = temp1 * 10000.0;
   deg = degree + minute / 60.0 + seconds / 3600.0;
   deg = sign * deg * g_dfDEG2RAD;
}

/******************************************************************************

  Convert from radians to DD.MMSSSSSSS

******************************************************************************/

void cCoorTrans::RAD2DDMMSS( double &rad )
{
   int sign = 1;
   if( rad < 0 ) sign = -1;
   rad *= sign;

   rad = rad * g_dfRAD2DEG;
   int degree = ( int ) ( rad );
   double temp1 = ( rad - degree ) * 60.0;
   int minute = ( int ) ( temp1 );
   double seconds = ( temp1 - minute ) * 60.0;
   rad = degree + minute / 100.0 + seconds / 10000.0;
   rad = rad * sign;
}

void cCoorTrans::DegreeFraction( double deg, 
								 int &Deg, int &Min, int &Sec, int &SecDec )
{
	int sign = 1;
	if( deg < 0 ) sign = -1;
	double DEG = deg * sign;

	Deg = ( int ) DEG;
	double temp1 = DEG - 1.0 * Deg;
	Min = ( int ) ( temp1 * 100.0 );
	temp1 = DEG - 1.0 * Deg - Min / 100.0;
	double seconds = temp1 * 10000.0;
	Sec= ( int ) seconds;

	SecDec = ( int ) ( ( seconds - Sec ) * 10000.0 );

	Deg = sign * Deg;
}

/******************************************************************************

  Convert from rectangular topocentric delta x, y, z to azimuth, elevation and 
  topocentric distance 

******************************************************************************/

void cCoorTrans::DeltaXYZ_to_EleAzi( double dX, double dY, double dZ,
									 double Lat, double Longi,
									 double &elevation, double &azimuth, 
									 double &distance )
{ 
	double height=   cos( Longi ) * cos( Lat ) * dX
				   + sin( Longi ) * cos( Lat ) * dY
				   + sin( Lat   ) * dZ;
    double north = - cos( Longi ) * sin( Lat ) * dX
				   - sin( Longi ) * sin( Lat ) * dY
				   + cos( Lat   ) * dZ;
	double east =  - sin( Longi ) * dX
				   + cos( Longi ) * dY;

    distance = sqrt( dX * dX + dY * dY + dZ * dZ );
	elevation = asin( height / distance );

    azimuth = atan2( east, north );
    if( azimuth < 0.0 ) azimuth += g_dfTWOPI;
}

void cCoorTrans::EleAziDis_to_DeltaXYZ( double azimuth, double elevation, double distance,
									    double Lat, double Longi,
									    double &dX, double &dY, double &dZ )
{
	double North = cos( elevation ) * cos( azimuth ) * distance;
	double East = cos( elevation ) * sin( azimuth ) * distance;
	double Up = sin( elevation ) * distance;

	double dfSinLat = sin( Lat ), dfCosLat = cos( Lat );
	double dfSinLon = sin( Longi ), dfCosLon = cos( Longi );
	dX = -dfCosLon * dfSinLat * North - dfSinLon * East + dfCosLon * dfCosLat * Up;
	dY = -dfSinLon * dfSinLat * North + dfCosLon * East + dfSinLon * dfCosLat * Up;
	dZ =  dfCosLat * North + dfSinLat * Up;
}

/******************************************************************************

  Convert from rectangular topocentric delta x, y, z to azimuth, elevation and 
  topocentric distance 

******************************************************************************/

void cCoorTrans::DeltaXYZ_to_EleAzi( double dX, double dY, double dZ,
									 double dfSinLat, double dfCosLat,
									 double dfSinLong, double dfCosLong,
									 double &elevation, double &azimuth, 
									 double &distance )
{ 
	double height=  dfCosLong * dfCosLat * dX + dfSinLong * dfCosLat * dY + dfSinLat * dZ;
    double north = -dfCosLong * dfSinLat * dX - dfSinLong * dfSinLat * dY + dfCosLat * dZ;
	double east =  -dfSinLong * dX + dfCosLong * dY;

    distance = sqrt( dX * dX + dY * dY + dZ * dZ );
	elevation = asin( height / distance );

    azimuth = atan2( east, north );
    if( azimuth < 0.0 ) azimuth += g_dfTWOPI;
}

/******************************************************************************

  Convert rectangular x, y, z to RA and DEC

******************************************************************************/

void cCoorTrans::DeltaXYZ_to_RAD( double dX,  double dY, double dZ,
				        		  double &RA, double &D, double &distance )
{
    distance = sqrt( dX * dX + dY * dY + dZ * dZ );
	D = asin( dZ / distance );

    RA = atan2( dY, dX );
    if( RA < 0.0 ) RA += g_dfTWOPI;
}

/******************************************************************************

  Compute the rotation matrix

******************************************************************************/

void cCoorTrans::RotatingAroundXAxis( double RotationAngle, 
									  double &X, double &Y, double &Z )
{
	double x = X, y = Y, z = Z;

	X =  x;
	Y =  cos( RotationAngle ) * y + sin( RotationAngle ) * z;
	Z = -sin( RotationAngle ) * y + cos( RotationAngle ) * z;
}


void cCoorTrans::RotatingAroundYAxis( double RotationAngle, 
									  double &X, double &Y, double &Z )
{
	double x = X, y = Y, z = Z;

	X = cos( RotationAngle ) * x - sin( RotationAngle ) * z;
	Y = y;
	Z = sin( RotationAngle ) * x + cos( RotationAngle ) * z;
}

void cCoorTrans::RotatingAroundZAxis( double RotationAngle, 
									  double &X, double &Y, double &Z )
{
	double x = X, y = Y, z = Z;

	X =  cos( RotationAngle ) * x + sin( RotationAngle ) * y;
	Y = -sin( RotationAngle ) * x + cos( RotationAngle ) * y;
	Z =  z;
}

void cCoorTrans::ComputedXdYdZFromAzimuthElevation( double Azimuth, 
												    double Elevation, 
										            double &dx, 
										            double &dy, 
										            double &dz )
{
	dx = cos( Elevation ) * cos( Azimuth );
	dy = cos( Elevation ) * sin( Azimuth );
	dz = sin( Elevation );
}

void cCoorTrans::ComputeAzimuthElevationFromdXdYdZ( double &Azimuth, 
												    double &Elevation, 
										            double dx, 
										            double dy, 
										            double dz )
{
	Azimuth = atan2( dy, dx );
	if( Azimuth < 0 ) Azimuth += g_dfTWOPI;

	double r = sqrt( dx * dx + dy * dy );

	Elevation = atan2( dz, r );
}

/***************************************************************************************

  Fortran subroutine: compute_ECEF_unit_vector_of_east_north_height

  compute unit vectors of East (E), North (N) and Vertical (V) in ECEF system
  see Eqs (5.2-1, 5.2-2, 5.2-3), p. 5-10, MicroCosm Vol .1

  dflat and dfLong are the geodetic lat and lon of the point of interest

  pdfEast, pdfNorth and pdfUp are the unit vector of the East, North and Up at the 
  point

***************************************************************************************/

void cCoorTrans::ComputeECEFUnitVector( double dfLat, double dfLong, 
		                                double *pdfEast, double *pdfNorth, double *pdfUp )
{
	double cb = cos( dfLat );
	double sb = sin( dfLat );
	double cl = cos( dfLong );
	double sl = sin( dfLong );

	// unit vector of East in ECEF system
	pdfEast[ 0 ] = -sl;		
    pdfEast[ 1 ] =  cl;
    pdfEast[ 2 ] =  0.0;

	// unit vector of North in ECEF system 
    pdfNorth[ 0 ] = -sb * cl;
    pdfNorth[ 1 ] = -sb * sl;
    pdfNorth[ 2 ] =  cb;

	// unit vector of Vertical in ECEF system 
    pdfUp[ 0 ] = cb * cl;
    pdfUp[ 1 ] = cb * sl;
    pdfUp[ 2 ] = sb;
}

/***************************************************************************************

  Conversion between cartesian positiona nd geodetic position

  bXYZ2BLH is true is conversion from xyz 2 blh is asked, otherwise false

***************************************************************************************/

void cCoorTrans::ConvertXYZ2BLH( double &x, double &y, double &z ,
								 double &b, double &l, double &h,
								 bool bXYZ2BLH, 
								 double semi_major, double eccentricity )
{
	double an;

	if( !bXYZ2BLH )	// from BLH to XYZ
	{
		double sb = sin( b );
		double cb = cos( b );

		an = semi_major / sqrt( 1.0 - eccentricity * sb * sb );
		x = ( an + h ) * cb * cos( l );
		y = ( an + h ) * cb * sin( l );
		z = ( an * ( 1.0 - eccentricity ) + h ) * sb;
	}
	else
	{
		l = atan2( y, x );
		if( l < 0.0 ) l += g_dfTWOPI;

		if( z == 0.0 ) 
		{
			b = 0;
			an = semi_major;
		}
		else if( sqrt( x * x + y * y ) == 0.0 )
		{
			if( z < 0 ) b = -g_dfHALFPI;
			else b = g_dfHALFPI;

			an = semi_major / sqrt( 1.0 - eccentricity );
		}
		else
		{
			double zxy = z / sqrt( x * x + y * y );
			double b0 = atan( zxy );
			an = 0;
			 
			for( int i = 0; i < 10; i++ )
			{
				an = semi_major / sqrt( 1.0 - eccentricity * sin( b0 ) * sin( b0 ) );
				b0 = b = atan( zxy * ( 1.0 + eccentricity * an * sin( b0 ) / z ) );
			}
		}
    
		h = sqrt ( x * x + y * y ) / cos( b ) - an;
	}	
}
	
/***************************************************************************************

  Compute the rotation matrix pdfRotationMatrix

  nAxis:	1 for x axis, 2 for y axis, 3 for z axis
  dfAngle:	rotation angle in radian

***************************************************************************************/

void cCoorTrans::ComputeRotationMatrix( int nAxis, double dfAngle, 
									    double *pdfRotationMatrix )
{
	double ca, sa;

	ca = cos( dfAngle );
	sa = sin( dfAngle );

	g_Matrix.Set_Vec_Zero( pdfRotationMatrix, 9 );
	pdfRotationMatrix[ 0 ] = pdfRotationMatrix[ 4 ] = pdfRotationMatrix[ 8 ] = ca;

	if( nAxis == 1 )	// about x axis
	{
		pdfRotationMatrix[ 0 ] = 1.0;
		pdfRotationMatrix[ 5 ] = sa;
		pdfRotationMatrix[ 7 ] = -sa;
	}

	if( nAxis == 2 )	// about y axis
	{
		pdfRotationMatrix[ 4 ] = 1.0;
		pdfRotationMatrix[ 2 ] = -sa;
		pdfRotationMatrix[ 6 ] = sa;
	}

	if( nAxis == 3 )	// about z axis
	{
		pdfRotationMatrix[ 8 ] = 1.0;
		pdfRotationMatrix[ 1 ] = sa;
		pdfRotationMatrix[ 3 ] = -sa;
	}
}	

/******************************************************************************

  Conversion between the local North/East/Up and Delta X, Y and Z given the
  latitide and longitude of the location of interest

  bNEU2DXYZ: true from NEU to XYZ, else XYZ to NEU

  latitude and longitude in radians

******************************************************************************/

void cCoorTrans::ConvertNEU2DXYZ( double &dX, double &dY, double &dZ,
		                          double &North, double &East, double &Up,
						          double latitude, double longitude, 
								  bool bNEU2DXYZ )
{
	double dfCosLong = cos( longitude );
	double dfSinLong = sin( longitude );
	double dfCosLat = cos( latitude );
	double dfSinLat = sin( latitude );

	if( !bNEU2DXYZ )	// from XYZ to NEU
	{
		Up =  dfCosLong * dfCosLat * dX + dfSinLong * dfCosLat * dY + dfSinLat * dZ;
		North = -dfCosLong * dfSinLat * dX - dfSinLong * dfSinLat * dY + dfCosLat * dZ;
		East =  -dfSinLong * dX + dfCosLong * dY;
	}
	else	// from NEU to XYZ
	{
		dX = -dfCosLong * dfSinLat * North - dfSinLong * East + dfCosLong * dfCosLat * Up;
		dY = -dfSinLong * dfSinLat * North + dfCosLong * East + dfSinLong * dfCosLat * Up;
		dZ =  dfCosLat * North + dfSinLat * Up;
	}
}

/******************************************************************************

  Conversion between the local North/East/Up and Delta X, Y and Z given the
  latitide and longitude of the location of interest

  pdfNEU[ 3 ] ir oder of north, east and up (height)

  bNEU2DXYZ: true from NEU to XYZ, else XYZ to NEU

  latitude and longitude in radians

******************************************************************************/

void cCoorTrans::ConvertNEU2DXYZ( double *pdfDeltaXYZ, double *pdfNEU, 
						          double latitude, double longitude, 
								  bool bNEU2DXYZ )
{
	ConvertNEU2DXYZ( pdfDeltaXYZ[ 0 ], pdfDeltaXYZ[ 1 ], pdfDeltaXYZ[ 2 ],
		             pdfNEU[ 0 ], pdfNEU[ 1 ], pdfNEU[ 2 ],
				     latitude, longitude, bNEU2DXYZ );
}

/******************************************************************************

  convert 3D positional difference into along, across and radial bias

  bias[ 0 ] = along track
  bias[ 1 ] = across (normal to) track
  bias[ 2 ] = radial (pointing to earth centre)

******************************************************************************/

void cCoorTrans::ConvertPosDif2OrbitBias( double *pos, double *vel, double *posDif, double *bias )
{
	// radial
	double r = sqrt( g_Matrix.Vec_Multiply_Vec( pos, pos, 3 ) );
	double dr = g_Matrix.Vec_Multiply_Vec( posDif, pos, 3 );
	
	bias[ 2 ] = dr / r;	

	// normal to orbit plane
	double normal[ 3 ];

	g_Matrix.VectorProduct( pos, vel, normal );

	r = sqrt( g_Matrix.Vec_Multiply_Vec( normal, normal, 3 ) );
	dr = g_Matrix.Vec_Multiply_Vec( posDif, normal, 3 );

	bias[ 1 ] = dr / r;	

	// along track
	double along[ 3 ];

	g_Matrix.VectorProduct( normal, pos, along );

	r = sqrt( g_Matrix.Vec_Multiply_Vec( along, along, 3 ) );
	dr = g_Matrix.Vec_Multiply_Vec( posDif, along, 3 );
	
	bias[ 0 ] = dr / r;	// along
}

void cCoorTrans::ConvertPosDif2OrbitBias( double *pos, double *vel, double *posDif, double *bias, double &flight_angle )
{
	// radial
	double r = sqrt( g_Matrix.Vec_Multiply_Vec( pos, pos, 3 ) );
	double dr = g_Matrix.Vec_Multiply_Vec( posDif, pos, 3 );
	
	bias[ 2 ] = dr / r;	

	// normal to orbit plane
	double normal[ 3 ];

	g_Matrix.VectorProduct( pos, vel, normal );

	r = sqrt( g_Matrix.Vec_Multiply_Vec( normal, normal, 3 ) );
	dr = g_Matrix.Vec_Multiply_Vec( posDif, normal, 3 );

	bias[ 1 ] = dr / r;	

	// along track
	double along[ 3 ];

	g_Matrix.VectorProduct( normal, pos, along );

	r = sqrt( g_Matrix.Vec_Multiply_Vec( along, along, 3 ) );
	dr = g_Matrix.Vec_Multiply_Vec( posDif, along, 3 );
	
	bias[ 0 ] = dr / r;	// along

	// flight angle
	double v = sqrt( g_Matrix.Vec_Multiply_Vec( vel, vel, 3 ) );
	flight_angle = g_Matrix.Vec_Multiply_Vec( vel, along, 3 ) / r / v;
	flight_angle = acos( flight_angle );
}

/******************************************************************************

  convert along, across and radial bias into 3D positional difference

  bias[ 0 ] = along track
  bias[ 1 ] = across (normal to) track
  bias[ 2 ] = radial (pointing to earth centre)

******************************************************************************/

void cCoorTrans::ConvertOrbitBias2PosDif( double *pos, double *vel, double *bias, double *posDif )
{
	double trans[ 9 ];

	// radial
	double r = sqrt( g_Matrix.Vec_Multiply_Vec( pos, pos, 3 ) );
	
	trans[ 2 ] = pos[ 0 ] / r;
	trans[ 5 ] = pos[ 1 ] / r;
	trans[ 8 ] = pos[ 2 ] / r;

	// normal to orbit plane
	double normal[ 3 ];

	g_Matrix.VectorProduct( pos, vel, normal );

	r = sqrt( g_Matrix.Vec_Multiply_Vec( normal, normal, 3 ) );

	trans[ 1 ] = normal[ 0 ] / r;
	trans[ 4 ] = normal[ 1 ] / r;
	trans[ 7 ] = normal[ 2 ] / r;

	// along track
	double along[ 3 ];

	g_Matrix.VectorProduct( normal, pos, along );

	r = sqrt( g_Matrix.Vec_Multiply_Vec( along, along, 3 ) );

	trans[ 0 ] = along[ 0 ] / r;
	trans[ 3 ] = along[ 1 ] / r;
	trans[ 6 ] = along[ 2 ] / r;
	
	g_Matrix.Mat_Multiply_Vec( trans, bias, posDif, 3, 3 );
}


// North = S cos( El ) cos( Az )
// East  = S cos( El ) cos( Az )
// Up    = S sin( El )
void cCoorTrans::ConvertDeltaAzElDis2DeltaNEU( double az, double el, double dis, 
		                                       double delta_az, double delta_el, double delta_dis, 
									           double &delta_north, double &delta_east, double &delta_up )
{
	double cE = cos( el ), sE = sin( el );
	double cA = cos( az ), sA = sin( az );

	delta_north = delta_dis * cE * cA - dis * sE * cA * delta_el - dis * cE * sA * delta_az;
	delta_east  = delta_dis * cE * sA - dis * sE * sA * delta_el + dis * cE * cA * delta_az;
	delta_up    = delta_dis * sE + dis * cE * delta_el;
}

// all variables are in unit of radian
void cCoorTrans::convertAzElBias2AlongCrossBias( double azBias, double elBias, double El, double azDot, double elDot,
		                                         double &along, double &cross )
{
	double c = cos( El );
	double c2 = c * c;
	double velocity = sqrt( azDot * azDot * c2 + elDot * elDot );	

	if( fabs( velocity ) < 1.0e-5 ) 
	{
		cross = 0;
		along = 0;
		return;
	}

	cross = ( elDot * azBias - azDot * elBias ) * c / velocity;	
	along = ( azDot * c2 * azBias + elDot * elBias ) / velocity;
}

// all variables are in unit of radian
void cCoorTrans::convertAlongCrossBias2AzElBias( double along, double cross, double El, double azDot, double elDot,
		                                         double &deltaAz, double &deltaEl )
{
	double c = cos( El );
	double c2 = c * c;
	double velocity = sqrt( azDot * azDot * c2 + elDot * elDot );	

	if( fabs( velocity ) < 1.0e-5 ) 
	{
		deltaAz = 0;
		deltaEl = 0;
		return;
	}

	deltaAz = ( along * azDot + cross * elDot / c ) / velocity;
	deltaEl = ( along * elDot - cross * azDot * c ) / velocity;

}
