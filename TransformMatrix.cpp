/***************************************************************************************

 compute transformation matrices between different systems. The systems considered are:

  EFEC XYZ (Mean),
  EECE XYZ (true),
  TOD (True of Date): Instantaneous Inertial System at TOD,
  J2000,
  TOR (True of Reference): the instantaneous Inertial System at TOR, 
                           where the integration is carried out. For example,
						   TOR can be chosen as the epoch within the observation 
						   periods.


 Reference: IERS TECH NOTE 21

***************************************************************************************/

#include <math.h>
#include <stdlib.h>

#include "Constant.h"
#include "Matrix.h"
#include "TransformMatrix.h"
#include "EOP.h"
#include "DateTimeZ.h"

#include "PolarMotion.h"
#include "GreenwichSiderealTime.h"
#include "NPElement.h"

#if defined(_WIN32) || defined(_WIN64)
		// Windows header
		#include <windows.h>
		// Windows API
#else
		// Linux header
		#include <sys/stat.h>
		#include <unistd.h>
		// Linux API
#endif

cTransMatrix g_TransMatrix;

cTransMatrix::cTransMatrix()
{

}


cTransMatrix::~cTransMatrix()
{

}

/***************************************************************************************

  This method is called when the TOR (True of Reference epoch is determined at the
  begin of the application running.

  Compute the transformation matrices between J2000 <--> TOR, and TOR <--> ECEF mean.

  dfIntJD and dfFractionJD are the integer and fractional parts of TOR epoch in ET

***************************************************************************************/

bool cTransMatrix::ComputeTransMatrixJ20002TOR( double dfIntJD, double dfFractionJD )
{
	try
	{
		double dfJD = dfIntJD + dfFractionJD;

		if( !J20002TOD( dfJD, m_pdfMatrixJ20002TOR ) ) 
		{
			return false;
		}

		GetTranspose( m_pdfMatrixJ20002TOR, m_pdfMatrixTOR2J2000 );

		// TOR to ECEF Mean
		if( !TOR2ECEFMean( dfIntJD, dfFractionJD, m_pdfMatrixTOR2ECEFMean ) ) 
		{
			return false;
		}

		// ECEF mean to TOR
		GetTranspose( m_pdfMatrixTOR2ECEFMean, m_pdfMatrixECEFMean2TOR );

		m_dfTORJD = dfIntJD + dfFractionJD;
		m_dfTORIntJD = dfIntJD; 
		m_dfTORFractionJD = dfFractionJD;

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

 FORTRAN Subroutine: compute_transformation_matrix_between_ECEF_TOR_TOD

 Compute the transformation matrices between ECEF Mean, TOR, TOD at the given epoch

 dfIntJD and dfFractionJD ar the integer and fractional parts of the JD of interest 
     epoch, it should be in ET system.

***************************************************************************************/

bool cTransMatrix::ComputeTransMatrixECEFTORTOD( double dfIntJD, double dfFractionJD )
{
	try
	{
		double dfJD = dfIntJD + dfFractionJD;

		// J2000 to TOD
		if( !J20002TOD( dfJD, m_pdfMatrixJ20002TOD ) ) 
		{
			return false;
		}

		GetTranspose( m_pdfMatrixJ20002TOD, m_pdfMatrixTOD2J2000 );

		// TOR to TOD
		Multiply( m_pdfMatrixJ20002TOD, m_pdfMatrixTOR2J2000, m_pdfMatrixTOR2TOD );
		GetTranspose( m_pdfMatrixTOR2TOD, m_pdfMatrixTOD2TOR );

		/************************************************************************************/
		// TOD to ECEF Mean
		if( !TOD2ECEFMean( dfIntJD, dfFractionJD, m_pdfMatrixTOD2ECEFMean ) ) 
		{
			return false;
		}

		GetTranspose( m_pdfMatrixTOD2ECEFMean, m_pdfMatrixECEFMean2TOD );

		// TOR to ECEF mean
		Multiply( m_pdfMatrixTOD2ECEFMean, m_pdfMatrixTOR2TOD, m_pdfMatrixTOR2ECEFMean );
		GetTranspose( m_pdfMatrixTOR2ECEFMean, m_pdfMatrixECEFMean2TOR );

		ECEFMean2TOR( dfIntJD, dfFractionJD );

		/************************************************************************************/

		/************************************************************************************
		// TOR to ECEF Mean
		if( !TOR2ECEFMean( dfIntJD, dfFractionJD, m_pdfMatrixTOR2ECEFMean ) ) 
		{
			return false;
		}

		// ECEF mean to TOR
		GetTranspose( m_pdfMatrixTOR2ECEFMean, m_pdfMatrixECEFMean2TOR );

		// TOD to ECEF mean
		Multiply( m_pdfMatrixTOR2TOD, m_pdfMatrixECEFMean2TOR, m_pdfMatrixECEFMean2TOD );
		GetTranspose( m_pdfMatrixECEFMean2TOD, m_pdfMatrixTOD2ECEFMean );

		************************************************************************************/
		
		m_dfTODJD = dfIntJD + dfFractionJD;
		m_dfTODIntJD = dfIntJD; 
		m_dfTODFractionJD = dfFractionJD;

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************
 
   Get matrix from ECEF mean to ECEF true

***************************************************************************************/

bool cTransMatrix::ECEFMean2True( double dfJD, double *pdfMatrix )
{
	return g_PolarMotion.GetPMMatrixMean2True( dfJD, pdfMatrix );
}


/***************************************************************************************

   Get matrix from ECEF true to ECEF mean

***************************************************************************************/

bool cTransMatrix::ECEFTrue2Mean( double dfJD, double *pdfMatrix )
{
	return g_PolarMotion.GetPMMatrixTrue2Mean( dfJD, pdfMatrix );
}


/***************************************************************************************

   Get matrix from ECEF mean to TOD

***************************************************************************************/

bool cTransMatrix::ECEFMean2TOD( double dfIntJD, double dfFractionJD, double *pdfMatrix )
{
	try
	{
		if( fabs( dfIntJD - m_dfTODIntJD ) < 1.0e-7 && 
			fabs( dfFractionJD - m_dfTODFractionJD ) < 1.0e-09 )
		{
			g_Matrix.Vec_Assign( m_pdfMatrixECEFMean2TOD, pdfMatrix, 9 );
			return true;
		}

		double pdfM1[ 9 ], pdfM2[ 9 ], dfJD = dfIntJD + dfFractionJD;

		if( !ECEFMean2True( dfJD, pdfM1 ) ) 
		{
			return false;
		}

		double dfEqn;

		if( !g_NPElement.GetEquinoxEquation( dfJD, dfEqn ) ) 
		{
			return false;
		}

		double dfIntJDUT, dfFraJDUT;

		if( !g_EOP.ConvertJDETUT( dfIntJD, dfFractionJD, dfIntJDUT, dfFraJDUT, true ) )
		{
			return false;
		}

		if( !cGreenwichST::ComputeGSTMatrixECEF2TOD( dfIntJDUT, 
													 dfFraJDUT, 
													 dfEqn, 
													 pdfM2     ) ) 
		{
			return false;
		}

		Multiply( pdfM2, pdfM1, pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get matrix from TOD to ECEF mean

***************************************************************************************/

bool cTransMatrix::TOD2ECEFMean( double dfIntJD, double dfFractionJD, double *pdfMatrix )
{
	try
	{
		if( fabs( dfIntJD - m_dfTODIntJD ) < 1.0e-7 && 
			fabs( dfFractionJD - m_dfTODFractionJD ) < 1.0e-09 )
		{
			g_Matrix.Vec_Assign( m_pdfMatrixTOD2ECEFMean, pdfMatrix, 9 );
			return true;
		}

		if( !ECEFMean2TOD( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		Transpose( pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

   Get matrix from ECEF true to TOD

***************************************************************************************/

bool cTransMatrix::ECEFTrue2TOD( double dfIntJD, double dfFractionJD, double *pdfMatrix )
{
	try
	{
		double dfEqn, dfJD = dfIntJD + dfFractionJD;

		if( !g_NPElement.GetEquinoxEquation( dfJD, dfEqn ) ) 
		{
			return false;
		}

		double dfIntJDUT, dfFraJDUT;

		if( !g_EOP.ConvertJDETUT( dfIntJD, dfFractionJD, dfIntJDUT, dfFraJDUT, true ) )
		{
			return false;
		}

		if( !cGreenwichST::ComputeGSTMatrixECEF2TOD( dfIntJDUT, dfFraJDUT, dfEqn, pdfMatrix ) ) 
		{
			return false;
		}

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Get matrix from TOD to ECEF true

***************************************************************************************/

bool cTransMatrix::TOD2ECEFTrue( double dfIntJD, double dfFractionJD, double *pdfMatrix )
{
	try
	{
		if( !ECEFTrue2TOD( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		Transpose( pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Fortran subroutine: compute_transformation_matrix_J2000_to_TOD

  Get the matrix from J2000 to TOD

  The matrix from J2000 to TOD at dfJD is N * P,
  where P is the precession matrix from J2000 to TOD, 
  and N is the nutation matrix from mean to true

  pdfMatrix is the required transform matrix vectorized by line

***************************************************************************************/

bool cTransMatrix::J20002TOD( double dfJD, double *pdfMatrix )
{
	try
	{
		if( fabs( dfJD - m_dfTODJD ) < 1.0e-7 )
		{
			g_Matrix.Vec_Assign( m_pdfMatrixJ20002TOD, pdfMatrix, 9 );
			return true;
		}

		double pdfP[ 9 ];

		// get precession matrix from J2000 to TOD
		if( !g_NPElement.GetPrecessionMatrix( dfJD, pdfP, false ) ) 
		{
			return false;
		}

		double pdfN[ 9 ];

		// get nuttaion matrix from mean to true
		if( !g_NPElement.GetNutationMatrix( dfJD, pdfN, false ) ) 
		{
			return false;
		}

		Multiply( pdfN, pdfP, pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

  Get the matrix from TOD to J2000

***************************************************************************************/

bool cTransMatrix::TOD2J2000( double dfJD, double *pdfMatrix )
{
	try
	{
		if( fabs( dfJD - m_dfTODJD ) < 1.0e-7 )
		{
			g_Matrix.Vec_Assign( m_pdfMatrixTOD2J2000, pdfMatrix, 9 );
			return true;
		}
		
		if( !J20002TOD( dfJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		Transpose( pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

 Get the matrix from TOR to TOD

 The matrix from TOR to TOD is M2 * M1,
 where M1 is the matrix from TOR to J2000, 
 and M2 is the matrix from J2000 to TOD

 Here, the matrix from J2000 to TOR is given

***************************************************************************************/

bool cTransMatrix::TOR2TOD( double dfJD, double *pdfMatrix )
{
	try
	{
		if( fabs( dfJD - m_dfTODJD ) < 1.0e-7 )
		{
			g_Matrix.Vec_Assign( m_pdfMatrixTOR2TOD, pdfMatrix, 9 );
			return true;
		}

		double pdfM2[ 9 ];

		if( !J20002TOD( dfJD, pdfM2 ) ) 
		{
			return false;
		}

		Multiply( pdfM2, m_pdfMatrixTOR2J2000, pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

 Get the matrix from TOD to TOR

***************************************************************************************/

bool cTransMatrix::TOD2TOR( double dfJD, double *pdfMatrix )
{
	try
	{
		if( fabs( dfJD - m_dfTODJD ) < 1.0e-7 )
		{
			g_Matrix.Vec_Assign( m_pdfMatrixTOD2TOR, pdfMatrix, 9 );
			return true;
		}

		if( !TOR2TOD( dfJD, pdfMatrix ) ) 
		{
			return false;
		}

		Transpose( pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

 Get the matrix from ECEF true to TOR

 This matrix is M2 * M1,
 where M1 is the matrix from ECEFTrue to TOD,
 and M2 is the matrix from TOD to TOR

 Here, the matrix from J2000 to TOR is given

***************************************************************************************/

bool cTransMatrix::ECEFTrue2TOR( double dfIntJD, double dfFractionJD, double *pdfMatrix )
{
	try 
	{
		double pdfM1[ 9 ], pdfM2[ 9 ], dfJD = dfIntJD + dfFractionJD;

		if( !ECEFTrue2TOD( dfIntJD, dfFractionJD, pdfM1 ) ) 
		{
			return false;
		}
		
		if( !TOD2TOR( dfJD, pdfM2 ) ) 
		{
			return false;
		}

		Multiply( pdfM2, pdfM1, pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

 Get the matrix from TOR to ECEF true

***************************************************************************************/

bool cTransMatrix::TOR2ECEFTrue( double dfIntJD, double dfFractionJD, double *pdfMatrix )
{
	try
	{
		if( !ECEFTrue2TOR( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}

		Transpose( pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

 Get the matrix from ECEF mean to TOR

 This matrix is M2 * M1,
 where M1 is the matrix from ECEF mean 2 ECEF true,
 and M2 is the matrix from ECEF true to TOR

 Here, the matrix from J2000 to TOR is given

***************************************************************************************/

bool cTransMatrix::ECEFMean2TOR( double dfIntJD, double dfFractionJD, double *pdfMatrix )
{
	try
	{
		if( fabs( dfIntJD - m_dfTORIntJD ) < 1.0e-7 && 
			fabs( dfFractionJD - m_dfTORFractionJD ) < 1.0e-09 )
		{
			g_Matrix.Vec_Assign( m_pdfMatrixECEFMean2TOR, pdfMatrix, 9 );
			return true;
		}

		double pdfM1[ 9 ], pdfM2[ 9 ], dfJD = dfIntJD + dfFractionJD;

		if( !ECEFMean2True( dfJD, pdfM1 ) ) 
		{
			return false;
		}

		if( !ECEFTrue2TOR( dfIntJD, dfFractionJD, pdfM2 ) ) 
		{
			return false;
		}

		Multiply( pdfM2, pdfM1, pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

 Get the matrix from TOR to ECEF true

***************************************************************************************/

bool cTransMatrix::TOR2ECEFMean( double dfIntJD, double dfFractionJD, double *pdfMatrix )
{
	try
	{
		if( !ECEFMean2TOR( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}

		Transpose( pdfMatrix );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

/***************************************************************************************
 
   State transformation from ECEF mean 2 true

***************************************************************************************/

bool cTransMatrix::TransECEFMean2True( double dfJD, double *pdfState )
{
	try
	{
		double pdfMatrix [ 9 ];

		if( !ECEFMean2True( dfJD, pdfMatrix ) ) 
		{
			return false;
		}

		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

   State transformation from ECEF true 2 mean

***************************************************************************************/

bool cTransMatrix::TransECEFTrue2Mean( double dfJD, double *pdfState )
{
	try
	{
		double pdfMatrix [ 9 ];

		if( !ECEFTrue2Mean( dfJD, pdfMatrix ) ) 
		{
			return false;
		}

		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

   State transformation from ECEF mean to TOD

***************************************************************************************/

bool cTransMatrix::TransECEFMean2TOD( double dfIntJD, double dfFractionJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];

		if( !ECEFMean2TOD( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  State transformation from TOD to ECEF mean

***************************************************************************************/

bool cTransMatrix::TransTOD2ECEFMean( double dfIntJD, double dfFractionJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];
		
		if( !TOD2ECEFMean( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Transform the velocity vector from ECEFMean to TOD

  the pdfState is the position of the object in TOD, in meter
  pdfVelocity is the velocity in m/s

***************************************************************************************/

bool cTransMatrix::TransECEFMean2TODVelocity( double dfIntJD, double dfFractionJD, 
											  double *pdfState, double *pdfVelocity )
{
	try
	{		
		if( !TransECEFMean2True( dfIntJD + dfFractionJD, pdfVelocity ) ) 
		{
			return false;
		}

		// from ECEF True to TOD
		if( !TransECEFTrue2TODVelocity( dfIntJD, dfFractionJD, pdfState, pdfVelocity ) ) 
		{
			return false;
		}
		
		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Transform the velocity vector from ECEFMean to TOD

  the pdfState is the position of the object in the ECEF, in meter
  pdfVelocity is the velocity in m/s

***************************************************************************************/
bool cTransMatrix::TransTOD2ECEFMeanVelocity( double dfIntJD, double dfFractionJD, 
											  double *pdfState, double *pdfVelocity )
{
	try
	{
		// from TOD to ECEF True
		// from mean to true ECEF
		double pdfStateTrue[ 3 ];
		g_Matrix.Vec_Assign( pdfState, pdfStateTrue, 3 );
		
		if( !TransTOD2ECEFTrueVelocity( dfIntJD, dfFractionJD, pdfStateTrue, pdfVelocity ) ) 
		{
			return false;
		}

		// from true to mean ECEF
		if( !TransECEFTrue2Mean( dfIntJD + dfFractionJD, pdfVelocity ) ) 
		{
			return false;
		}

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

   State transformation from ECEF true to TOD

***************************************************************************************/

bool cTransMatrix::TransECEFTrue2TOD( double dfIntJD, double dfFractionJD, 
									  double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];
		
		if( !ECEFTrue2TOD( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  State transformation from TOD to ECEF true

***************************************************************************************/

bool cTransMatrix::TransTOD2ECEFTrue( double dfIntJD, double dfFractionJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];
		
		if( !TOD2ECEFTrue( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Transform the velocity vector from ECEFTrue to TOD

  the pdfState is the position of the object in TOD, in meter
  pdfVelocity is the velocity in m/s

***************************************************************************************/

bool cTransMatrix::TransECEFTrue2TODVelocity( double dfIntJD, double dfFractionJD, 
											  double *pdfState, double *pdfVelocity )
{
	try
	{
		double dfGST, dfSinGST, dfCosGST, dfIntJDUT, dfFraJDUT;
		double dfJD = dfIntJD + dfFractionJD, dfEqn;

		if( !g_EOP.ConvertJDETUT( dfIntJD, dfFractionJD, dfIntJDUT, dfFraJDUT, true ) )
		{
			return false;
		}

		if( !g_NPElement.GetEquinoxEquation( dfJD, dfEqn ) )
		{
			return false;
		}

		cGreenwichST::ComputeGST( dfIntJDUT, dfFraJDUT, dfEqn, dfGST, dfSinGST, dfCosGST );

		double pdf[ 3 ];

		pdf[ 2 ] =  pdfVelocity[ 2 ];	// z axis
		pdf[ 0 ] =  pdfVelocity[ 0 ] * dfCosGST - 
					pdfVelocity[ 1 ] * dfSinGST - 
					pdfState[ 1 ] * g_dfEarthAngVelocity;
		pdf[ 1 ] =  pdfVelocity[ 0 ] * dfSinGST + 
					pdfVelocity[ 1 ] * dfCosGST + 
					pdfState[ 0 ] * g_dfEarthAngVelocity;

		g_Matrix.Vec_Assign( pdf, pdfVelocity, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Transform the velocity vector from ECEFTrue to TOD

  the pdfState is the position of the object in ECEF, in meter
  pdfVelocity is the velocity in m/s

***************************************************************************************/

bool cTransMatrix::TransTOD2ECEFTrueVelocity( double dfIntJD, double dfFractionJD, 
											  double *pdfState, double *pdfVelocity  )
{
	try
	{
		double dfGST, dfSinGST, dfCosGST, dfIntJDUT, dfFraJDUT;
		double dfJD = dfIntJD + dfFractionJD, dfEqn;

		if( !g_EOP.ConvertJDETUT( dfIntJD, dfFractionJD, dfIntJDUT, dfFraJDUT, true ) )
		{
			return false;
		}

		if( !g_NPElement.GetEquinoxEquation( dfJD, dfEqn ) )
		{
			return false;
		}

		cGreenwichST::ComputeGST( dfIntJDUT, dfFraJDUT, dfEqn, dfGST, dfSinGST, dfCosGST );

		double pdf[ 3 ];

		pdf[ 2 ] =  pdfVelocity[ 2 ];	// z axis
		pdf[ 0 ] =  pdfVelocity[ 0 ] * dfCosGST + 
					pdfVelocity[ 1 ] * dfSinGST + 
					pdfState[ 1 ] * g_dfEarthAngVelocity;
		pdf[ 1 ] = -pdfVelocity[ 0 ] * dfSinGST + 
					pdfVelocity[ 1 ] * dfCosGST - 
					pdfState[ 0 ] * g_dfEarthAngVelocity;

		g_Matrix.Vec_Assign( pdf, pdfVelocity, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Transform the velocity vector from ECEFTrue to TOR

  the pdfState is the position of the object in TOR, in meter
  pdfVelocity is the velocity in m/s

***************************************************************************************/

bool cTransMatrix::TransECEFTrue2TORVelocity( double dfIntJD, double dfFractionJD, 
											  double *pdfState, double *pdfVelocity )
{
	try
	{
		// ECEF true to TOD
		if( !TransECEFTrue2TODVelocity( dfIntJD, dfFractionJD, pdfState, pdfVelocity ) )
		{
			return false;
		}

		// then TOD to TOR
		if( !TransTOD2TOR( dfIntJD+dfFractionJD, pdfVelocity ) )
		{
			return false;
		}

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  Transform the velocity vector from ECEFTrue to TOD

  the pdfState is the position of the object in ECEF, in meter
  pdfVelocity is the velocity in m/s

***************************************************************************************/

bool cTransMatrix::TransTOR2ECEFTrueVelocity( double dfIntJD, double dfFractionJD, 
											  double *pdfState, double *pdfVelocity  )
{
	try
	{
		// TOR to TOD first
		if( !TransTOR2TOD( dfIntJD+dfFractionJD, pdfVelocity ) )
		{
			return false;
		}

		// TOD to ECEF true
		if( !TransTOD2ECEFTrueVelocity( dfIntJD, dfFractionJD, pdfState, pdfVelocity ) )
		{
			return false;
		}

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  State transformation from J2000 to TOR

***************************************************************************************/

bool cTransMatrix::TransJ20002TOR( double *pdfState )
{
	try
	{	
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( m_pdfMatrixJ20002TOR, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  State transformation from J2000 to TOD

***************************************************************************************/

bool cTransMatrix::TransJ20002TOD( double dfJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];
		
		if( !J20002TOD( dfJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

  State transformation from TOD to J2000

***************************************************************************************/

bool cTransMatrix::TransTOD2J2000( double dfJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];
		
		if( !TOD2J2000( dfJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

  State transformation from TOR to TOD

***************************************************************************************/

bool cTransMatrix::TransTOR2TOD( double dfJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];

		if( !TOR2TOD( dfJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

  State transformation from TOD to TOR

***************************************************************************************/

bool cTransMatrix::TransTOD2TOR( double dfJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];
		
		if( !TOD2TOR( dfJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  State transformation from ECEF True to TOR

***************************************************************************************/

bool cTransMatrix::TransECEFTrue2TOR( double dfIntJD, double dfFractionJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];

		if( !ECEFTrue2TOR( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

 State transformation from TOR to ECEF true

***************************************************************************************/

bool cTransMatrix::TransTOR2ECEFTrue( double dfIntJD, double dfFractionJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];
		
		if( !TOR2ECEFTrue( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

  State transformation from ECEF mean to TOR

***************************************************************************************/

bool cTransMatrix::TransECEFMean2TOR( double dfIntJD, double dfFractionJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];

		if( !ECEFMean2TOR( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

 State transformation from TOR to ECEF mean

***************************************************************************************/

bool cTransMatrix::TransTOR2ECEFMean( double dfIntJD, double dfFractionJD, double *pdfState )
{
	try
	{
		double pdfMatrix[ 9 ];
		
		if( !TOR2ECEFMean( dfIntJD, dfFractionJD, pdfMatrix ) ) 
		{
			return false;
		}
		
		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );

		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}


/***************************************************************************************

 State transformation from ECEF mean to TOR

***************************************************************************************/

//bool cTransMatrix::TransECEFMean2TOR( double *pdfState )
//{
//	double pdf[ 3 ];
//
//	g_Matrix.Mat_Multiply_Vec( m_pdfMatrixECEFMean2TOR, pdfState, pdf, 3, 3 );
//
//	g_Matrix.Vec_Assign( pdf, pdfState, 3 );
//
//	return true;
//}
//
//
///***************************************************************************************
//
// Get the matrix from TOR to ECEF mean
//
//***************************************************************************************/
//
//bool cTransMatrix::TransTOR2ECEFMean( double *pdfState )
//{
//	double pdf[ 3 ];
//
//	g_Matrix.Mat_Multiply_Vec( m_pdfMatrixTOR2ECEFMean, pdfState, pdf, 3, 3 );
//
//	g_Matrix.Vec_Assign( pdf, pdfState, 3 );
//
//	return true;
//}
//
/***************************************************************************************

 Convert the given state vector from the mean at dfJD1 to the mean at dfJD2 considering 
 only the precession

***************************************************************************************/

bool cTransMatrix::PrecessionTrans( double dfJD1, double dfJD2, double *pdfState )
{
	try
	{
		cPrecession Precession;
		
		double pdfPrecessionMatrix[ 9 ];
		
		if( !Precession.Precession( dfJD1, dfJD2, pdfPrecessionMatrix ) ) 
		{
			return false;
		}

		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfPrecessionMatrix, pdfState, pdf, 3, 3 );
		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}

/***************************************************************************************

 Convert the given state vector from the mean at the given dfJD to the true at the given
 dfJD or vice versa, considering only the nutation

  bMean2True: true from mean to true, false from true to mean

***************************************************************************************/

bool cTransMatrix::NutationTrans( double dfJD, double *pdfState, bool bMean2True )
{
	try
	{
		cNutation Nutation;

		double pdfMatrix[ 9 ];

		if( !Nutation.Nutation( dfJD, !bMean2True, pdfMatrix ) )
		{
			return false;
		}

		double pdf[ 3 ];

		g_Matrix.Mat_Multiply_Vec( pdfMatrix, pdfState, pdf, 3, 3 );
		g_Matrix.Vec_Assign( pdf, pdfState, 3 );

		return true;
	}
	catch( ... )
	{
		return false;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////


/***************************************************************************************

 Multilication of two matrix of 3*3

 All the matrices are vectorized by row

***************************************************************************************/

void cTransMatrix::Multiply( double *pdfMatrix1, double *pdfMatrix2, double *pdfMatrix3 )
{
	try
	{
		double tem1; 
		int k1, k2, k3, l1, l2;

		for ( k1 = 1 ; k1 <= 3 ; k1++ )
		{
			for ( k2 = 1 ; k2 <= 3 ; k2++ )
			{
				tem1 = 0;
			
				for ( k3 = 1 ; k3 <= 3 ; k3++ )
				{
					l1 = (k1-1) * 3 + k3 -1;
					l2 = (k3-1) * 3 + k2 -1;
					tem1 = tem1 + pdfMatrix1[ l1 ] * pdfMatrix2[ l2 ];
				}
			
				l1 = (k1-1) * 3 + k2 -1;
				pdfMatrix3[ l1 ] = tem1;
			}
		}	

		return;
	}
	catch( ... )
	{
	}
}


/***************************************************************************************

 transpose of matrix of 3*3

***************************************************************************************/

void cTransMatrix::Transpose( double *pdfMatrix )
{
	swap( pdfMatrix[ 1 ], pdfMatrix[ 3 ] );
	swap( pdfMatrix[ 2 ], pdfMatrix[ 6 ] );
	swap( pdfMatrix[ 5 ], pdfMatrix[ 7 ] );
}


/***************************************************************************************

 transpose of matrix of 3*3

***************************************************************************************/

void cTransMatrix::GetTranspose( double *pdfSource, double *pdfResult )
{
	for( int i = 0; i < 9; i++ ) pdfResult[ i ] = pdfSource[ i ];
	
	Transpose( pdfResult );
}

void cTransMatrix::swap( double &df1, double &df2 )
{
	double t = df1;
	df1 = df2;
	df2 = t;
}


void cTransMatrix::ECEFMean2TOR( double dfIntJD, double dfFraJD )
{
	double dfJD = dfIntJD + dfFraJD;
	
	// mean to true pole
	double f1[ 9 ];
	ECEFMean2True( dfJD, f1 );

	// true pole to TOD
	double f2[ 9 ];
	ECEFTrue2TOD( dfIntJD, dfFraJD, f2 );

	// TOD to J2000
	double f3[ 9 ];
	TOD2J2000( dfJD, f3 );

	double f4[ 9 ];

	Multiply( f2, f1, f4 );
	Multiply( f3, f4, f1 );
	Multiply( m_pdfMatrixJ20002TOR, f1, f4 );
}

bool cTransMatrix::TransTeme2Tod(double jd_utc, double *pdfState)
{
    // TOD pv to TEME pv
	double tt{0};
	g_EOP.ConvertETUTC(tt, jd_utc, false);
	double equinox{0};
	g_NPElement.GetEquinoxEquation(tt, equinox);

	// tod to teme matrix
	double teme_to_tod[9] = {cos(equinox), -sin(equinox), 0.0,
							 sin(equinox), cos(equinox),  0.0,
							 0.0,          0.0,           1.0};

	double teme_state[3];
	g_Matrix.Mat_Multiply_Vec(teme_to_tod, pdfState, teme_state, 3, 3);
	g_Matrix.Vec_Assign(teme_state, pdfState, 3);

	return true;
}