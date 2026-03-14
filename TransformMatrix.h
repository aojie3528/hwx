/***************************************************************************************

 compute transformation matrices between different systems. The systems considered are:

  EFEC XYZ (Mean),
  EECE XYZ (TRUE),
  TOD (True of Date): Instantaneous Inertial System
  J2000
  TOR (True of Reference): the Inertial System where the integration is carried out.


 Reference: IERS TECH NOTE 21

***************************************************************************************/

#if!defined INC_TRANSFORMATIONMATRICES
#define INC_TRANSFORMATIONMATRICES

extern class cTransMatrix g_TransMatrix;

class cTransMatrix
{
	// the JD of the TOR epoch for m_pdfMatrixTOR2ECEFMean and m_pdfMatrixECEFMean2TOR
	double m_dfTORJD, m_dfTORIntJD, m_dfTORFractionJD;

	// the JD of the matrices ,
	//                        m_pdfMatrixTOD2ECEF  and m_pdfMatrixECEF2TOD,
	//                        m_pdfMatrixTOR2TOD   and m_pdfMatrixTOD2TOR,
	//                        m_pdfMatrixJ20002TOD and m_pdfMatrixTOD2J2000
	double m_dfTODJD, m_dfTODIntJD, m_dfTODFractionJD;

public:

//  TOR to ECEF mean
	double m_pdfMatrixTOR2ECEFMean[ 9 ];	// vectorized by line
	double m_pdfMatrixECEFMean2TOR[ 9 ];	// vectorized by line

//  TOD to ECEF mean
	double m_pdfMatrixTOD2ECEFMean[ 9 ];	// vectorized by line
	double m_pdfMatrixECEFMean2TOD[ 9 ];	// vectorized by line

//  TOR to TOD and TOD to TOR
	double m_pdfMatrixTOR2TOD[ 9 ];		// vectorized by line
	double m_pdfMatrixTOD2TOR[ 9 ];		// vectorized by line

//  J2000 to TOR and TOR to J2000
	double m_pdfMatrixJ20002TOR[ 9 ];	// vectorized by line
	double m_pdfMatrixTOR2J2000[ 9 ];	// vectorized by line

//  J2000 to TOD
	double m_pdfMatrixJ20002TOD[ 9 ];	// vectorized by line
	double m_pdfMatrixTOD2J2000[ 9 ];	// vectorized by line

public:

	cTransMatrix();
	~cTransMatrix();

	bool ComputeTransMatrixJ20002TOR( double dfIntJD, double dfFractionJD );

	bool ComputeTransMatrixECEFTORTOD( double dfIntJD, double dfFractionJD );

	bool ECEFMean2True( double dfJD, double *pdfMatrix );
	bool ECEFTrue2Mean( double dfJD, double *pdfMatrix );

	bool ECEFMean2TOD( double dfIntJD, double dfFractionJD, double *pdfMatrix );
	bool TOD2ECEFMean( double dfIntJD, double dfFractionJD, double *pdfMatrix );

	bool ECEFTrue2TOD( double dfIntJD, double dfFractionJD, double *pdfMatrix );
	bool TOD2ECEFTrue( double dfIntJD, double dfFractionJD, double *pdfMatrix );

	bool J20002TOD( double dfJD, double *pdfMatrix ); 
	bool TOD2J2000( double dfJD, double *pdfMatrix ); 

	bool TOR2TOD( double dfJD, double *pdfMatrix );
	bool TOD2TOR( double dfJD, double *pdfMatrix );

	bool ECEFTrue2TOR( double dfIntJD, double dfFractionJD, double *pdfMatrix );
	bool TOR2ECEFTrue( double dfIntJD, double dfFractionJD, double *pdfMatrix );

	bool ECEFMean2TOR( double dfIntJD, double dfFractionJD, double *pdfMatrix );
	bool TOR2ECEFMean( double dfIntJD, double dfFractionJD, double *pdfMatrix );

	// state transformation 
	bool TransECEFMean2True( double dfJD, double *pdfState );
	bool TransECEFTrue2Mean( double dfJD, double *pdfState );

	bool TransECEFMean2TOD( double dfIntJD, double dfFractionJD, double *pdfState );
	bool TransTOD2ECEFMean( double dfIntJD, double dfFractionJD, double *pdfState );

	bool TransECEFMean2TODVelocity( double dfIntJD, double dfFractionJD, 
		                            double *pdfState, double *pdfVelocity );
	bool TransTOD2ECEFMeanVelocity( double dfIntJD, double dfFractionJD, 
		                            double *pdfState, double *pdfVelocity );

	bool TransECEFTrue2TOD( double dfIntJD, double dfFractionJD, double *pdfState );
	bool TransTOD2ECEFTrue( double dfIntJD, double dfFractionJD, double *pdfState );

	bool TransECEFTrue2TODVelocity( double dfIntJD, double dfFractionJD, 
		                            double *pdfState, double *pdfVelocity );
	bool TransTOD2ECEFTrueVelocity( double dfIntJD, double dfFractionJD, 
		                            double *pdfState, double *pdfVelocity );

	bool TransECEFTrue2TORVelocity( double dfIntJD, double dfFractionJD, 
		                            double *pdfState, double *pdfVelocity );
	bool TransTOR2ECEFTrueVelocity( double dfIntJD, double dfFractionJD, 
		                            double *pdfState, double *pdfVelocity );

	bool TransJ20002TOR( double *pdfState );
	bool TransJ20002TOD( double dfJD, double *pdfState ); 
	bool TransTOD2J2000( double dfJD, double *pdfState ); 

	bool TransTOR2TOD( double dfJD, double *pdfState );
	bool TransTOD2TOR( double dfJD, double *pdfState );

	bool TransECEFTrue2TOR( double dfIntJD, double dfFractionJD, double *pdfState );
	bool TransTOR2ECEFTrue( double dfIntJD, double dfFractionJD, double *pdfState );

	bool TransECEFMean2TOR( double dfIntJD, double dfFractionJD, double *pdfState );
	bool TransTOR2ECEFMean( double dfIntJD, double dfFractionJD, double *pdfState );

//	bool TransECEFMean2TOR( double *pdfState );
//	bool TransTOR2ECEFMean( double *pdfState );

	bool PrecessionTrans( double dfJD1, double dfJD2, double *pdfState );
	bool NutationTrans( double dfJD, double *pdfState, bool bMean2True );

	bool TransTeme2Tod(double jd_utc, double *pdfState);

private:

	void Multiply( double *pdfMatrix1, double *pdfMatrix2, double *pdfMatrix3 );
	void Transpose( double *pdfMatrix );
	void GetTranspose( double *pdfSource, double *pdfResult );
	void swap( double &df1, double &df2 );

	void ECEFMean2TOR( double dfIntJD, double dfFraJD );
};


#endif

