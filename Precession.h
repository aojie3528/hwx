/***************************************************************************************

 compute precession parameters


***************************************************************************************/

#if!defined INC_PRECESSION
#define INC_PRECESSION



class cPrecession
{
	// three precession parameters, in radians
	double m_dfZeta, m_dfZee, m_dfTheta;

	double m_pdfPrecessionMatrixTOR2J2000[ 9 ];	// vectorized by line
	double m_pdfPrecessionMatrixJ20002TOR[ 9 ];	// vectorized by line

public:

	cPrecession();
	~cPrecession();

	bool ComputeParameters( double dfJD, double &dfZeta, double &dfZee, double &dfTheta );

	bool Precession( double dfJD, bool bTOR2J2000 );
	bool Precession( double dfJD, bool bTOR2J2000, double *dfPrecessionMatrix );
	bool Precession( double dfJD1, double dfJD2, double *pdfPrecessionMatrix );

	bool GetParameters( double &dfZeta, double &dfZee, double &dfTheta );
	bool GetPrecessionMatrixTOR2J2000( double dfZeta, double dfZee, double dfTheta,
		                               double *pdfPrecessionMatrix );
	bool GetPrecessionMatrixJ20002TOR( double dfZeta, double dfZee, double dfTheta,
		                               double *pdfPrecessionMatrix );
	bool GetPrecessionMatrixTOR2J2000( double *pdfPrecessionMatrix );
	bool GetPrecessionMatrixJ20002TOR( double *pdfPrecessionMatrix );

private:

	void ComputePrecessionMatrixTOR2J2000();
	void ComputePrecessionMatrixJ20002TOR();
};

#endif