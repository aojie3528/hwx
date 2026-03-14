/***************************************************************************************

 compute Greenwich sidereal time


 Reference: IERS TECH NOTE 21

***************************************************************************************/

#if!defined INC_GREENWICHST
#define INC_GREENWICHST

#include <windows.h>

class cGreenwichST
{
public:

	cGreenwichST();
	~cGreenwichST();

	static double ComputeGST( double dfTJD );

	static BOOL ComputeGST( double dfIntJD, double dfFractionJD, 
		                    double dfEquinoxEquation,  double &dfGST );
	static BOOL ComputeGST( double dfIntJD, double dfFractionJD, 
		                    double dfEquinoxEquation, double &dfGST,  
							double &dfSinGST, double &dfCosGST );
	static BOOL ComputeGSTMatrixECEF2TOD( double dfIntJD, double dfFractionJD, 
		                                  double dfEquinoxEquation, double *pdfGSTMatrix );
	static BOOL ComputeGSTMatrixTOD2ECEF( double dfIntJD, double dfFractionJD, 
		                                  double dfEquinoxEquation, double *pdfGSTMatrix );
};

#endif

