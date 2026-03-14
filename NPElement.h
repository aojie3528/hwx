/***************************************************************************************

 about nutation and precession parameters


***************************************************************************************/

#if!defined INC_NPELEMENT
#define INC_NPELEMENT


#include <vector>
#include "Precession.h"
#include "Nutation.h"
#include "LagrangianInterpolation.h"

struct stNPElement
{
	double dfJD, dfMJD;

	double dfZeta, dfZee, dfTheta;	// in radian, precession elements
	double dfPSI, dfEPS;			// in radian, nutation elements
	double dfOblm;					// in radian mean obliquity

	double dfEqn;					// in radian, Equinox Equation
};

extern class cNPElement g_NPElement;

class cNPElement
{
	cPrecession m_Precession;
	cNutation m_Nutation;
	cLagInterpolation m_LagInterpolation;

	int m_nNumberEpochs;
	double m_dfStartJD, m_dfEndJD;
	std::vector <stNPElement> m_pstNPElement;

public:

	cNPElement();
	~cNPElement();

    int GetNumberEpochs() { return m_nNumberEpochs; }
	void SetNumberEpochs( int n ) { m_nNumberEpochs = n; }

	bool ComputeNPElement( double dfJDStart, double dfJDEnd );

	bool ComputeNPElement( double dfJD, stNPElement &tNPElement );
	bool ComputeNPMatrix( double dfJD, double *matrix, bool bJ2K2TOD );
	bool ComputeAndSetNPElement( double dfJD );

	bool GetNPElement( double dfJD, stNPElement &tNPElement );

	bool GetPrecessionElement( double dfJD, double &dfZeta, double &dfZee, double &dfTheta );
	bool GetPrecessionMatrix( double dfJD, double *pdfMatrix, bool bTOR2J2000 );

	bool GetNutationElement( double dfJD, double &dfPSI, double &dfEPS, double &dfOblm );
	bool GetNutationMatrix( double dfJD, double *pdfMatrix, bool bTrue2Mean );

	bool GetMatrixTOD2J2000( double dfJD, double *pdfMatrix );
	bool GetMatrixJ20002TOD( double dfJD, double *pdfMatrix );

	bool GetEquinoxEquation( double dfJD, double &dfEqn );
};

#endif
