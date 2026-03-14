/***************************************************************************************

 polar motion parameters and rotation matrix

  Jizhang Sang, EOS

 (C) COPYRIGHT -Electro Optic Systems Pty. Ltd. Australia (EOS) 2002,
 All rights reserved. No part of this program may be photocopied, 
 reproduced or translated to another program language without the prior
 written consent of EOS

 Reference: IERS TECH NOTE 21

***************************************************************************************/

#if!defined INC_POLARMOTION
#define INC_POLARMOTION


#include <vector>

#include "EOP.h"

#include "LagrangianInterpolation.h"

struct stPolarMotion
{
	double dfJD;		
	double dfX, dfY;	// in radians

};

extern class cPolarMotion g_PolarMotion;

class cPolarMotion
{
	int m_nNumber;
	std::vector <stPolarMotion> m_pstPM;

	double m_dfJD, m_dfX, m_dfY;
	double m_pdfPMMatrixMean2True[ 9 ];
	double m_pdfPMMatrixTrue2Mean[ 9 ];

	cLagInterpolation m_LagInterpolation;

public:

	cPolarMotion();
	~cPolarMotion();

	bool SetPM();
	bool SetPMLOCAL(cEOP&);
	bool SetPM( double dfJD, double dfX, double dfY );

	bool GetPM( double dfJD, double &dfX, double &dfY, bool bMean2True );
	bool GetPMMatrixMean2True( double *pdfPMMatrix );
	bool GetPMMatrixMean2True( double dfJD, double *pdfPMMatrix );
	bool GetPMMatrixTrue2Mean( double *pdfPMMatrix );
	bool GetPMMatrixTrue2Mean( double dfJD, double *pdfPMMatrix );

private:

	bool ComputePMMatrixMean2True();
	bool ComputePMMatrixMean2True( double dfX, double dfY, double *pdfMatrix );
	bool ComputePMMatrixTrue2Mean();
	bool ComputePMMatrixTrue2Mean( double dfX, double dfY, double *pdfMatrix );
};

#endif