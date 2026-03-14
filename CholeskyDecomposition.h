/***************************************************************************************

 Cholesky Decomposition


 Reference: Section 2.9, pp. 89-91, Numerical Recipes in Fortran, William Press, et al.

***************************************************************************************/

#if!defined INC_CHOLESKYSECOM
#define INC_CHOLESKYSECOM

#include <windows.h>

class cCholeskyDecom
{
public:

	cCholeskyDecom();
	~cCholeskyDecom();

	BOOL Decompose( double *pdfMatrix, int nDim );
	BOOL Decompose( double *pdfMatrix, int nDim, int nDim0 );
	BOOL LinearEquation( double *pdfMatrix, double *pdfB, int nDim, double *pdfX );
	BOOL InverseLowTriangular( double *pdfMatrix, int nDim );
	BOOL MatrixInversing( double *pdfMatrix, int nDim, int nCase );
	BOOL MatrixInversing( double *pdfMatrix, int nDim, int nDim0, int nCase );
};

#endif