/************************************************************************

  Least squares fit of a set of observation/times in the form of a 
  polynomial:

    obs = c0 + t * c1 + t^2 * c2 + t^3 * c3 + ... + t^n * cn

	where n is the order of the polynomial

  In order for the fitting to be successful, the number of observations 
  must be equal to or more than (n+1).

  There is no outlier detection and removal measures used in the current 
  version. 

*************************************************************************/

#ifndef POLYFIT_H
#define POLYFIT_H

#include "Matrix.h"

class cPolyFit
{
	bool m_bInit;

	cMatrix m_Matrix;

	double *m_pdfObs;
	double *m_pdfTime;
	double *m_pdfResidual;
	int m_nNumberObs;
	int m_nNumberCoefficients;

	double *m_pdfCoe;

	double m_dfSigma;

	bool SolveCoefficients( bool *mark );
	bool SolveCoefficients(double SigmaThreshold, bool* mark);
	bool SolveCoefficients();
	bool SolveCoefficientsAzEl();
	void ComputeCoeff( double t, double *a );
	void ComputeCoeffAzEl( double t, double *a );


public:

	cPolyFit();
	~cPolyFit();
	
	bool SetObsAndTime( double *Obs, double *Time, int NumberObs, int PolyOrder );
	bool Fitting( bool *mark);
	bool Fitting(double SigmaThreshold, bool* mark);
	bool Fitting();
	bool FittingAzEl();

	double ComputeFittedValueAtGivenTime( double t );
	double ComputeFittedAzElValueAtGivenTime( double t );
	void GetCoefficients(int &numberCoefficients, double *pdfCoeffs );
	void SetCoefficients(int numberCoefficients, double *pdfCoeffs );
	double GetSigma();
	void GetResiduls( double *residual );

	double ComputeDotAtTimeT( double t );
};

#endif
