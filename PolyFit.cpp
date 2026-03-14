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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include "Constant.h"
#include "PolyFit.h"

const int MaxNumberObs = 50000;
const int MaxNumberUnknowns = 100;

cPolyFit::cPolyFit()
{	
	m_pdfObs = new double[ MaxNumberObs ];
	m_pdfTime = new double[ MaxNumberObs ];
	m_pdfResidual = new double[ MaxNumberObs ];
	m_pdfCoe = new double[ MaxNumberUnknowns ];

	if( !m_pdfObs || !m_pdfTime || !m_pdfResidual || !m_pdfCoe ) m_bInit = false;
	else m_bInit = true;
}

cPolyFit::~cPolyFit()
{
	if( m_pdfObs ) delete []m_pdfObs;
	if( m_pdfTime ) delete []m_pdfTime;
	if( m_pdfCoe ) delete []m_pdfCoe;
	if( m_pdfResidual ) delete []m_pdfResidual;
}

/******************************************************************************

  compute the coefficients of the polynomial, ie 1.0, t, t^2, t^3, ..., t^n

******************************************************************************/

void cPolyFit::ComputeCoeff( double t, double *a )
{
    a[ 0 ] = 1.0; 
	for( int i = 1; i < m_nNumberCoefficients; i++ ) a[ i ] = a[ i - 1 ] * t;
}

void cPolyFit::ComputeCoeffAzEl( double t, double *a )
{
    a[ 0 ] = 1.0; 
	a[ 1 ] = t;
	a[ 2 ] = t * t;
	a[ 3 ] = pow( t, 1.2 );
	a[ 4 ] = pow( t, 1.4 );
	a[ 5 ] = pow( t, 1.6 );
	a[ 6 ] = pow( t, 1.8 );
}


/******************************************************************************

  Set the observations, time and order of polynomial.

******************************************************************************/

bool cPolyFit::SetObsAndTime( double *Obs, double *Time, int NumberObs, int PolyOrder )
{
	if( !m_bInit ) 
	{
		return false;
	}

	if( NumberObs > MaxNumberObs ) 	// two many obs 
	{
		return false;
	}

	m_nNumberObs = NumberObs;

    m_Matrix.Vec_Assign( Obs, m_pdfObs, m_nNumberObs );
    m_Matrix.Vec_Assign( Time, m_pdfTime, m_nNumberObs );

	m_nNumberCoefficients = PolyOrder + 1;	

	// to chech the number of observations are more than the number of coefficients
	if( m_nNumberCoefficients > m_nNumberObs ) 
	{
		return false;
	}

	return true;
}

/******************************************************************************

  after calling SetObsAndTime() method, User calls this function to perform
  the least squares polynomial fitting of the observations

  ******************************************************************************/

bool cPolyFit::Fitting( bool *mark )
{
    return SolveCoefficients( mark );
}

bool cPolyFit::Fitting(double SigmaThreshold, bool* mark)
{
	return SolveCoefficients(SigmaThreshold, mark);
}

bool cPolyFit::Fitting()
{
	return SolveCoefficients();
}

bool cPolyFit::FittingAzEl()
{
    return SolveCoefficientsAzEl();
}


/******************************************************************************

  Least squares polynomial fitting of a set of observations and times

******************************************************************************/

bool cPolyFit::SolveCoefficients()
{
	if (!m_bInit)
	{
		return false;
	}

	double *ata = new double[MaxNumberUnknowns * MaxNumberUnknowns];
	double *a = new double[MaxNumberUnknowns];
	double *atl = new double[MaxNumberUnknowns];
	double *coe_approximate = new double[MaxNumberUnknowns];

	bool *pointValid = new bool[m_nNumberObs];

	bool retVal = true;


	if (!ata || !a || !atl || !coe_approximate || !pointValid)
	{
		if (ata) delete[]ata;
		if (a) delete[]a;
		if (atl) delete[]atl;
		if (coe_approximate) delete[]coe_approximate;
		if (pointValid) delete[] pointValid;
		return false;
	}


	double vv = 0, v, com_obs;
	int i, j, k, k0, m;
	int totalPointsNotUsed = 0;
	int pointsNotUsedThisRound = 0;

	// Initialise
	for (i = 0; i < m_nNumberObs; i++)  pointValid[i] = true;

	do
	{

		m_Matrix.Set_Vec_Zero(coe_approximate, m_nNumberCoefficients);

		for (m = 0; m < 2; m++)
		{
			m_Matrix.Set_Vec_Zero(atl, MaxNumberUnknowns);
			m_Matrix.Set_Vec_Zero(ata, MaxNumberUnknowns * MaxNumberUnknowns);

			for (i = 0; i < m_nNumberObs; i++)
			{
				if (pointValid[i])
				{
					ComputeCoeff(m_pdfTime[i], a); // coeffcients of error equation
					com_obs = m_Matrix.Vec_Multiply_Vec(coe_approximate, a, m_nNumberCoefficients);

					v = m_pdfObs[i] - com_obs;   /* constant term of error equation */

					k0 = 0;
					for (j = 0; j < m_nNumberCoefficients; j++)
					{
						atl[j] += v * a[j];
						for (k = 0; k <= j; k++)
						{
							ata[k0] += a[j] * a[k];
							k0++;
						}
					}  // normal equations formation
				}
			}

			// solve normal equations
			m_Matrix.TriMat_Inverse(ata, m_nNumberCoefficients);
			m_Matrix.TriMat_Multiply_Vec(ata, atl, a, m_nNumberCoefficients);

			for (j = 0; j < m_nNumberCoefficients; j++) coe_approximate[j] += a[j];
		}

		for (j = 0; j < m_nNumberCoefficients; j++) m_pdfCoe[j] = coe_approximate[j];

		int redudant_points = m_nNumberObs - m_nNumberCoefficients - totalPointsNotUsed;
		pointsNotUsedThisRound = 0;

		if (redudant_points > 0)
		{
			vv = 0;
			double max_res = 0;
			int max_res_index = -1;

			for (i = 0; i < m_nNumberObs; i++)
			{
				if (pointValid[i])
				{
					ComputeCoeff(m_pdfTime[i], a); // compute coeffcients of chebyshev polynomial

					com_obs = m_Matrix.Vec_Multiply_Vec(m_pdfCoe, a, m_nNumberCoefficients);	// fitted value at time[ i ]

					m_pdfResidual[i] = m_pdfObs[i] - com_obs;
					vv += m_pdfResidual[i] * m_pdfResidual[i];

					if (fabs(m_pdfResidual[i]) > max_res)
					{
						max_res = fabs(m_pdfResidual[i]);
						max_res_index = i;
					}
				}
			}

			m_dfSigma = sqrt(vv / (redudant_points));

			if (max_res > 3.0 * m_dfSigma)
			{
				totalPointsNotUsed++;
				pointsNotUsedThisRound++;
				pointValid[max_res_index] = false;
			}
		}

		if (m_nNumberObs - m_nNumberCoefficients - totalPointsNotUsed < 0)
		{
			retVal = false;
			break;
		}

		/*        double notUsedPercent = ( 100.0 * totalPointsNotUsed ) / m_nNumberObs;
				if ( notUsedPercent > 9.5 )
				{
					retVal = false;
					break;
				}*/

		if (0 == pointsNotUsedThisRound) break;	// no outlier is detected and return successfully

	} while (true);

	delete[]ata;
	delete[]a;
	delete[]atl;
	delete[]coe_approximate;
	delete[]pointValid;

	return retVal;
}


bool cPolyFit::SolveCoefficients( bool *mark )
{
	if( !m_bInit ) 
	{
		return false;
	}

    double *ata = new double[ MaxNumberUnknowns * MaxNumberUnknowns ];
    double *a = new double[ MaxNumberUnknowns ];
	double *atl = new double[ MaxNumberUnknowns ];
	double *coe_approximate = new double[ MaxNumberUnknowns ];

    bool *pointValid = new bool[m_nNumberObs];
    
	bool retVal = true;


	if( !ata || !a || !atl || !coe_approximate || !pointValid )
	{
		if( ata ) delete []ata;
		if( a ) delete []a;
		if( atl ) delete []atl;
		if( coe_approximate ) delete []coe_approximate;
        if ( pointValid ) delete [] pointValid;
		return false;
	}

                                                  
	double vv = 0, v, com_obs;
    int i, j, k, k0, m;
    int totalPointsNotUsed = 0;
    int pointsNotUsedThisRound =0;

    // Initialise
    for( i = 0; i < m_nNumberObs; i++ )  pointValid[ i ] = true;

    do
    {
	    m_Matrix.Set_Vec_Zero( coe_approximate, m_nNumberCoefficients );

	    for( m = 0; m < 2; m++ )
	    {
		    m_Matrix.Set_Vec_Zero( atl, MaxNumberUnknowns );
		    m_Matrix.Set_Vec_Zero( ata, MaxNumberUnknowns * MaxNumberUnknowns );

		    for( i = 0; i < m_nNumberObs; i++ )
		    {
                if ( pointValid[ i ] )
                {
			        ComputeCoeff( m_pdfTime[ i ], a ); // coefficients of error equation
			        com_obs = m_Matrix.Vec_Multiply_Vec( coe_approximate, a, m_nNumberCoefficients );

			        v = m_pdfObs[ i ] - com_obs;   /* constant term of error equation */

			        k0 = 0;
			        for( j = 0; j < m_nNumberCoefficients; j++ )
			        {
				        atl[ j ] += v * a[ j ];
				        for( k = 0; k <= j; k++ )              
				        { 
					        ata[ k0 ] += a[ j ] * a[ k ]; 
					        k0++; 
				        }                    
			        }  // normal equations formation
                }
		    }

		    // solve normal equations
		    m_Matrix.TriMat_Inverse( ata, m_nNumberCoefficients );
		    m_Matrix.TriMat_Multiply_Vec( ata, atl, a, m_nNumberCoefficients );

		    for( j = 0; j < m_nNumberCoefficients; j++ ) coe_approximate[ j ] += a[ j ];
	    }

	    for( j = 0; j < m_nNumberCoefficients; j++ ) m_pdfCoe[ j ] = coe_approximate[ j ];

        int redudant_points = m_nNumberObs - m_nNumberCoefficients - totalPointsNotUsed;
		pointsNotUsedThisRound = 0;

		if( redudant_points > 0 )
		{
			vv = 0;
			double max_res = 0;
			int max_res_index = -1;

			for( i = 0; i < m_nNumberObs; i++ )
			{
				if ( pointValid[i] )
				{
					ComputeCoeff( m_pdfTime[ i ], a ); // compute coefficients of polynomial
        	    
			        com_obs = m_Matrix.Vec_Multiply_Vec( m_pdfCoe, a, m_nNumberCoefficients );	// fitted value at time[ i ]
        		
			        m_pdfResidual[ i ] = m_pdfObs[ i ] - com_obs;
			        vv += m_pdfResidual[ i ] * m_pdfResidual[ i ];

					if( fabs( m_pdfResidual[ i ] ) > max_res )
					{
						max_res = fabs( m_pdfResidual[ i ] );
						max_res_index = i;
					}
				}
			}

		    m_dfSigma = sqrt( vv / ( redudant_points ) );

			if( max_res > 6.0 * m_dfSigma )
			{
				totalPointsNotUsed++;
				pointsNotUsedThisRound++;
				pointValid[ max_res_index ] = false;
			}        
		}

		if( m_nNumberObs - m_nNumberCoefficients - totalPointsNotUsed <= 0 )
		{
            retVal = false;
            break;
		}

/*        double notUsedPercent = ( 100.0 * totalPointsNotUsed ) / m_nNumberObs;
        if ( notUsedPercent > 9.5 )
        {
            retVal = false;
            break;
        }*/

        if ( 0 == pointsNotUsedThisRound ) break;	// no outlier is detected and return successfully

    } while( true ); 

	if( retVal )
	{
		for( int i = 0; i < m_nNumberObs; i++ ) mark[ i ] = pointValid[ i ];
	}

	delete []ata;
	delete []a;
	delete []atl;
	delete []coe_approximate;
    delete []pointValid;

	return retVal;
}

bool cPolyFit::SolveCoefficients(double SigmaThreshold, bool* mark)
{
	if (!m_bInit)
	{
		return false;
	}

	double* ata = new double[MaxNumberUnknowns * MaxNumberUnknowns];
	double* a = new double[MaxNumberUnknowns];
	double* atl = new double[MaxNumberUnknowns];
	double* coe_approximate = new double[MaxNumberUnknowns];

	bool* pointValid = new bool[m_nNumberObs];

	bool retVal = true;


	if (!ata || !a || !atl || !coe_approximate || !pointValid)
	{
		if (ata) delete[]ata;
		if (a) delete[]a;
		if (atl) delete[]atl;
		if (coe_approximate) delete[]coe_approximate;
		if (pointValid) delete[] pointValid;
		return false;
	}


	double vv = 0, v, com_obs;
	int i, j, k, k0, m;
	int totalPointsNotUsed = 0;
	int pointsNotUsedThisRound = 0;

	// Initialise
	for (i = 0; i < m_nNumberObs; i++)  pointValid[i] = true;

	do
	{
		m_Matrix.Set_Vec_Zero(coe_approximate, m_nNumberCoefficients);

		for (m = 0; m < 2; m++)
		{
			m_Matrix.Set_Vec_Zero(atl, MaxNumberUnknowns);
			m_Matrix.Set_Vec_Zero(ata, MaxNumberUnknowns * MaxNumberUnknowns);

			for (i = 0; i < m_nNumberObs; i++)
			{
				if (pointValid[i])
				{
					ComputeCoeff(m_pdfTime[i], a); // coefficients of error equation
					com_obs = m_Matrix.Vec_Multiply_Vec(coe_approximate, a, m_nNumberCoefficients);

					v = m_pdfObs[i] - com_obs;   /* constant term of error equation */

					k0 = 0;
					for (j = 0; j < m_nNumberCoefficients; j++)
					{
						atl[j] += v * a[j];
						for (k = 0; k <= j; k++)
						{
							ata[k0] += a[j] * a[k];
							k0++;
						}
					}  // normal equations formation
				}
			}

			// solve normal equations
			m_Matrix.TriMat_Inverse(ata, m_nNumberCoefficients);
			m_Matrix.TriMat_Multiply_Vec(ata, atl, a, m_nNumberCoefficients);

			for (j = 0; j < m_nNumberCoefficients; j++) coe_approximate[j] += a[j];
		}

		for (j = 0; j < m_nNumberCoefficients; j++) m_pdfCoe[j] = coe_approximate[j];

		int redudant_points = m_nNumberObs - m_nNumberCoefficients - totalPointsNotUsed;
		pointsNotUsedThisRound = 0;

		if (redudant_points > 0)
		{
			vv = 0;
			double max_res = 0;
			int max_res_index = -1;

			for (i = 0; i < m_nNumberObs; i++)
			{
				if (pointValid[i])
				{
					ComputeCoeff(m_pdfTime[i], a); // compute coefficients of polynomial

					com_obs = m_Matrix.Vec_Multiply_Vec(m_pdfCoe, a, m_nNumberCoefficients);	// fitted value at time[ i ]

					m_pdfResidual[i] = m_pdfObs[i] - com_obs;
					vv += m_pdfResidual[i] * m_pdfResidual[i];

					if (fabs(m_pdfResidual[i]) > max_res)
					{
						max_res = fabs(m_pdfResidual[i]);
						max_res_index = i;
					}
				}
			}

			m_dfSigma = sqrt(vv / (redudant_points));

			if (max_res > SigmaThreshold * m_dfSigma)
			{
				totalPointsNotUsed++;
				pointsNotUsedThisRound++;
				pointValid[max_res_index] = false;
			}
		}

		if (m_nNumberObs - m_nNumberCoefficients - totalPointsNotUsed <= 0)
		{
			retVal = false;
			break;
		}

		/*        double notUsedPercent = ( 100.0 * totalPointsNotUsed ) / m_nNumberObs;
				if ( notUsedPercent > 9.5 )
				{
					retVal = false;
					break;
				}*/

		if (0 == pointsNotUsedThisRound) break;	// no outlier is detected and return successfully

	} while (true);

	if (retVal)
	{
		for (int i = 0; i < m_nNumberObs; i++) mark[i] = pointValid[i];
	}

	delete[]ata;
	delete[]a;
	delete[]atl;
	delete[]coe_approximate;
	delete[]pointValid;

	return retVal;
}


/******************************************************************************

  Least squares polynomial fitting of a set of observations and times

******************************************************************************/

bool cPolyFit::SolveCoefficientsAzEl()
{
	if( !m_bInit ) 
	{
		return false;
	}

    double *ata = new double[ MaxNumberUnknowns * MaxNumberUnknowns ];
    double *a = new double[ MaxNumberUnknowns ];
	double *atl = new double[ MaxNumberUnknowns ];
	double *coe_approximate = new double[ MaxNumberUnknowns ];

    bool *pointValid = new bool[m_nNumberObs];
    
	bool retVal = true;

	if( !ata || !a || !atl || !coe_approximate || !pointValid )
	{
		if( ata ) delete []ata;
		if( a ) delete []a;
		if( atl ) delete []atl;
		if( coe_approximate ) delete []coe_approximate;
        if ( pointValid ) delete [] pointValid;
		return false;
	}

                                                  
	double vv = 0, v, com_obs;
    int i, j, k, k0, m;
    int totalPointsNotUsed = 0;
    int pointsNotUsedThisRound =0;

	m_nNumberCoefficients = 7;

    // Initialise
    for( i = 0; i < m_nNumberObs; i++ )  pointValid[ i ] = true;

    do
    {

	    m_Matrix.Set_Vec_Zero( coe_approximate, m_nNumberCoefficients );

	    for( m = 0; m < 2; m++ )
	    {
		    m_Matrix.Set_Vec_Zero( atl, MaxNumberUnknowns );
		    m_Matrix.Set_Vec_Zero( ata, MaxNumberUnknowns * MaxNumberUnknowns );

		    for( i = 0; i < m_nNumberObs; i++ )
		    {
                if ( pointValid[ i ] )
                {
			        ComputeCoeffAzEl( m_pdfTime[ i ], a ); // coeffcients of error equation
			        com_obs = m_Matrix.Vec_Multiply_Vec( coe_approximate, a, m_nNumberCoefficients );

			        v = m_pdfObs[ i ] - com_obs;   /* constant term of error equation */

			        k0 = 0;
			        for( j = 0; j < m_nNumberCoefficients; j++ )
			        {
				        atl[ j ] += v * a[ j ];
				        for( k = 0; k <= j; k++ )              
				        { 
					        ata[ k0 ] += a[ j ] * a[ k ]; 
					        k0++; 
				        }                    
			        }  // normal equations formation
                }
		    }

		    // solve normal equations
		    m_Matrix.TriMat_Inverse( ata, m_nNumberCoefficients );
		    m_Matrix.TriMat_Multiply_Vec( ata, atl, a, m_nNumberCoefficients );

		    for( j = 0; j < m_nNumberCoefficients; j++ ) coe_approximate[ j ] += a[ j ];
	    }

	    for( j = 0; j < m_nNumberCoefficients; j++ ) m_pdfCoe[ j ] = coe_approximate[ j ];

        int redudant_points = m_nNumberObs - m_nNumberCoefficients - totalPointsNotUsed;
		pointsNotUsedThisRound = 0;

		if( redudant_points > 0 )
		{
			vv = 0;
			double max_res = 0;
			int max_res_index = -1;

			for( i = 0; i < m_nNumberObs; i++ )
			{
				if ( pointValid[i] )
				{
					ComputeCoeffAzEl( m_pdfTime[ i ], a ); // compute coeffcients of chebyshev polynomial
        	    
			        com_obs = m_Matrix.Vec_Multiply_Vec( m_pdfCoe, a, m_nNumberCoefficients );	// fitted value at time[ i ]
        		
			        m_pdfResidual[ i ] = m_pdfObs[ i ] - com_obs;
			        vv += m_pdfResidual[ i ] * m_pdfResidual[ i ];

					if( fabs( m_pdfResidual[ i ] ) > max_res )
					{
						max_res = fabs( m_pdfResidual[ i ] );
						max_res_index = i;
					}
				}
			}

		    m_dfSigma = sqrt( vv / ( redudant_points ) );

			if( max_res > 3.0 * m_dfSigma )
			{
				totalPointsNotUsed++;
				pointsNotUsedThisRound++;
				pointValid[ max_res_index ] = false;
			}        
		}

		if( m_nNumberObs - m_nNumberCoefficients - totalPointsNotUsed <= 0 )
		{
            retVal = false;
            break;
		}

/*        double notUsedPercent = ( 100.0 * totalPointsNotUsed ) / m_nNumberObs;
        if ( notUsedPercent > 9.5 )
        {
            retVal = false;
            break;
        }*/

        if ( 0 == pointsNotUsedThisRound ) break;	// no outlier is detected and return successfully

    } while( true ); 

	delete []ata;
	delete []a;
	delete []atl;
	delete []coe_approximate;
    delete []pointValid;

	return retVal;
}


/******************************************************************************

  Get the polynomial coefficients 

******************************************************************************/

void cPolyFit::GetCoefficients( int &NumCoefficients, double *coeff )
{
    NumCoefficients = m_nNumberCoefficients; 
     
    m_Matrix.Vec_Assign( m_pdfCoe, coeff, m_nNumberCoefficients );
}

/******************************************************************************

  Set the coefficients of a polynomial.
  
  numberCoefficients: number of coefficients in the array pdfCoeffs. The number
                      of coefficients is the order of the polynomial plus 1.

  pdfCoeffs:

  After calling this function, one then can call ComputeFittedValueAtGivenTime() 
  and ComputeDotAtTimeT() to compute the value of the polynomial at the given 
  time as well as the first derivative of teh polynomial at the given time.

******************************************************************************/

void cPolyFit::SetCoefficients(int numberCoefficients, double *pdfCoeffs )
{
    m_nNumberCoefficients = numberCoefficients; 
     
    m_Matrix.Vec_Assign( pdfCoeffs, m_pdfCoe, m_nNumberCoefficients );
}

/******************************************************************************

  Get the RMS of the observation residuals after the fitting

******************************************************************************/

double cPolyFit::GetSigma()
{
	return m_dfSigma;
}

/******************************************************************************

  get the observation residuals after the polynomial fitting 

******************************************************************************/

void cPolyFit::GetResiduls( double *residual )
{
	m_Matrix.Vec_Assign( m_pdfResidual, residual, m_nNumberObs );
}

/******************************************************************************

  compute the fitted value at the given epoch from the fitted polynomial

******************************************************************************/

double cPolyFit::ComputeFittedValueAtGivenTime( double t )
{
    double *a = new double[ MaxNumberUnknowns ];

	if( !a ) return 1.0e20;
	
	ComputeCoeff( t, a ); 

	double p = m_Matrix.Vec_Multiply_Vec( a, m_pdfCoe, m_nNumberCoefficients );

	delete []a;

	return p;
}

/******************************************************************************

  compute the fitted value at the given epoch from the fitted polynomial

******************************************************************************/

double cPolyFit::ComputeFittedAzElValueAtGivenTime( double t )
{
    double *a = new double[ MaxNumberUnknowns ];

	if( !a ) return 1.0e20;
	
	ComputeCoeffAzEl( t, a ); 

	double p = m_Matrix.Vec_Multiply_Vec( a, m_pdfCoe, m_nNumberCoefficients );

	delete []a;

	return p;
}


/******************************************************************************

  compute the velocity at the given epoch from the fitted polynomial

******************************************************************************/

double cPolyFit::ComputeDotAtTimeT( double t )
{
	double deltaT = 0.001;

	double a1 = ComputeFittedValueAtGivenTime( t );
	double a2 = ComputeFittedValueAtGivenTime( t + deltaT );

	double dot = ( a2 - a1 ) / deltaT;

	return dot;
}


