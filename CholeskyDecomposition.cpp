/***************************************************************************************

 Cholesky Decomposition


 Reference: Section 2.9, pp. 89-91, Numerical Recipes in Fortran, William Press, et al.

***************************************************************************************/

#include <math.h>
#include "CholeskyDecomposition.h"
#include <vector>

using namespace std;

cCholeskyDecom::cCholeskyDecom()
{

} 

cCholeskyDecom::~cCholeskyDecom()
{

}


/***************************************************************************************

 See P. 90, Numerical Recipes in Fortran, William Press, et al.
 
 A = L x L(T)
 A is a symmetric and positive definite square matrix of dimension nDim, then
 L is a low triangular matrix. This is the Cholesky Decomposition (CD).

 The CD method is extremely stable numerically. Failure of CD simply indicates
 that A (possible with roundoff error) is not positive definite
 
 Here, pdfMatrix is the low triangluar part of A as input and L as output

 Fortran Subroutine: cholesky_decomposition

***************************************************************************************/

BOOL cCholeskyDecom::Decompose( double *pdfMatrix, int nDim )
{
	if( !pdfMatrix )
	{
		return FALSE;
	}

	try
	{
		double sum;
		int i, ii, j, ji, k, ik, jk;

		k = nDim * ( nDim + 1 ) / 2;
		vector <double> v_cholesky( 10 );

        if ( (int)v_cholesky.size() < k )
        {
            v_cholesky.resize( k );
        }

		for( i = 0; i < k; i++ )  v_cholesky[ i ] = 0.0;

		for( i = 0; i < nDim; i++ )
		{
			ii = ( i + 1 ) * ( i + 2 ) / 2 - 1;

			for( j = i; j < nDim; j++ )
			{
				ji = j * ( j + 1 ) / 2 + i;
				sum = pdfMatrix[ ji ];

				if( i > 0 )
				{
					for( k = i - 1; k >= 0; k-- )
					{
						ik = i * ( i + 1 ) / 2 + k;
						jk = j * ( j + 1 ) / 2 + k;

						sum = sum - v_cholesky[ ik ] * v_cholesky[ jk ];
					}
				}

				if( i == j )
				{
					if( sum <= 0.0 )
					{
						return FALSE;
					}

					v_cholesky[ ii ] = sqrt( sum );
				}
				else
				{
					v_cholesky[ ji ] = sum / v_cholesky[ ii ];
				}
			}
		}

		k = nDim * ( nDim + 1 ) / 2;
		for( i = 0; i < k; i++ ) pdfMatrix[ i ] = v_cholesky[ i ];

		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}
}

BOOL cCholeskyDecom::Decompose( double *pdfMatrix, int nDim, int nDim0 )
{
	if( !pdfMatrix )
	{
		return FALSE;
	}

	try
	{
		double sum;
		int i, ii, j, ji, k, ik, jk;

		double *v = new double[ nDim0 ];
		if( !v ) 
		{
			return FALSE;
		}

		for( i = 0; i < nDim0; i++ )  v[ i ] = 0.0;

		for( i = 0; i < nDim; i++ )
		{
			ii = ( i + 1 ) * ( i + 2 ) / 2 - 1;

			for( j = i; j < nDim; j++ )
			{
				ji = j * ( j + 1 ) / 2 + i;
				sum = pdfMatrix[ ji ];

				if( i > 0 )
				{
					for( k = i - 1; k >= 0; k-- )
					{
						ik = i * ( i + 1 ) / 2 + k;
						jk = j * ( j + 1 ) / 2 + k;

						sum = sum - v[ ik ] * v[ jk ];
					}
				}

				if( i == j )
				{
					if( sum <= 0.0 )
					{
						delete []v;
						return FALSE;
					}

					v[ ii ] = sqrt( sum );
				}
				else
				{
					v[ ji ] = sum / v[ ii ];
				}
			}
		}

		for( i = 0; i < nDim0; i++ ) pdfMatrix[ i ] = v[ i ];

		delete []v;
		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}
}

/***************************************************************************************

 use Cholesky decomposition to solve linear equations
 A x = b
 A is a symmetric and positive definite square matrix.
 
 because A = L L(T)

 A x = b -->  L L(T) x = b --> L y = b, so that y is determined
 then y = L(T) x, x is determined. Because L is a low triangular matrix,
 solving y and x is very easy.

 In this software, pdfMatrix is the low triangular part of A

 Fortran Subroutine: cholesky_decomposition_solution

***************************************************************************************/

BOOL cCholeskyDecom::LinearEquation( double *pdfMatrix, double *pdfB, 
									 int nDim, double *pdfX )
{
	if( !pdfMatrix || !pdfB || !pdfX )
	{
		return FALSE;
	}

	if( !Decompose( pdfMatrix, nDim, nDim * ( nDim + 1 ) / 2 ) ) 
	{
		return FALSE;
	}

	try
	{
		int i, ii, k, ik, ki;
		double sum;

		for( i = 0; i < nDim; i++ )	// L y = b
		{
			ii = ( i + 1 ) * ( i + 2 ) / 2 - 1;
			sum = pdfB[ i ];

			if( i > 0 ) 
			{
				for( k = i - 1; k >= 0; k-- )
				{
					ik = i * ( i + 1 ) / 2 + k;
					sum = sum - pdfMatrix[ ik ] * pdfX[ k ];
				}
			}

			pdfX[ i ] = sum / pdfMatrix[ ii ];
		}

		for( i = nDim - 1; i >= 0; i-- )	// L(T) x = y
		{
			ii = ( i + 1 ) * ( i + 2 ) / 2 - 1;
			sum = pdfX[ i ];

			if( i < nDim - 1 )
			{
				for( k = i + 1; k < nDim; k++ )
				{
					ki = k * ( k + 1 ) / 2 + i;
					sum = sum - pdfMatrix[ ki ] * pdfX[ k ];
				}
			}
			pdfX[ i ] = sum / pdfMatrix[ ii ];
		}

		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}
}

/***************************************************************************************

 inversion of a low triangular matrix L
 
   L = * 0 0 0 0 0 0 0 
       * * 0 0 0 0 0 0 
       * * * 0 0 0 0 0
       * * * * 0 0 0 0
       * * * * * 0 0 0
       * * * * * * 0 0
       * * * * * * * 0
       * * * * * * * *

 It is easily seen that the inversion of a low triangular matrix is still a
 low triangular matrix. 

 Here, pdfMatrix is the **** part of L as input, and the *** part of L(-1)

 Fortran Subroutine: inversion_cholesky_low_triangular_matrix

***************************************************************************************/

BOOL cCholeskyDecom::InverseLowTriangular( double *pdfMatrix, int nDim )
{
	if( !pdfMatrix )
	{
		return FALSE;
	}

	try
	{
		int i, ii, j, ji, k, jk, ik;
		double c;

		for( i = 0; i < nDim - 1; i++ )
		{
			ii = ( i + 1 ) * ( i + 2 ) / 2 - 1;

			if( fabs( pdfMatrix[ ii ] ) < 1.0e-50 )
			{
				return FALSE;
			}

			for( j = i + 1; j < nDim; j++ )
			{
				ji = j * ( j + 1 ) / 2  + i;
				
				c = -pdfMatrix[ ji ] / pdfMatrix[ ii ];

				for( k = 0; k <= i; k++ )
				{
					jk = j * ( j + 1 ) / 2 + k;
					ik = i * ( i + 1 ) / 2 + k;

					if( k == i ) pdfMatrix[ jk ] = c;
					else pdfMatrix[ jk ] += c * pdfMatrix[ ik ];

				}
			}
		}

		for( i = 0; i < nDim; i++ )
		{
			ii = ( i + 1 ) * ( i + 2 ) / 2 - 1;
			c = pdfMatrix[ ii ];
			pdfMatrix[ ii ] = 1.0 / c;

			ii = i * ( i + 1 ) / 2;
			
			for( j = 0; j < i; j++ )
			{
				pdfMatrix[ ii + j ] /= c;
			}
		}

		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}
}

/***************************************************************************************

 Compute (L(-1))(T) * L(-1), where L(-1) is a low triangular matrix
 
 This method is for computing the inversion of a symmetric and positive
 definite square matrix using the Cholesky Decomposition. Because 
 A = L * L(T), so A(-1) = ( L(-1)(T) ) * ( L(-1) )


 here as input, pdfMatrix is the low triangular part of A if nCase = 1
                             L                               nCase = 2
                             L(-1)                           nCase = 3

 as output, pdfMatrix is the low triangular part of A(-1)

 Fortran Subroutine: inversion_cholesky_decomposition

 nCase: 1 when the input matrix is the origianl coefficients
        2 when the input matrix is the decomposed low triangle matrix of the original matrix

***************************************************************************************/

BOOL cCholeskyDecom::MatrixInversing( double *pdfMatrix, int nDim, int nCase )
{
	if( !pdfMatrix )
	{
		return FALSE;
	}

	switch( nCase )
	{
		case 1:
			if( !Decompose( pdfMatrix, nDim ) ) 
			{
				return FALSE;
			}
			
			if( !InverseLowTriangular( pdfMatrix, nDim ) ) 
			{
				return FALSE;
			}
			
			break;

		case 2:
			
			if( !InverseLowTriangular( pdfMatrix, nDim ) ) 
			{
				return FALSE;
			}
			
			break;
	}

	try
	{
		// the following computes ( L(-1)(T) ) * ( L(-1) )
		int i, j, k, ik, kj, ij;
		double sum;

		k = nDim * ( nDim + 1 ) / 2;

		vector <double> v_temp( 10 );

        if ( (int)v_temp.size() < k )
        {
            v_temp.resize( k );
        }

		for( i = 0; i < k; i++ ) v_temp[ i ] = 0.0;

		for( i = 0; i < nDim; i++ )
		{
			for( j = 0; j <= i; j++ )
			{
				sum = 0.0;

				for( k = i; k < nDim; k++ )
				{
					ik = k * ( k + 1 ) / 2 + i;
					kj = k * ( k + 1 ) / 2 + j;
					sum += pdfMatrix[ ik ] * pdfMatrix[ kj ];
				}
				
				ij = i * ( i + 1 ) / 2 + j;
				v_temp[ ij ] = sum;
			}
		}

		k = nDim * ( nDim + 1 ) / 2;

		for( i = 0; i < k; i++ ) pdfMatrix[ i ] = v_temp[ i ];

		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}
}


BOOL cCholeskyDecom::MatrixInversing( double *pdfMatrix, int nDim, int nDim0, int nCase )
{
	if( !pdfMatrix )
	{
		return FALSE;
	}

	switch( nCase )
	{
		case 1:
			if( !Decompose( pdfMatrix, nDim, nDim0 ) ) 
			{
				return FALSE;
			}
			
			if( !InverseLowTriangular( pdfMatrix, nDim ) ) 
			{
				return FALSE;
			}
			
			break;

		case 2:
			
			if( !InverseLowTriangular( pdfMatrix, nDim ) ) 
			{
				return FALSE;
			}
			
			break;
	}

	try
	{
		// the following computes ( L(-1)(T) ) * ( L(-1) )
		int i, j, k, ik, kj, ij;
		double sum;

		double *v = new double[ nDim0 ];
		if( !v ) 
		{
			return FALSE;
		}

		for( i = 0; i < nDim0; i++ ) v[ i ] = 0.0;

		for( i = 0; i < nDim; i++ )
		{
			for( j = 0; j <= i; j++ )
			{
				sum = 0.0;

				for( k = i; k < nDim; k++ )
				{
					ik = k * ( k + 1 ) / 2 + i;
					kj = k * ( k + 1 ) / 2 + j;
					sum += pdfMatrix[ ik ] * pdfMatrix[ kj ];
				}
				
				ij = i * ( i + 1 ) / 2 + j;
				v[ ij ] = sum;
			}
		}

		for( i = 0; i < nDim0; i++ ) pdfMatrix[ i ] = v[ i ];

		delete []v;

		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}
}
