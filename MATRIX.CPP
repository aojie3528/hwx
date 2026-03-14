/*****************************************************************************

  Matrix computation implementation.

******************************************************************************/

//----------------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------------

#include <math.h>
#include "Matrix.h"

#include "CholeskyDecomposition.h"

class cMatrix g_Matrix;

/****************************************************************************

  constructor

****************************************************************************/

cMatrix::cMatrix()
{
	;
}


cMatrix::~cMatrix()
{
	;
}


/****************************************************************************

  copy a vector of double into another 

****************************************************************************/

void cMatrix::Vec_Assign( double *Source, double *Result, int n )
{
	if( !Source || !Result )
	{
		return;
	}

	for( int i = 0; i < n; i++ ) Result[ i ] = Source[ i ];
}


/****************************************************************************

  set elements of a vector of double to zero 

****************************************************************************/

void cMatrix::Set_Vec_Zero( double *Vec, int n )
{
	if( !Vec )
	{
		return;
	}

	for ( int i = 0; i < n; i++ ) Vec[ i ] = 0;
}


/****************************************************************************

  multiply a vector of double by a constant 

****************************************************************************/

void cMatrix::Vec_Multiply_By_Constant( double *Vec, int n, double c )
{
	if( !Vec )
	{
		return;
	}

	for ( int i = 0; i < n; i++ ) Vec [ i ] *= c;
}

/****************************************************************************

  multiply a triangulised matrix of double by a constant 

****************************************************************************/

void cMatrix::TriMat_Multiply_By_Constant( double *TriMat, int n, 
										   double c )
{
	if( !TriMat )
	{
		return;
	}

	int i = 0, k = n * (n+1) / 2;
    
    for( i = 0; i < k; i++ ) TriMat[ i ] *= c;
}

/****************************************************************************

  multiply a matrix of double by a constant 

****************************************************************************/

void cMatrix::Mat_Multiply_By_Constant( double *Mat, int m, int n, 
									    double c )
{
	if( !Mat )
	{
		return;
	}

    int i = 0, k = m * n;
    
    for( i = 0; i < k; i++ ) Mat[ i ] *= c;
}

/****************************************************************************

  matrix transpose 

  As input, Mat is a vectorized matrix by line
  As output, Mat is the transpose of the input, vectorized by line

****************************************************************************/

BOOL cMatrix::Mat_Transpose( double *Mat, int m, int n )
{
	if( !Mat )
	{
		return FALSE;
	}

	double *ttt = new double[ m * n ];
	int i, j, ij, ji;

	if( !ttt )
	{
		return FALSE;
	}

	try
	{
		for( i = 0; i < m * n; i++ ) ttt[ i ] = 0;

		for( i = 0; i < m; i++ )
		for( j = 0; j < n; j++ )
		{    
			ij = i * n + j;   
			ji = j * m + i;   
			ttt[ ji ] = Mat[ ij ];    
		}

		ij = m * n;

		for( i = 0; i < ij; i++ ) Mat[ i ] = ttt[ i ];

		delete []ttt;

		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}
} 

/****************************************************************************

  transpose of an upper triangular matrix of order n
  
  As input, TriMat is the upper triangular matrix vectorized by line in the 
  form of ( 0s are not included)

     * * * * *
	 0 * * * *
     0 0 * * *
	 0 0 0 * *
	 0 0 0 0 *

  As output, TriMat becomes a lower triangular matrix vectorized by line
  ( 0s are not included)

     * 0 0 0 0
	 * * 0 0 0
	 * * * 0 0
	 * * * * 0
	 * * * * *

****************************************************************************/

BOOL cMatrix::UpperTriMat_Transpose( double *TriMat, int n )
{
	double *a = new double[ n * ( n + 1 ) / 2 ];

	if( !a )
	{
		return FALSE;
	}

	int l = 0;

	for( int i = 1; i <= n; i++ )	// the column
	{
		for( int j = i; j <= n; j++ )	// the row
		{
			int ji = ( j - 1 ) * j / 2 + i - 1;
			a[ ji ] = TriMat[ l ];
			l++;
		}
	}

	Vec_Assign( a, TriMat, n * ( n + 1 ) / 2 );

	delete []a;

	return TRUE;
}

/****************************************************************************

  transpose of a lower triangular matrix of order n
  
  As input, TriMat is the lower triangular matrix vectorized by line 
  ( 0s are not included)

     * 0 0 0 0
	 * * 0 0 0
	 * * * 0 0
	 * * * * 0
	 * * * * *

  As output, TriMat becomes an upper triangular matrix vectorized by line
  ( 0s are not included)

     * * * * *
	 0 * * * *
     0 0 * * *
	 0 0 0 * *
	 0 0 0 0 *

	 ****************************************************************************/

BOOL cMatrix::LowerTriMat_Transpose( double *TriMat, int n )
{
	double *a = new double[ n * ( n + 1 ) / 2 ];

	if( !a )
	{
		return FALSE;
	}

	int l = 0;

	for( int i = 1; i <= n; i++ )	// the column
	{
		for( int j = i; j <= n; j++ )	// the row
		{
			int ji = ( j - 1 ) * j / 2 + i - 1;
			a[ l ] = TriMat[ ji ];
			l++;
		}
	}

	Vec_Assign( a, TriMat, n * ( n + 1 ) / 2 );

	delete []a;

	return TRUE;
}


/****************************************************************************

  inversing of a symmetric matrix
  
  As input, TriMat the low triangle of a symmetric matrix of dimension n, 
            vectorized by line
  As output, TriMat the low triangle of the inverse of input, 
            vectorized by line

****************************************************************************/

BOOL cMatrix::TriMat_Inverse( double *TriMat, int n )
{
    double c, bs, bf, bh, *b;
	double *ea;
    int i, k, j, ii;
	
	b = new double[ n ];

	if( !b )
	{
		return FALSE;
	}

	try
	{
		Set_Vec_Zero( b, n );

		ea = TriMat;

		for( k = 1; k < n; k++ )
		{
			ii = k * ( k + 1 ) / 2 - 1;
			c = *( ea + ii );

			for ( j = k + 1; j < n + 1; j++ )
			{
				ii = j * ( j - 1 ) / 2 + k - 1;
				bf = *( ea + ii ) / c;
				
				if( fabs( bf) > 0.0 )
				{
					for( i = j; i < n + 1; i++ )
					{
						ii = i * ( i - 1 ) / 2 + k - 1;
						bh = *( ea + ii );
						ii = ii + j - k;
						*( ea + ii ) = *( ea + ii ) - bh * bf;
					}
				}
			}   
		}

		for( k = n; k > 0; k-- )
		{
			ii = k * ( k + 1 ) / 2 - 1;
			c = *( ea + ii );
	        
			for( i = n; i > k - 1; i-- )
			{
				if( i != k ) 
				{
					ii = i * ( i - 1 ) / 2 + k - 1;
					*( b + i - 1 ) = *( ea + ii );
				}
			}

			for( j = n; j > k - 1; j-- )
			{
				bs = 0.0;
				if( k == j ) bs = 1.0;
				
				for( i = k + 1; i < n + 1; i++ )
				{
					if( i == n + 1 ) goto RET5;
					ii = j * ( j - 1 ) / 2 + i - 1;
					if( i > j ) ii = i * ( i - 1 ) / 2 + j - 1;
					if( i <= n) bs = bs - *( ea + ii ) * *( b + i - 1 );
				}
	RET5:
				bs = bs / c;
				ii = j * ( j - 1 ) / 2 + k - 1 ;
				*( ea + ii ) = bs;
			}
		}

		delete []b;

		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}
}

/****************************************************************************

  inversing of a square matrix
  
  As input, Mat the square matrix of dimension n, vectorized by line
  As output, Mat the inverse of the original input matrix

  n: the order of the square matrix

****************************************************************************/

BOOL cMatrix::Mat_Inverse( double *Mat, int n )
{
	double *b = new double [  n * n ];

	if( !b )
	{
		return FALSE;
	}

	try
	{
		Set_Vec_Zero( b, n * n );

		int i;
		for( i = 0; i < n; i++ ) b[ i * n + i ] = 1.0;

		for( i = 0; i < n - 1; i++ )
		{
			int ii = i * n + i;
			double c = Mat[ ii ];

			for ( int j = i + 1; j < n; j++ )
			{
				int ij = j * n + i;
				double bf = Mat[ ij ] / c;

				for( int k = 0; k < n; k++ )
				{
					ii = i * n + k;
					int ik = j * n + k;
					
					Mat[ ik ] -= Mat[ ii ] * bf;
					b[ ik ] -= b[ ii ] * bf;					
				}				
			}   
		}

		for( i = 0; i < n; i++ )
		{
			int ii = i * n + i;
			double c = Mat[ ii ];

			for( int j = 0; j < n; j++ )
			{
				int ij = i * n + j;
				Mat[ ij ] /= c;
				b[ ij ] /= c;
			}
		}

		for( i = n - 1; i > 0; i-- )
		{
			int ii = i * n + i;

			for ( int j = 0; j < i; j++ )
			{
				int ij = j * n + i;
				double bf = Mat[ ij ];

				for( int k = 0; k < n; k++ )
				{
					ii = i * n + k;
					int ik = j * n + k;
					
					b[ ik ] -= b[ ii ] * bf;					
				}				
			}   
		}

		Vec_Assign( b, Mat, n * n );

		delete []b;

		return TRUE;
	}
	catch( ... )
	{
		return FALSE;
	}	
}

/****************************************************************************

  dot product of two vectors 

****************************************************************************/

double cMatrix::Vec_Multiply_Vec( double *Vec1, double *Vec2, int n )
{
	if( !Vec1 || !Vec2 )
	{
		return 1.0e20;
	}

    double res = 0;
  
    for( int i = 0; i < n; i++ ) res = res + Vec1[ i ] * Vec2[ i ];

    return res;
}

/****************************************************************************

  multiply a vector with another to form a matrix 

****************************************************************************/

void cMatrix::Vec_Multiply_Vec_To_Mat( double *Vec1, double *Vec2, 
									   double *Mat, int n1, int n2 )
{
	if( !Vec1 || !Vec2 || !Mat )
	{
		return;
	}

	try
	{
		int i, j, k = 0;

		for( i = 0; i < n1; i++ )
		{
			for( j = 0; j < n2; j++ )
			{
				Mat[ k ]= Vec1[ i ] * Vec2[ j ];
				k++;
			}
		}
		
		return;
	}
	catch( ... )
	{
	}
}


/****************************************************************************

  multiply a vector with itself to form a symmetric matrix 

  TriMat is the low triangle of a symmetric matrix of dimension n, 
  vectorized by line

****************************************************************************/

void cMatrix::Vec_Multiply_Vec_To_TriMat( double *Vec, 
										  double *TriMat, int n )
{
	if( !Vec || !TriMat )
	{
		return;
	}

	try
	{
		int i, j, k = 0;

		for( i = 0; i < n; i++ )
		{
			for( j = 0; j <= i; j++ )
			{
				TriMat[ k ]= Vec[ i ] * Vec[ j ];
				k++;
			}
		}
		return;
	}
	catch( ... )
	{
	}
}

/****************************************************************************

  multiply a vector with a symmetric matrix 

  TriMat is the low triangle of a symmetric matrix of dimension n, 
  vectorized by line


****************************************************************************/

void cMatrix::Vec_Multiply_TriMat( double *Vec, double *TriMat, 
								   double *Result, int n )
{
	if( !Vec || !TriMat || !Result )
	{
		return;
	}

	try
	{
		int i, j, ii, ij;
		double res;

		for( i = 0; i < n; i++ )
		{
			res = 0; 
			ii = i * ( i + 1 ) / 2;
			
			for( j = 0; j < n; j++ )
			{
				if( j <= i ) ij = ii + j;
				else 
					ij = j * ( j + 1 ) / 2 + i;

				res = res + Vec[ j ] * TriMat[ ij ];
			}

			Result[ i ] = res;
		}
		return;
	}
	catch( ... )
	{
	}
}


/****************************************************************************

  multiply a vector with a matrix 

  Vec is a vector of dimenaion n
  Mat is a matrix of dimension n * m, vectorized by line

****************************************************************************/

void cMatrix::Vec_Multiply_Mat( double *Vec, double *Mat, double *Result, 
							    int n, int m )
{
	if( !Vec || !Mat || !Result )
	{
		return;
	}

	try
	{
		int i, j, ij;
		double res;

		for( i = 0; i < m; i++ )
		{
			res = 0;
	        
			for( j = 0; j < n; j++ )
			{
				ij = j * m + i; 
				res = res + Vec[ j ] * Mat[ ij ];
			}
			
			Result[ i ] = res;
		}
		return;
	}
	catch( ... )
	{
	}
}


/****************************************************************************

  multiply a symmetric matrix with a vector

  TriMat is the low triangle of a symmetric matrix of dimension n, 
  vectorized by line

****************************************************************************/

void cMatrix::TriMat_Multiply_Vec( double *TriMat, double *Vec, 
								   double *Result, int n )
{
	if( !Vec || !TriMat || !Result )
	{
		return;
	}

	Vec_Multiply_TriMat( Vec, TriMat, Result, n );
}


/****************************************************************************

  multiply a matrix with a vector

  Mat is a matrix of dimension m * n, vectorized by line
  Vec is a vector of dimenaion n

****************************************************************************/

void cMatrix::Mat_Multiply_Vec( double *Mat, double *Vec, 
							    double *Result, int m, int n )
{
	if( !Vec || !Mat || !Result )
	{
		return;
	}

	try
	{
		int i, j, ij;
		double res;

		for( i = 0; i < m; i++ )
		{
    		res = 0;
			
			for( j = 0; j < n; j++ )
			{
				ij = i * n + j;
				res = res + Vec[ j ] * Mat[ ij ];
			}

			Result[ i ] = res;
		}
		return;
	}
	catch( ... )
	{
	}
}


/****************************************************************************

  multiply a symmetric matrix with another

  TriMat1 and TriMat2 are the low triangles of two symmetric matries of 
  dimension n, vectorized by line

  Result is a matrix of dimension n, vectorized by line

****************************************************************************/

void cMatrix::TriMat_Multiply_TriMat( double *TriMat1, double *TriMat2, 
									  double *Result, int n )
{
	if( !TriMat1 || !TriMat2 || !Result )
	{
		return;
	}

	try
	{
		int i, j, ii;
		double res[ 1000 ], vec[ 1000 ];
		
		for( i = 0; i < n; i++ )
		{
   			ii = i * ( i + 1 ) / 2;
			
			for( j = 0; j <= i; j++ )  vec[ j ]= TriMat1 [ ii + j ];
			
			if( i + 1 < n )
			{
				for( j = i + 1; j < n; j++ )  
				vec[ j ]= TriMat1[ j * ( j + 1 ) / 2 + i ];
			}

			Vec_Multiply_TriMat( vec, TriMat2, res, n );

			ii = i * n;
			
			for( j = 0; j < n; j++ ) Result[ ii + j ] =  res[ j ];
		}
		return;
	}
	catch( ... )
	{
	}
}


/****************************************************************************

  multiply a symmetric matrix with a matrix
  
  TriMat is the low triangle of a symmetric matrix of dimension n, 
  vectorized by line

  Mat is a matrix of dimension  n * m, vectorized by line

  Result is a mtrix of n * m, vectorized by line

****************************************************************************/

void cMatrix::TriMat_Multiply_Mat( double *TriMat, double *Mat, 
								   double *Result, int n, int m )
{
	if( !TriMat || !Mat || !Result )
	{
		return;
	}

	try
	{
		int i, j, ij;
		double vec[ 1000 ], res[ 1000 ];

		for( i = 0; i < m; i++ )
		{
			for( j = 0; j < n; j++ )
			{
				ij = j * m + i;
				vec[ j ] = Mat[ ij ];
			}

			TriMat_Multiply_Vec( TriMat, vec, res, n );

			for( j = 0; j < n; j++ )
			{
				ij = j * m + i;
				Result[ ij ]= res[ j ];
			}
		}
		return;
	}
	catch( ... )
	{
	}
}


/****************************************************************************

  multiply a  matrix with a symmetric matrix

  TriMat is the low triangle of a symmetric matrix of dimension n, 
  vectorized by line

  Mat is a matrix of dimension  m * n, vectorized by line

  Result is a mtrix of m * n, vectorized by line

****************************************************************************/

void cMatrix::Mat_Multiply_TriMat( double *Mat, double *TriMat, 
								   double *Result, int m, int n )
{
	if( !TriMat || !Mat || !Result )
	{
		return;
	}

	try
	{
		int i, j, ij;
		double vec[ 1000 ], res[ 1000 ];
			
		for( i = 0; i < m; i++ )
		{
			for( j = 0; j < n; j++ )
			{
				ij = i * n + j;
				vec[ j ]= Mat[ ij ];
			}

			Vec_Multiply_TriMat( vec, TriMat, res, n );

			for( j = 0; j < n; j++ )
			{
				ij = i * n + j;
				Result [ ij ] = res[ j ];
			}
		}
		return;
	}
	catch( ... )
	{
	}
}


/****************************************************************************

  multiply a  matrix with another 

  Mat is a matrix of dimension  m * n, vectorized by line
  Mat is a matrix of dimension  n * p, vectorized by line

  Result is a matrix of dimension  m * p, vectorized by line

****************************************************************************/

void cMatrix::Mat_Multiply_Mat( double *Mat1, double *Mat2, 
							    double *Result, int m, int n, int p )
{
	if( !Mat1 || !Mat2 || !Result )
	{
		return;
	}

	try
	{
		double tem1; 
		int k1, k2, k3, l1, l2;

		for ( k1 = 1 ; k1 <= m ; k1++ )
		{
			for ( k2 = 1 ; k2 <= p ; k2++ )
			{
				tem1 = 0;
			
				for ( k3 = 1 ; k3 <= n ; k3++ )
				{
					l1 = (k1-1) * n + k3 -1;
					l2 = (k3-1) * p + k2 -1;
					tem1 = tem1 + Mat1[ l1 ] * Mat2[ l2 ];
				}
			
				l1 = (k1-1) * p + k2 -1;
				Result[ l1 ] = tem1;
			}
		}	
		return;
	}
	catch( ... )
	{
	}
}

/***************************************************************************************

  multiplication of two lower triangular matrices of order n

  two lower triangular matrices are vectorized by line, respectively (0s are not included)

     * 0 0 0 0
	 * * 0 0 0
	 * * * 0 0
	 * * * * 0
	 * * * * *

  the result matrix is also a lower triangular matrix with the same storage arrangement

***************************************************************************************/
BOOL cMatrix::LowerTriMat_Multiply_LowerTriMat( double *LowerTriMat1, 
		                                        double *LowerTriMat2, 
										        double *LowerTriResult, 
		                                        int n )
{
	if( !LowerTriMat1 || !LowerTriMat2 || !LowerTriResult )
	{
		return FALSE;
	}

	int l = 0;
	for( int i = 1; i <= n; i++ ) // the row
	{
		for( int j = 1; j <= i; j++ )	// the column
		{
			LowerTriResult[ l ] = 0;

			int ii = ( i - 1 ) * i / 2 - 1;
			for( int k = j; k <= i; k++ )
			{
				int kj = ( k - 1 ) * k / 2 + j - 1;
				LowerTriResult[ l ] += LowerTriMat1[ ii + k ] * LowerTriMat2[ kj ];
			}

			l++;
		}
	}

	return TRUE;
}


/***************************************************************************************

  multiplication of two upper triangular matrices of order n

  two lower triangular matrices are vectorized by line, respectively (0s are not included)

     * * * * *
	 0 * * * *
     0 0 * * *
	 0 0 0 * *
	 0 0 0 0 *

  the result matrix is also a lower triangular matrix with the same storage arrangement

***************************************************************************************/

BOOL cMatrix::UpperTriMat_Multiply_UpperTriMat( double *UpperTriMat1, 
		                                        double *UpperTriMat2, 
										        double *UpperTriResult, 
		                                        int n )
{
	if( !UpperTriMat1 || !UpperTriMat2 || !UpperTriResult )
	{
		return FALSE;
	}

	int n1 = n * ( n + 1 ) / 2;
	double *U1 = new double[ n1 ];
	double *U2 = new double[ n1 ];

	if( !U1 || !U2 ) 
	{
		if( U1 ) delete []U1;
		if( U2 ) delete []U2;

		return FALSE;
	}


	Vec_Assign( UpperTriMat1, U1, n1 );
	Vec_Assign( UpperTriMat2, U2, n1 );

	if( UpperTriMat_Transpose( U1, n ) && UpperTriMat_Transpose( U2, n ) &&
        LowerTriMat_Multiply_LowerTriMat( U2, U1, UpperTriResult, n )    &&
		LowerTriMat_Transpose( UpperTriResult, n )                          )
	{
		delete []U1;
		delete []U2;

		return TRUE;
	}

	delete []U1;
	delete []U2;

	return FALSE;
}

/****************************************************************************

  add one matrix with another 

  Mat1, Mat2 is two matrics of dimension  m * n, vectorized by line
  

  Result is a matrix of dimension  m * n, vectorized by line

****************************************************************************/

void cMatrix::Mat_Add_Mat( double *Mat1, double *Mat2, double *Result, int m, int n )
{
	try
	{
		if( !Mat1 || !Mat2 || !Result ) 
		{
			return;
		}

		for ( int k = 0; k < m * n ; k++ )
		{
			Result[ k ] = Mat1[ k ] + Mat2[ k ];
		}	
		return;
	}
	catch( ... )
	{
	}
}

/****************************************************************************

  substract one matrix by another 

  Mat1, Mat2 is two matrices of dimension  m * n, vectorized by line
  

  Result is a matrix of dimension  m * n, vectorized by line

****************************************************************************/
void cMatrix::Mat_Sub_Mat( double *Mat1, double *Mat2, double *Result, int m, int n )
{
	try
	{
		if( !Mat1 || !Mat2 || !Result ) 
		{
			return;
		}

		for ( int k = 0; k < m * n ; k++ )
		{
			Result[ k ] = Mat1[ k ] - Mat2[ k ];
		}	
		return;
	}
	catch( ... )
	{
	}
}

/****************************************************************************

  triangulise a symmetric square matrix, the result is the low triangle part,
  
  Mat is a symmetric square matrix of dimension  n * n, vectorized by line

  Result_TriMat is the low triangle part of the input, vectorized by line

****************************************************************************/

void cMatrix::Mat_Trianglise( double *Mat, double *Result_TriMat, int n )
{
	if( !Mat || !Result_TriMat )
	{
		return;
	}

	try
	{
		int l = 0;

		for( int i = 0; i < n; i++ )
		{
			int k = i * n;
			for( int j = 0; j <= i; j++ )
			{
				Result_TriMat[ l ] = Mat[ k + j ];
				l++;
			}
		}
		return;
	}
	catch( ... )
	{
	}
}

/****************************************************************************

  Recover the full symmetric matrix from the lower triangular part
  
  TriMat: is the low triangle part of the full symmetric matrix of order n, 
          vectorized by line
  Result_Mat is the symmetric square matrix of order n, vectorized by line  

****************************************************************************/

void cMatrix::TriMat2Mat( double *TriMat, double *Result_Mat, int n )
{
	if( !TriMat || !Result_Mat )
	{
		return;
	}

	// first fill the lower triangular part
	int l = 0;
	for( int i = 1; i <= n; i++ )	// row
	{
		for( int j = 1; j <= i; j++ )	// column
		{
			int ij = ( i - 1 ) * n + j - 1;
			Result_Mat[ ij ] = TriMat[ l ];
			l++;
		}
	}

	// then the upper triangular part
	for( int i = 1; i < n; i++ )	// row
	{
		for( int j = i + 1; j <= n; j++ )	// column
		{
			int ij = ( i - 1 ) * n + j - 1;
			int ji = ( j - 1 ) * j / 2 + i - 1;
			Result_Mat[ ij ] = TriMat[ ji ];
		}
	}
}

/****************************************************************************

  The full matrix is a symmetrical and positive definite matrix.

   The matrix is partitioned in the form of

	| m11   m12 |
	| m21   m22	|

   where m11 is a symmetric square matrix of order n1, m22 a symmetric square
   matrix of order n2, m12 is of dimension n1 x n2, m21 = m12(T)

   The inverse of this matrix can be partitioned into 

	| p11   p12 |
	| p21   p22	|

   where p11 is a symmetric square matrix of order n1, p22 a symmetric square
   matrix of order n2, p12 is of dimension n1 x n2, p21 = p12(T)

   Reference: Eq 6.3.35 - 6.3.40, page 394, Tapley et al
   
****************************************************************************/

BOOL cMatrix::PartitionedMatrixInverse( double *p11, double *p12, double *p22, 
									    int n1, int n2 )
{
	if( !p11 || !p12 || !p22 )
	{
		return FALSE;
	}

	int n = __max( n1 * n1, n2 * n2 );
	double *p11o = new double[ n ];
	double *p22o = new double[ n ];
	double *t1 = new double[ n ];
	double *t2 = new double[ n ];
	double *t3 = new double[ n ];
	double *p21 = new double[ n ];

	if( !p11o || !p22o || !t1 || !t2 || !t3 || !p21 )
	{
		if( p11o ) delete []p11o;
		if( p22o ) delete []p22o;

		if( t1 ) delete []t1;
		if( t2 ) delete []t2;
		if( t3 ) delete []t3;
		if( p21 ) delete []p21;

		return FALSE;
	}

	Vec_Assign( p11, p11o, n1 * n1 );
	Vec_Assign( p22, p22o, n2 * n2 );

	cCholeskyDecom CD;

	Mat_Trianglise( p11o, t1, n1 );
	if( !CD.MatrixInversing( t1, n1, 1 ) )
	{
		if( p11o ) delete []p11o;
		if( p22o ) delete []p22o;

		if( t1 ) delete []t1;
		if( t2 ) delete []t2;
		if( t3 ) delete []t3;
		if( p21 ) delete []p21;

		return FALSE;
	}
	TriMat2Mat( t1, p11, n1 );

	Mat_Trianglise( p22o, t2, n2 );
	if( !CD.MatrixInversing( t2, n2, 1 ) )
	{
		if( p11o ) delete []p11o;
		if( p22o ) delete []p22o;

		if( t1 ) delete []t1;
		if( t2 ) delete []t2;
		if( t3 ) delete []t3;
		if( p21 ) delete []p21;

		return FALSE;
	}
	TriMat2Mat( t2, p22, n2 );
	
	Vec_Assign( p12, p21, n1 * n2 );
	Mat_Transpose( p21, n1, n2 );

	// 6.3.38
    Mat_Multiply_Mat( p12, p22, t1, n1, n2, n2 );
	Mat_Multiply_Mat( t1, p21, t2, n1, n2, n1 );
	Mat_Sub_Mat( p11o, t2, t1, n1, n1 );

	Mat_Trianglise( t1, t2, n1 );
	if( !CD.MatrixInversing( t2, n1, 1 ) )
	{
		if( p11o ) delete []p11o;
		if( p22o ) delete []p22o;

		if( t1 ) delete []t1;
		if( t2 ) delete []t2;
		if( t3 ) delete []t3;
		if( p21 ) delete []p21;

		return FALSE;
	}
	TriMat2Mat( t2, t1, n1 );

	// 6.3.39
    Mat_Multiply_Mat( p21, p11, t3, n2, n1, n1 );
	Mat_Multiply_Mat( t3, p12, t2, n2, n1, n2 );
	Mat_Sub_Mat( p22o, t2, t3, n2, n2 );
	Mat_Trianglise( t3, t2, n2 );
	if( !CD.MatrixInversing( t2, n2, 1 ) )
	{
		if( p11o ) delete []p11o;
		if( p22o ) delete []p22o;

		if( t1 ) delete []t1;
		if( t2 ) delete []t2;
		if( t3 ) delete []t3;
		if( p21 ) delete []p21;

		return FALSE;
	}
	TriMat2Mat( t2, t3, n2 );

	// 6.3.37
    Mat_Multiply_Mat( p11, p12, t2, n1, n1, n2 );
	Mat_Multiply_Mat( t2, t3, p12, n1, n2, n2 );
	Vec_Multiply_By_Constant( p12, n1 * n2, -1.0 );

	// assignment
	Vec_Assign( t1, p11, n1 * n1 );
	Vec_Assign( t3, p22, n2 * n2 );

	if( p11o ) delete []p11o;
	if( p22o ) delete []p22o;
	if( t1 ) delete []t1;
	if( t2 ) delete []t2;
	if( t3 ) delete []t3;
	if( p21 ) delete []p21;

	return TRUE;
}

/****************************************************************************

  solve a linear equation system with Gauss-Jordan algorithm

  a is a matrix of dimension n * n, vectorized by line

****************************************************************************/

void cMatrix::GaussianJordan( double *a, double *b, int n )
{
	if( !a || !b )
	{
		return;
	}

	double *A = new double[ n * n ];

	if( !A )
	{
		return;
	}

	try
	{
		Vec_Assign( a, A, n*n );

		for( int i = 0; i < n - 1; i++ )
		{
			int ii = i * n + i;
			double aa = A[ ii ];

			for( int j = i + 1; j < n; j++ )
			{
				int ij = j * n + i;
				double bb = - A[ ij ] / aa;

				for( int k = i + 1; k < n; k++ )
				{
					int ik = i * n + k;
					int jk = j * n + k;
					A[ jk ] += A[ ik ] * bb;
				}

				b[ j ] += b[ i ] * bb;
			}
		}

		b[ n - 1 ] /= A[ n * n - 1 ];

		for( int k = n - 2; k > -1; k-- )
		{
			double cc = b[ k ];
			for( int j = n - 1; j > k; j-- )
			{
				int jk = k * n + j;
				cc -= A[ jk ] * b[ j ]; 
			}

			b[ k ] = cc / A[ k * n + k ];
		}

		delete []A;

		return;
	}
	catch( ... )
	{
	}
}


/****************************************************************************

  To get the vector product of two vectors of 3 elements

****************************************************************************/

void cMatrix::VectorProduct( double *pdfVec1, double *pdfVec2, double *pdfVec3 )
{
	if( !pdfVec1 || !pdfVec2 || !pdfVec3 )
	{
		return;
	}

	try 
	{
		pdfVec3[ 0 ] = pdfVec1[ 1 ] * pdfVec2[ 2 ] - pdfVec1[ 2 ] * pdfVec2[ 1 ];
		pdfVec3[ 1 ] = pdfVec1[ 2 ] * pdfVec2[ 0 ] - pdfVec1[ 0 ] * pdfVec2[ 2 ];
		pdfVec3[ 2 ] = pdfVec1[ 0 ] * pdfVec2[ 1 ] - pdfVec1[ 1 ] * pdfVec2[ 0 ];
		return;
	}
	catch( ... )
	{
	}
}

/****************************************************************************

  to determine eigenvalues and corresponding eigenvectors of a real, 
  symmetric matrix.

  This method is a translation of the Fortran Subroutines of tred2, 
  on pp 467-468, and tqli, on pp. 473-474, of Numerical Recipes in Fortran,
  by Press et al.

  TriMat: the low triangle of a symmetric matrix of dimension n, 
          vectorized by line (0-based)

  evalues: determined eigenvalues (0-based)
  evectors: determined eigenvectors vectoried by line (0-based), with the 
            following meaning 

				e11  e21  e31  e41 .... en1
				e12  e22  e32  e42 .... en2
				e13  e23  e33  e43 .... en3
				....
				e1n  e2n  e3n  e4n .... enn
            
			where e11, e12, e13 ... e1n is the normalized eigenvector corresponding to
			eigenvalue evalues[ 0 ], e21, e22, e23 ... e2n is the normalized eigenvector 
			corresponding to eigenvalue evalues[ 1 ], etc

  bValueOnly: TRUE if only eigenvalues are required.

****************************************************************************/

BOOL cMatrix::EigensSymmMatrix( double *TriMat, int n, 
							    double *evalues, double *evectors, 
								BOOL bValueOnly )

{
	if( !TriMat || !evalues || !evectors )
	{
		return FALSE;
	}

	double *a = new double[ n * n + 1 ];
	double *d = new double[ n + 1 ];
	double *e = new double[ n + 1 ];

	if( !a || !d || !e )
	{
		if( a ) delete []a;
		if( d ) delete []d;
		if( e ) delete []e;

		return FALSE;
	}

	int i, j, k, l;
	double f, g, h, hh, scale;

	// expand the given TriMat into a full square matrix which is 1 based index

	l = 1;
	for( i = 0; i < n; i++ )
	{
		k = ( i + 1 ) * ( i + 2 ) / 2 - 1 - i;
		for( j = 0; j <= i; j++ )
		{
			a[ l ] = TriMat[ k + j ];
			l++;
		}

		for( j = i + 1; j < n; j++ )
		{
			k = ( j + 1 ) * ( j + 2 ) / 2 - j + i - 1;
			a[ l ] = TriMat[ k ];
			l++;
		}
	}

	// copy of tred2
	for( i = n; i >= 2; i--)
	{
		l = i - 1;
		h = 0.0;
		scale = 0.0;

		if( l > 1 )
		{
			for( k = 1; k <= l; k++ ) scale += fabs( a[ i * n - n + k ] );

			if( scale == 0.0 ) // skip transformation
				e[ i ] = a[ i * n - n + l ];
			else
			{
				for( k = 1; k <= l; k++ )
				{
					a[ i * n - n + k ] /= scale;
					h += a[ i * n - n + k ] * a[ i * n - n + k ];	// form sigma in h
				}
				f = a[ i * n - n + l ];
				g = ( f >= 0.0 ? -sqrt( h ) : sqrt( h ) );
				e[ i ] = scale * g;
				h -= f * g;
				a[ i * n - n + l ] = f - g;
				f = 0.0;
				for( j = 1; j <= l; j++ )
				{
					a[ j * n - n + i ] = a[ i * n - n + j ] / h;
					g = 0.0;
					for( k = 1; k <= j; k++ ) g += a[ j * n - n + k ] * a[ i * n - n + k ];
					for( k = j + 1; k <= l; k++ ) g += a[ k * n - n + j ] * a[ i * n - n + k ];
					e[ j ] = g / h;
					f += e[ j ] * a[ i * n - n + j ];
				}

				hh = f / ( h + h );
				for( j = 1; j <= l; j++ ) 
				{
					f = a[ i * n - n + j ];
					g = e[ j ] - hh * f;
					e[ j ] = g;
					for( k = 1; k <= j; k++ ) a[ j * n - n + k ] -= f * e[ k ] + g * a[ i * n - n + k ];
				}
			}
		}
		else
			e[ i ] = a[ i * n - n + l ];

		d[ i ] = h;
	}

	d[ 1 ] = e[ 1 ] = 0.0;

	if( bValueOnly )
	{
		for( i = 1; i <= n; i++ ) d[ i ] = a[ i * n - n + i ];
	}
	else
	{
		for( i = 1; i <= n; i++ ) 
		{
			l = i - 1;
			if( d[ i ] != 0.0 ) 
			{
				for( j = 1; j <= l; j++ )
				{
					g = 0.0;
					for( k = 1; k <= l; k++ ) g += a[ i * n - n + k ] * a[ k * n - n + j ];
					for( k = 1; k <= l; k++ ) a[ k * n - n + j ] -= g * a[ k * n - n + i ];
				}
			}
			d[ i ] = a[ i * n - n + i ];
			a[ i * n - n + i ] = 1.0;
			for( j = 1; j <= l; j++ ) a[ j * n - n + i ] = a[ i * n - n + j ] = 0.0;
		}
	}

	// copy of tqli
	int iter, m;
	double b, c, dd, p, r, s;

	for( i = 2; i <= n; i++ ) e[ i - 1 ] = e[ i ];
	e[ n ] = 0.0;

	for( l = 1; l <= n; l++ )
	{
		iter = 0;
ONE:
		for( m = l; m <= n - 1; m++ )
		{
			dd = fabs( d[ m ] ) + fabs( d[ m + 1 ] );
			if( fabs( e[ m ] ) + dd == dd ) goto TWO;
		}
		m = n;

TWO:
		if( m != l )
		{
			if( iter == 30 )
			{
				delete []a;
				delete []d;
				delete []e;

				return FALSE;
			}

			iter++;

			g = ( d[ l + 1 ] - d[ l ] ) / ( 2.0 * e[ l ] );
			r = pythag( g, 1.0 );
			g = d[ m ] - d[ l ] + e[ l ] / ( g + ( g >= 0.0 ? r : -r ) );
			s = c = 1.0;
			p = 0.0;
			for( i = m - 1; i >= l; i-- )
			{
				f = s * e[ i ];
				b = c * e[ i ];
				r = pythag( f, g );
				e[ i + 1 ] = r;
				if( r == 0.0 ) 
				{
					d[ i + 1 ] -= p;
					e[ m ] = 0.0;
					goto ONE;
				}
				
				s = f / r;
				c = g / r;
				g = d[ i + 1 ] - p;
				r = ( d[ i ] - g ) * s + 2.0 * c * b;
				p = s * r;
				d[ i + 1 ] = g + p;
				g = c * r - b;

				if( !bValueOnly )
				{
					for( k = 1; k <= n; k++ )
					{
						f = a[ k * n - n + i + 1 ];
						a[ k * n - n + i + 1 ] = s * a[ k * n - n + i ] + c * f;
						a[ k * n - n + i ] = c * a[ k * n - n + i ] - s * f;
					}
				}
			}
				
			d[ l ] -= p;
			e[ l ] = g;
			e[ m ] = 0.0;
			goto ONE;
		}
	}

	for( i = 0; i < n; i++ ) evalues[ i ] = d[ i + 1 ];

	if( !bValueOnly ) for( i = 0; i < n * n; i++ ) evectors[ i ] = a[ i + 1 ];

	delete []a;
	delete []d;
	delete []e;

	return TRUE;
}

/******************************************************************************

  compute sqrt( a * a + b * b ) without destructive underflow or overflow

  subroutine pythag on pp. 62-63 of Numerical Recipe of Press et al.

******************************************************************************/
double cMatrix::pythag( double a, double b )
{
	double r;

	double absa = fabs( a );
	double absb = fabs( b );

	if( absa > absb ) 
	{
		r = absb / absa;
		r = r * r;
		r = absa * sqrt( 1.0 + r );
	}
	else
	{
		if( absb == 0.0 ) r = 0;
		else
		{
			r = absa / absb;
			r = r * r;
			r = absb * sqrt( 1.0 + r );
		}
	}

	return r;
}

//compute vector norm2
double cMatrix::Vec_norm2( double *Vec, int n)
{
	double sum = 0;

	for( int i = 0; i < n; ++i )
	{
		sum += pow( Vec[ i ], 2.0 );
	}

	double normA = sqrt( sum );

	return normA;
}


double cMatrix::Vec_Norm(double* Vec1, int n)
{
	if (!Vec1)
	{
		return 1.0e20;
	}

	double res = 0;

	for (int i = 0; i < n; i++) res = res + Vec1[i] * Vec1[i];

	res = sqrt(res);
	return res;
}



void cMatrix::Vec_Multiply_By_Constant(double *Vec, int n, double c, double *ResultVec)
{
	if (!Vec)
	{
		return;
	}

	for (int i = 0; i < n; i++) ResultVec[i] = c * Vec[i];
}


void cMatrix::Vec_Minus_Vec(double *Vec1, double *Vec2, int n, double *ResultVec)
{
	if (!Vec1 || !Vec2)
	{
		return;
	}

	for (int i = 0; i < n; i++) ResultVec[i] = Vec1[i] - Vec2[i];
}
