
/*****************************************************************************

         Routine Matrix computations
         --------------------------- 

******************************************************************************/

#if ! defined( _INC_MATRIXCOM )
#define _INC_MATRIXCOM

#include <windows.h>

extern class cMatrix g_Matrix;

class cMatrix
{

public:
    cMatrix();
    ~cMatrix();

    void Vec_Assign( double *Source, double *Result, int n );
    void Set_Vec_Zero( double *Vec, int n );

    void Vec_Multiply_By_Constant( double *Vec, int n, double c );
    void TriMat_Multiply_By_Constant( double *TriMat, int n,double c );
    void Mat_Multiply_By_Constant( double *Mat, int m, int n, double c );

    BOOL Mat_Transpose( double *Mat, int m, int n );
	BOOL UpperTriMat_Transpose( double *TriMat, int n );
	BOOL LowerTriMat_Transpose( double *TriMat, int n );

    BOOL TriMat_Inverse( double *TriMat, int n );
	BOOL Mat_Inverse( double *Mat, int n );
     
    double Vec_Multiply_Vec( double *Vec1, double *Vec2, int n );
    void Vec_Multiply_Vec_To_Mat( double *Vec1, double *Vec2, double *Mat, 
		                          int n1, int n2 );
    void Vec_Multiply_Vec_To_TriMat( double *Vec, double *TriMat, int n );
    void Vec_Multiply_TriMat( double *Vec, double *TriMat, double *Result, 
		                      int n );
    void Vec_Multiply_Mat( double *Vec, double *Mat, double *Result, 
		                   int n, int m );

    void TriMat_Multiply_Vec( double *TriMat, double *Vec, double *Result, 
		                      int n );
	void Mat_Multiply_Vec( double *Mat, double *Vec, double *Result, 
		                   int m, int n );

    void TriMat_Multiply_TriMat( double *TriMat1, double *TriMat2, 
		                         double *Result, int n );
    void TriMat_Multiply_Mat( double *TriMat, double *Mat, double *Result, 
		                      int n, int m );
    void Mat_Multiply_TriMat( double *Mat, double *TriMat, double *Result, 
		                      int m, int n );
    void Mat_Multiply_Mat( double *Mat1, double *Mat2, double *Result, 
		                   int m, int n, int p );
    BOOL LowerTriMat_Multiply_LowerTriMat( double *LowerTriMat1, 
		                                   double *LowerTriMat2, 
										   double *LowerTriResult, 
		                                   int n );
    BOOL UpperTriMat_Multiply_UpperTriMat( double *UpperTriMat1, 
		                                   double *UpperTriMat2, 
										   double *UpperTriResult, 
		                                   int n );
	void Mat_Add_Mat( double *Mat1, double *Mat2, double *Result, int m, int n );
	void Mat_Sub_Mat( double *Mat1, double *Mat2, double *Result, int m, int n );

    void Mat_Trianglise( double *Mat, double *Result_TriMat, int n );
    void TriMat2Mat( double *TriMat, double *Result_Mat, int n );

	BOOL PartitionedMatrixInverse( double *p11, double *p12, double *p22, int n1, int n2 );

	void GaussianJordan( double *a, double *b, int n );

	void VectorProduct( double *pdfVec1, double *pdfVec2, double *pdfVec3 );

	BOOL EigensSymmMatrix( double *TriMat, int n, double *evalues, double *evectors, BOOL bValueOnly );

	double pythag( double a, double b );
	double Vec_norm2( double *Vec, int n);

	double Vec_Norm(double* Vec1, int n);
	void Vec_Multiply_By_Constant(double *Vec, int n, double c, double *ResultVec);
	void Vec_Minus_Vec(double *Vec1, double *Vec2, int n, double *ResultVec);

};

#endif