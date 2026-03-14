/***************************************************************************************

 Lagrangian Interpolation

  Jizhang Sang, EOS

 (C) COPYRIGHT -Electro Optic Systems Pty. Ltd. Australia (EOS) 2002,
 All rights reserved. No part of this program may be photocopied, 
 reproduced or translated to another program language without the prior
 written consent of EOS

 Reference: Eq (3.1.1), P. 102, Numerical Recipes in Fortran, William Press, et al.

***************************************************************************************/
// #include "stdafx.h"

#include "LagrangianInterpolation.h"
#include <vector>


using namespace std;

cLagInterpolation::cLagInterpolation()
{

}

cLagInterpolation::~cLagInterpolation()
{

}

/***************************************************************************************

  Fortran subroutine: Lagrangian_interpolation

  Lagrange interpolation by given a number of time/obs.

***************************************************************************************/

bool cLagInterpolation::Interpolate( int nNumber, double *pdfTime, double *pdfObs, 
									 double dfTime, double &dfObs )
{
	double dfTemp;

	vector <double> v_time;

	int i, j;

	for( i = 0; i < nNumber; i++ )
	{
		dfTemp = dfTime - pdfTime[ i ];
		v_time.emplace_back(dfTemp);
	}

	dfObs = 0.0;

	for( i = 0; i < nNumber; i++ )
	{
		dfTemp = pdfObs[ i ];

		for( j = 0; j < nNumber; j++ )
		{
			if( i != j )
				dfTemp = dfTemp * v_time[ j ] / ( pdfTime[ i ] - pdfTime[ j ] );
		}

		dfObs += dfTemp;
	}

	return true;
}
