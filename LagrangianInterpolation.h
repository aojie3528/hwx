/***************************************************************************************
 
 Lagrangian Interpolation

  Jizhang Sang, EOS

 (C) COPYRIGHT -Electro Optic Systems Pty. Ltd. Australia (EOS) 2002,
 All rights reserved. No part of this program may be photocopied, 
 reproduced or translated to another program language without the prior
 written consent of EOS

 Reference: Eq (3.1.1), P. 102, Numerical Recipes in Fortran, William Press, et al.

***************************************************************************************/

#ifndef INC_LAGINTERPOLATION
#define INC_LAGINTERPOLATION



class cLagInterpolation
{

public:

	cLagInterpolation();
	~cLagInterpolation();

	bool Interpolate( int nNumber, double *pdfTime, double *pdfObs, 
		              double dfTime, double &dfObs );

};

#endif
