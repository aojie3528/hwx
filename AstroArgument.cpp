/***************************************************************************************

 compute astronomical arguments

***************************************************************************************/

#include "AstroArgument.h"
#include "Constant.h"

cAstroArg::cAstroArg()
{

}

cAstroArg::~cAstroArg()
{

}

/***************************************************************************************

 Compute 5 fundamental elements at given time

 dfT is the centuries of epoch since J2000
 see p. 23, IERS Tech Note 21

***************************************************************************************/

void cAstroArg::ComputeFundamentals( double dfT, double &dfL, double &dfLp,
                                     double &dfF, double &dfD, double &dfOm )
{
	double dfT2 = dfT * dfT,
		   dfT3 = dfT2 * dfT,
		   dfT4 = dfT3 * dfT;

	// computation of fundamental arguments

	// mean longitude of Moon
	dfL =   485868.249036 + 1717915923.2178 * dfT +
		    31.8792 * dfT2 + 0.051635 * dfT3 -
			0.00024470 * dfT4;	// in arc seconds

	// mean longitude of Sun
	dfLp=   1287104.793048 + 129596581.0481 * dfT -
		    0.5532 * dfT2 - 0.000136 * dfT3 -
			0.00001149 * dfT4;	// in arc seconds


	dfF =   335779.526232 + 1739527262.8478 * dfT -
		    12.7512 * dfT2 - 0.001037 * dfT3 +
			0.00000417 * dfT4;	// in arc seconds

	dfD =   1072260.703692  + 1602961601.2090 * dfT -
		    6.3706 * dfT2 + 0.006593 * dfT3 -
			0.00003169 * dfT4;	// in arc seconds

	dfOm=   450160.398036 - 6962890.2665 * dfT +
		    7.4722 * dfT2 + 0.007702 * dfT3 -
			0.00005939 * dfT4;	// in arc seconds

	// in radians
	dfL *= g_dfSEC2RAD;
	dfLp *= g_dfSEC2RAD;
	dfF *= g_dfSEC2RAD;
	dfD *= g_dfSEC2RAD;
	dfOm *= g_dfSEC2RAD;
}
