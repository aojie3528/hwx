/***************************************************************************************

 compute astronomical arguments

 Reference: IERS TECH NOTE 21

***************************************************************************************/

#if!defined INC_ASTROARGUMENT
#define INC_ASTROARGUMENT

class cAstroArg
{
public:

	cAstroArg();
	~cAstroArg();

	static void ComputeFundamentals( double dfT, double &dfL, double &dfLp,
	                                 double &dfF, double &dfD, double &dfOm );

};

#endif