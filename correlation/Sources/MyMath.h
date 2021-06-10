#ifndef USED_MYMATH
#define USED_MYMATH
#include "Math.h"

__forceinline int round(float x) 
{
	#ifdef WIN32
		int retval;
		__asm fld x
		__asm fistp retval
		return retval;
	#else
		return (((x - floor(x)) < 0.5) ? (floorf(x)) : (ceilf(x)));
	#endif
}

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#endif