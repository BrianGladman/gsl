#ifndef _GSL_MATH_H_
#include <math.h>
#include <limits.h>
#include <float.h>

#ifndef M_SQRT2
#define M_SQRT2		1.41421356237309504880
#endif

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#ifndef M_PI_4
#define M_PI_4        	0.78539816339744830962      /* pi/4 */
#endif

#ifndef M_LN10
#define M_LN10     	      	2.30258509299404568402	    /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2 	      	      	0.69314718055994530942	    /* ln(2) */
#endif

#ifndef M_EULER
#define M_EULER                 0.5772156649015328606065    /* Euler constant */
#endif


/* need to determine this stuff at configure time, just use placeholder
   values for now */

#define GSL_MACH_EPS		1.e-14
#define GSL_SQRT_MACH_EPS	1.e-7
#define GSL_ROOT3_MACH_EPS      2.154e-5
#define GSL_ROOT4_MACH_EPS      0.0003162

#define GSL_SQRT_DBL_MIN        2.e-154
#define GSL_LOG_DBL_MIN       	(DBL_MIN_10_EXP * M_LN10)
#define GSL_LOG_DBL_MAX       	(DBL_MAX_10_EXP * M_LN10)


/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)


#endif /* !_GSL_MATH_H_ */
