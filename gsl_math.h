#ifndef _GSL_MATH_H_
#include <math.h>
#include <limits.h>
#include <float.h>

#ifndef M_E
#define M_E             2.71828182845904523536028747135 /* e */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.41421356237309504880168872421
#endif

#ifndef M_SQRT3
#define M_SQRT3		1.73205080756887729352744634151
#endif

#ifndef M_PI
#define M_PI		3.14159265358979323846264338328
#endif

#ifndef M_PI_4
#define M_PI_4        	0.785398163397448309661566084582    /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI        1.77245385090551602729816748334     /* sqrt(pi) */
#endif

#ifndef M_LN10
#define M_LN10     	2.30258509299404568401799145468	    /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2 	      	0.69314718055994530941723212146	    /* ln(2) */
#endif

#ifndef M_EULER
#define M_EULER         0.57721566490153286060651209008     /* Euler constant */
#endif


/* need to determine this stuff at configure time, just use placeholder
   values for now */

#define GSL_MACH_EPS		1.e-14
#define GSL_SQRT_MACH_EPS	1.e-7
#define GSL_ROOT3_MACH_EPS      2.154e-5
#define GSL_ROOT4_MACH_EPS      0.0003162
#define GSL_ROOT5_MACH_EPS      0.00158489
#define GSL_ROOT6_MACH_EPS      0.00464159
#define GSL_LOG_MACH_EPS      	-32.2362

#define GSL_SQRT_DBL_MIN        2.e-154
#define GSL_SQRT_DBL_MAX        5.e+153
#define GSL_ROOT3_DBL_MIN       3.42e-103
#define GSL_ROOT3_DBL_MAX       2.92e+102
#define GSL_ROOT4_DBL_MIN       1.414213e-77
#define GSL_ROOT4_DBL_MAX       7.071067e+76
#define GSL_LOG_DBL_MIN       	(DBL_MIN_10_EXP * M_LN10)
#define GSL_LOG_DBL_MAX       	(DBL_MAX_10_EXP * M_LN10)


/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)


#endif /* !_GSL_MATH_H_ */
