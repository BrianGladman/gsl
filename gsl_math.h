#ifndef _GSL_MATH_H
#define _GSL_MATH_H
#include <math.h>

#ifndef M_E
#define M_E        2.71828182845904523536028747135	/* e */
#endif

#ifndef M_SQRT2
#define M_SQRT2	   1.41421356237309504880168872421
#endif

#ifndef M_SQRT3
#define M_SQRT3	   1.73205080756887729352744634151
#endif

#ifndef M_PI
#define M_PI	   3.14159265358979323846264338328
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830966156608458	/* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334	/* sqrt(pi) */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468	/* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2 	   0.69314718055994530941723212146	/* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135	/* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008	/* Euler constant */
#endif


/* magic constants; mostly for the benefit of the implementation */
#include <gsl_machine.h>
#include <gsl_precision.h>


/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)


/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

extern inline int GSL_MAX_INT (int a, int b);
extern inline int GSL_MIN_INT (int a, int b);
extern inline double GSL_MAX_DBL (double a, double b);
extern inline double GSL_MIN_DBL (double a, double b);
extern inline long double GSL_MAX_LDBL (long double a, long double b);
extern inline long double GSL_MIN_LDBL (long double a, long double b);

extern inline int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

extern inline int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

extern inline double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

extern inline double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

extern inline long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

extern inline long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */


#endif /* !_GSL_MATH_H */
