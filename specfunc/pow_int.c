/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_pow_2(const double x) { return x*x;   }
double gsl_sf_pow_3(const double x) { return x*x*x; }
double gsl_sf_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
double gsl_sf_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
double gsl_sf_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
double gsl_sf_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
double gsl_sf_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
double gsl_sf_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }


double gsl_sf_pow_int(double x, int n)
{
  double value = 1.0;

  if(n < 0) {
    if(x == 0.0) return 0.0; /* FIXME: should be Inf */
    x = 1.0/x;
    n = -n;
  }

  /* repeated squaring method 
   * returns 0.0^0 = 1.0, so continuous in x
   */
  do {
     if(GSL_IS_ODD(n)) value *= x;
     n >>= 1;
     x *= x;
  } while (n);

  return value; 
}
