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


int gsl_sf_pow_int_impl(double x, int n, gsl_sf_result * result)
{
  double value = 1.0;
  int count = 0;

  if(result == 0) {
    return GSL_EFAULT;
  }

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
     ++count;
  } while (n);

  result->val = value;
  result->err = 2.0 * GSL_DBL_EPSILON * (count + 1.0) * fabs(value); 

  return GSL_SUCCESS;
}


int gsl_sf_pow_int_e(const double x, const int n, gsl_sf_result * result)
{
  int status = gsl_sf_pow_int_impl(x, n, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_pow_int", status);
  }
  return status;
  
}

double gsl_sf_pow_int(const double x, const int n)
{
  gsl_sf_result p;
  int status = gsl_sf_pow_int_impl(x, n, &p);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_pow_int", status);
  }
  return p.val;
}
