/* Author:  G. Jungman
 * RCS:     $Id$
 */

#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"

extern int gsl_sf_bessel_Y0_impl(const double, double *);
extern int gsl_sf_bessel_Y1_impl(const double, double *);


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Yn_impl(int n, const double x, double * result)
{
  int sign = 1;
  const double xmax = n * 500.;

  if(n < 0) {
    /* reduce to case n >= 0 */
    n = -n;
    if(GSL_IS_ODD(n)) sign = -1;
  }
  
  if(n == 0) {
    double b0;
    int status = gsl_sf_bessel_Y0_impl(x, &b0);
    *result = sign * b0;
    return status;
  }
  else if(n == 1) {
    double b0;
    int status = gsl_sf_bessel_Y1_impl(x, &b0);
    *result = sign * b0;
    return status;
  }
  else {
    if(x <= 0.) {
      return GSL_EDOM;
    }
    if(x < 1.999999999) {
      double den = M_PI * gsl_sf_pow_int(0.5*x, n);
      if(den == 0.) {
	return GSL_EDOM;
      }
    }
    else if(x < xmax){
      int j;
      double by, bym, byp;
      double two_over_x = 2.0/x;
      gsl_sf_bessel_Y1_impl(x, &by);
      gsl_sf_bessel_Y0_impl(x, &bym);
      for(j=1; j<n; j++) { 
	byp = j*two_over_x*by - bym;
	bym = by;
	by  = byp;
      }
      *result = sign * by;
      return GSL_SUCCESS;/* GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW */
    }
    else {
      double ampl  = gsl_sf_bessel_asymp_Mnu(n, x);
      double theta = gsl_sf_bessel_asymp_thetanu(n, x);
      *result = sign * ampl * sin(theta);
      return GSL_SUCCESS;
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Yn_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_Yn_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Yn_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_Yn(const int n, const double x)
{
  double y;
  int status = gsl_sf_bessel_Yn_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Yn", status);
  }
  return y;
}
