/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_expint.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC shi.f, W. Fullerton

 series for shi  on the interval  0.00000e+00 to  1.40625e-01
					with weighted error   4.67e-20
					 log weighted error  19.33
			       significant figures required  17.07
				    decimal places required  19.75
*/
static double shi_data[7] = {
   0.0078372685688900950695,
   0.0039227664934234563973,
   0.0000041346787887617267,
   0.0000000024707480372883,
   0.0000000000009379295591,
   0.0000000000000002451817,
   0.0000000000000000000467
};
static gsl_sf_cheb_series shi_cs = {
  shi_data,
  6,
  -1, 1,
  (double *)0,
  (double *)0,
  6
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_Shi_impl(const double x, double * result)
{
  const double xsml = GSL_SQRT_MACH_EPS;  /* sqrt (r1mach(3)) */
  const double ax   = fabs(x);

  if(ax < xsml) {
    *result = x;
    return GSL_SUCCESS;
  }
  else if(ax <= 0.375) {
    *result = x * (1.0 + gsl_sf_cheb_eval(&shi_cs, 128.0*x*x/9.0-1.0));
    return GSL_SUCCESS;
  }
  else {
    double Ei, E1;
    int status_Ei = gsl_sf_expint_Ei_impl(x, &Ei);
    int status_E1 = gsl_sf_expint_E1_impl(x, &E1);
    *result = 0.5*(Ei + E1);
    if(status_Ei == GSL_EUNDRFLW && status_E1 == GSL_EUNDRFLW) {
      return GSL_EUNDRFLW;
    }
    else if(status_Ei == GSL_EOVRFLW || status_E1 == GSL_EOVRFLW) {
      return GSL_EOVRFLW;
    }
    else {
      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_Chi_impl(const double x, double * result)
{
  double Ei, E1;
  int status_Ei = gsl_sf_expint_Ei_impl(x, &Ei);
  int status_E1 = gsl_sf_expint_E1_impl(x, &E1);
  if(status_Ei == GSL_EDOM || status_E1 == GSL_EDOM) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(status_Ei == GSL_EUNDRFLW && status_E1 == GSL_EUNDRFLW) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(status_Ei == GSL_EOVRFLW || status_E1 == GSL_EOVRFLW) {
    *result = 0.0;
    return GSL_EOVRFLW;
  }
  else {
    *result = 0.5 * (Ei - E1);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_Shi_e(const double x, double * result)
{
  int status = gsl_sf_Shi_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_Shi_e", status);
  }
  return status;
}

int gsl_sf_Chi_e(const double x, double * result)
{
  int status = gsl_sf_Chi_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_Chi_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_Shi(const double x)
{
  double y;
  int status = gsl_sf_Shi_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_Shi", status);
  }
  return y;
}

double gsl_sf_Chi(const double x)
{
  double y;
  int status = gsl_sf_Chi_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_Chi", status);
  }
  return y;
}
