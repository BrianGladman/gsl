/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel_amp_phase.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#define ROOT_EIGHT (2.0*M_SQRT2)


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* based on SLATEC besj1, 1983 version, w. fullerton */

/* chebyshev expansions

 series for bj1        on the interval  0.	    to  1.60000d+01
					with weighted error   4.48e-17
					 log weighted error  16.35
			       significant figures required  15.77
				    decimal places required  16.89

*/
static double bj1_data[12] = {
  -0.11726141513332787,
  -0.25361521830790640,
   0.050127080984469569,
  -0.004631514809625081,
   0.000247996229415914,
  -0.000008678948686278,
   0.000000214293917143,
  -0.000000003936093079,
   0.000000000055911823,
  -0.000000000000632761,
   0.000000000000005840,
  -0.000000000000000044,
};
static gsl_sf_cheb_series bj1_cs = {
  bj1_data,
  11,
  -1, 1,
  (double *)0,
  (double *)0,
  8
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_J1_impl(const double x, double * result)
{ 
  const double xmin = ROOT_EIGHT * GSL_SQRT_MACH_EPS;
  const double xmax = 1.0/GSL_SQRT_MACH_EPS;
  double y = fabs(x);
  
  if(y == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(y < 2.*DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(y < xmin) {
    *result = 0.5*x;
    return GSL_SUCCESS;
  }
  else if(y < 4.0) {
    *result = x * (0.25 + gsl_sf_cheb_eval(&bj1_cs, 0.125*y*y-1.0));
    return GSL_SUCCESS;
  }
  else if(y < xmax) {
    double z     = 32.0/(y*y) - 1.0;
    double ampl  = (0.75 + gsl_sf_cheb_eval(&_bessel_amp_phase_bm1_cs, z)) / sqrt(y);
    double theta = y - 3.0*M_PI_4 + gsl_sf_cheb_eval(&_bessel_amp_phase_bth1_cs, z) / y;
    *result = (x < 0.0 ? -ampl : ampl) /* fortran_sign(ampl, x) */ * cos (theta);
    return GSL_SUCCESS;
  }
  else {
    double ampl  = gsl_sf_bessel_asymp_Mnu(1.0, y);
    double theta = gsl_sf_bessel_asymp_thetanu(1.0, y);
    *result = ampl * cos(theta);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_J1_e(const double x, double * result)
{
  int status = gsl_sf_bessel_J1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_J1_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_J1(const double x)
{
  double y;
  int status = gsl_sf_bessel_J1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_J1", status);
  }
  return y;
}
