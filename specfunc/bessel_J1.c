/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#include "bessel_amp_phase.h"

#define ROOT_EIGHT 2.82842712474619009760337744842


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besj1, 1983 version, w. fullerton */

/* chebyshev expansions

 series for bj1        on the interval  0.	    to  1.60000d+01
					with weighted error   4.48e-17
					 log weighted error  16.35
			       significant figures required  15.77
				    decimal places required  16.89

*/
static double bj1_data[12] = {
  -.11726141513332787,
  -.25361521830790640,
   .050127080984469569,
  -.004631514809625081,
   .000247996229415914,
  -.000008678948686278,
   .000000214293917143,
  -.000000003936093079,
   .000000000055911823,
  -.000000000000632761,
   .000000000000005840,
  -.000000000000000044,
};
static struct gsl_sf_ChebSeries bj1_cs = {
  bj1_data,
  11,
  -1, 1,
  (double *)0,
  (double *)0
};

int gsl_sf_bessel_J1_impl(const double x, double * result)
{ 
  const double xmin = ROOT_EIGHT * GSL_SQRT_MACH_EPS;
  const double xmax = 1./GSL_SQRT_MACH_EPS;
  double y = fabs(x);
  
  if(y == 0.) {
    *result = 0.;
    return GSL_SUCCESS;
  }
  else if(y < 2.*DBL_MIN) {
    *result = 0.;
    return GSL_EUNDRFLW;
  }
  else if(y < xmin) {
    *result = 0.5*x;
    return GSL_SUCCESS;
  }
  else if(y < 4.0) {
    *result = x * (.25 + gsl_sf_cheb_eval(.125*y*y-1., &bj1_cs));
    return GSL_SUCCESS;
  }
  else if(y < xmax) {
    double z     = 32.0/(y*y) - 1.0;
    double ampl  = (0.75 + gsl_sf_cheb_eval(z, &_bessel_amp_phase_bm1_cs)) / sqrt(y);
    double theta = y - 3.0*M_PI_4 + gsl_sf_cheb_eval(z, &_bessel_amp_phase_bth1_cs) / y;
    *result = (x < 0. ? -ampl : ampl) /* fortran_sign(ampl, x) */ * cos (theta);
    return GSL_SUCCESS;
  }
  else {
    double ampl  = gsl_sf_bessel_asymp_Mnu(1., y);
    double theta = gsl_sf_bessel_asymp_thetanu(1., y);
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
