/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel_amp_phase.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besy1, 1977 version, w. fullerton */

/* chebyshev expansions

 series for by1        on the interval  0.	    to  1.60000d+01
					with weighted error   1.87e-18
					 log weighted error  17.73
			       significant figures required  17.83
				    decimal places required  18.30
*/

static double by1_data[14] = {
  0.03208047100611908629,
  1.262707897433500450,
  0.00649996189992317500,
 -0.08936164528860504117,
  0.01325088122175709545,
 -0.00089790591196483523,
  0.00003647361487958306,
 -0.00000100137438166600,
  0.00000001994539657390,
 -0.00000000030230656018,
  0.00000000000360987815,
 -0.00000000000003487488,
  0.00000000000000027838,
 -0.00000000000000000186
};
static gsl_sf_cheb_series by1_cs = {
  by1_data,
  13,
  -1, 1,
  (double *)0,
  (double *)0,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Y1_impl(const double x, double * result)
{
  const double two_over_pi = 2.0/M_PI;
  const double xmin = 1.571*DBL_MIN; /*exp ( amax1(alog(r1mach(1)), -alog(r1mach(2)))+.01) */
  const double x_small = 2.0 * GSL_SQRT_MACH_EPS;
  const double xmax    = 1.0/GSL_MACH_EPS;
  
  if(x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < xmin) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(x < x_small) {
    double J1 = 0.0;
    int status = gsl_sf_bessel_J1_impl(x, &J1);
    *result = two_over_pi * log(0.5*x) * J1 + (0.5 + gsl_sf_cheb_eval(&by1_cs, -1.0))/x;
    return status;
  }
  else if(x < 4.0) {
    double J1 = 0.0;
    int status = gsl_sf_bessel_J1_impl(x, &J1);
    *result = two_over_pi * log(0.5*x) * J1 + (0.5 + gsl_sf_cheb_eval(&by1_cs, 0.125*x*x-1.0))/x;
    if(status == GSL_EUNDRFLW)
      return GSL_ELOSS;
    else
      return status;
  }
  else if(x < xmax) {
    double z     = 32.0/(x*x) - 1.0;
    double ampl  = (0.75 + gsl_sf_cheb_eval(&_bessel_amp_phase_bm1_cs, z)) / sqrt(x);
    double theta = x - 3.0*M_PI_4 + gsl_sf_cheb_eval(&_bessel_amp_phase_bth1_cs, z) / x;
    *result = ampl * sin (theta);
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Y1_e(const double x, double * result)
{
  int status = gsl_sf_bessel_Y1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Y1_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_Y1(const double x)
{
  double y;
  int status = gsl_sf_bessel_Y1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Y1", status);
  }
  return y;
}
