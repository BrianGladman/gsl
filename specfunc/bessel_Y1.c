/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#include "bessel_amp_phase.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besy1, 1977 version, w. fullerton */

/* chebyshev expansions

 series for by1        on the interval  0.	    to  1.60000d+01
					with weighted error   1.87e-18
					 log weighted error  17.73
			       significant figures required  17.83
				    decimal places required  18.30
*/

static double by1_data[14] = {
   .03208047100611908629,
  1.262707897433500450,
   .00649996189992317500,
  -.08936164528860504117,
   .01325088122175709545,
  -.00089790591196483523,
   .00003647361487958306,
  -.00000100137438166600,
   .00000001994539657390,
  -.00000000030230656018,
   .00000000000360987815,
  -.00000000000003487488,
   .00000000000000027838,
  -.00000000000000000186
};
   
static struct gsl_sf_ChebSeries by1_cs = {
  by1_data,
  13,
  -1, 1
};

int gsl_sf_bessel_Y1_impl(const double x, double * result)
{
  const double two_over_pi = 2./M_PI;
  const double xmin = 1.571*DBL_MIN; /*exp ( amax1(alog(r1mach(1)), -alog(r1mach(2)))+.01) */
  const double x_small = 2. * GSL_SQRT_MACH_EPS;
  const double xmax    = 1./GSL_MACH_EPS;
  
  if(x <= 0.) {
    return GSL_EDOM;
  }
  else if(x < xmin) {
    *result = 0.; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(x < x_small) {
    *result = two_over_pi * log(0.5*x) * gsl_sf_bessel_J1(x) +
              (0.5 + gsl_sf_cheb_eval(-1., &by1_cs))/x;
    return GSL_SUCCESS;
  }
  else if(x < 4.0) {
    *result = two_over_pi * log(0.5*x) * gsl_sf_bessel_J1(x) +
              (0.5 + gsl_sf_cheb_eval(0.125*x*x-1., &by1_cs))/x;
    return GSL_SUCCESS;
  }
  else if(x < xmax) {
    double z     = 32.0/(x*x) - 1.0;
    double ampl  = (0.75 + gsl_sf_cheb_eval(z, &_bessel_amp_phase_bm1_cs)) / sqrt(x);
    double theta = x - 3.0*M_PI_4 + gsl_sf_cheb_eval(z, &_bessel_amp_phase_bth1_cs) / x;
    *result = ampl * sin (theta);
    return GSL_SUCCESS;
  }
  else {
    *result = 0.;
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
    GSL_WARNING("gsl_sf_bessel_Y1");
  }
  return y;
}
