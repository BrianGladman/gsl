/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besi1(), besi1e() */

/* chebyshev expansions

 series for bi1        on the interval  0.	    to  9.00000d+00
					with weighted error   2.40e-17
					 log weighted error  16.62
			       significant figures required  16.23
				    decimal places required  17.14

 series for ai1        on the interval  1.25000d-01 to  3.33333d-01
					with weighted error   6.98e-17
					 log weighted error  16.16
			       significant figures required  14.53
				    decimal places required  16.82

 series for ai12       on the interval  0.	    to  1.25000d-01
				       with weighted error   3.55e-17
					log weighted error  16.45
			      significant figures required  14.69
				   decimal places required  17.12
*/

static double bi1_data[11] = {
  -.001971713261099859,
   .407348876675464810,
   .034838994299959456,
   .001545394556300123,
   .000041888521098377,
   .000000764902676483,
   .000000010042493924,
   .000000000099322077,
   .000000000000766380,
   .000000000000004741,
   .000000000000000024
};

static struct gsl_sf_cheb_series bi1_cs = {
  bi1_data,
  10,
  -1, 1,
  (double *)0,
  (double *)0
};

static double ai1_data[21] = {
  -.02846744181881479,
  -.01922953231443221,
  -.00061151858579437,
  -.00002069971253350,
   .00000858561914581,
   .00000104949824671,
  -.00000029183389184,
  -.00000001559378146,
   .00000001318012367,
  -.00000000144842341,
  -.00000000029085122,
   .00000000012663889,
  -.00000000001664947,
  -.00000000000166665,
   .00000000000124260,
  -.00000000000027315,
   .00000000000002023,
   .00000000000000730,
  -.00000000000000333,
   .00000000000000071,
  -.00000000000000006
};

static struct gsl_sf_cheb_series ai1_cs = {
  ai1_data,
  20,
  -1, 1,
  (double *)0,
  (double *)0
};

static double ai12_data[22] = {
   .02857623501828014,
  -.00976109749136147,
  -.00011058893876263,
  -.00000388256480887,
  -.00000025122362377,
  -.00000002631468847,
  -.00000000383538039,
  -.00000000055897433,
  -.00000000001897495,
   .00000000003252602,
   .00000000001412580,
   .00000000000203564,
  -.00000000000071985,
  -.00000000000040836,
  -.00000000000002101,
   .00000000000004273,
   .00000000000001041,
  -.00000000000000382,
  -.00000000000000186,
   .00000000000000033,
   .00000000000000028,
  -.00000000000000003
};

struct gsl_sf_cheb_series ai12_cs = {
  ai12_data,
  21,
  -1, 1,
  (double *)0,
  (double *)0
};


#define ROOT_EIGHT 2.82842712474619


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_I1_scaled_impl(const double x, double * result)
{
  static double xmin    = 2.0 * DBL_MIN;
  static double x_small = ROOT_EIGHT * GSL_SQRT_MACH_EPS;
  double y = fabs(x);

  if(y == 0.) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(y < xmin) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(y < x_small) {
    *result = 0.5*x;
    return GSL_SUCCESS;
  }
  else if(y <= 3.0) {
    *result = x * exp(-y) * (0.875 + gsl_sf_cheb_eval(&bi1_cs, y*y/4.5-1.0));
    return GSL_SUCCESS;
  }
  else if(y <= 8.0) {
    double b = (0.375 + gsl_sf_cheb_eval(&ai1_cs, (48.0/y-11.0)/5.0)) / sqrt(y);
    *result = ( x > 0.0 ? b : -b) /* fortran_sign(b, x) */;
    return GSL_SUCCESS;
  }
  else {
    double b = (.375 + gsl_sf_cheb_eval(&ai12_cs, 16.0/y-1.0)) / sqrt(y);
    *result = ( x > 0. ? b : -b) /* fortran_sign(b, x) */;
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_I1_impl(const double x, double * result)
{
  static double xmin    = 2.0 * DBL_MIN;
  static double x_small = ROOT_EIGHT * GSL_SQRT_MACH_EPS;
  static double xmax    = GSL_LOG_DBL_MAX; /* alog (r1mach(2)) */
  double y = fabs(x);

  if(y == 0.) {
    *result = 0.;
    return GSL_SUCCESS;
  }
  else if(y < xmin) {
    *result = 0.;
    return GSL_EUNDRFLW;
  }
  else if(y < x_small) {
    *result = 0.5*x;
    return GSL_SUCCESS;
  }
  else if(y <= 3.0) {
    *result = x * (.875 + gsl_sf_cheb_eval(&bi1_cs, y*y/4.5-1.0));
    return GSL_SUCCESS;
  }
  else if(y < xmax) {
    *result = exp(y) * gsl_sf_bessel_I1_scaled(x);
    return GSL_SUCCESS;
  }
  else {
    *result = 0.; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_I1_scaled_e(const double x, double * result)
{
  int status = gsl_sf_bessel_I1_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_I1_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_I1_e(const double x, double * result)
{
  int status = gsl_sf_bessel_I1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_I1_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_I1_scaled(const double x)
{
  double y;
  int status = gsl_sf_bessel_I1_scaled_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_I1_scaled", status);
  }
  return y;
}

double gsl_sf_bessel_I1(const double x)
{
  double y;
  int status = gsl_sf_bessel_I1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_I1", status);
  }
  return y;
}
