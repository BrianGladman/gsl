/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#define ROOT_EIGHT (2.0*M_SQRT2)


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
  -0.001971713261099859,
   0.407348876675464810,
   0.034838994299959456,
   0.001545394556300123,
   0.000041888521098377,
   0.000000764902676483,
   0.000000010042493924,
   0.000000000099322077,
   0.000000000000766380,
   0.000000000000004741,
   0.000000000000000024
};
static gsl_sf_cheb_series bi1_cs = {
  bi1_data,
  10,
  -1, 1,
  (double *)0,
  (double *)0,
  10
};

static double ai1_data[21] = {
  -0.02846744181881479,
  -0.01922953231443221,
  -0.00061151858579437,
  -0.00002069971253350,
   0.00000858561914581,
   0.00000104949824671,
  -0.00000029183389184,
  -0.00000001559378146,
   0.00000001318012367,
  -0.00000000144842341,
  -0.00000000029085122,
   0.00000000012663889,
  -0.00000000001664947,
  -0.00000000000166665,
   0.00000000000124260,
  -0.00000000000027315,
   0.00000000000002023,
   0.00000000000000730,
  -0.00000000000000333,
   0.00000000000000071,
  -0.00000000000000006
};
static gsl_sf_cheb_series ai1_cs = {
  ai1_data,
  20,
  -1, 1,
  (double *)0,
  (double *)0,
  11
};

static double ai12_data[22] = {
   0.02857623501828014,
  -0.00976109749136147,
  -0.00011058893876263,
  -0.00000388256480887,
  -0.00000025122362377,
  -0.00000002631468847,
  -0.00000000383538039,
  -0.00000000055897433,
  -0.00000000001897495,
   0.00000000003252602,
   0.00000000001412580,
   0.00000000000203564,
  -0.00000000000071985,
  -0.00000000000040836,
  -0.00000000000002101,
   0.00000000000004273,
   0.00000000000001041,
  -0.00000000000000382,
  -0.00000000000000186,
   0.00000000000000033,
   0.00000000000000028,
  -0.00000000000000003
};
static gsl_sf_cheb_series ai12_cs = {
  ai12_data,
  21,
  -1, 1,
  (double *)0,
  (double *)0,
  9
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_I1_scaled_impl(const double x, double * result)
{
  const double xmin    = 2.0 * GSL_DBL_MIN;
  const double x_small = ROOT_EIGHT * GSL_SQRT_DBL_EPSILON;
  const double y = fabs(x);

  if(y == 0.0) {
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
    *result = ( x > 0.0 ? b : -b);
    return GSL_SUCCESS;
  }
  else {
    double b = (0.375 + gsl_sf_cheb_eval(&ai12_cs, 16.0/y-1.0)) / sqrt(y);
    *result = ( x > 0.0 ? b : -b);
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_I1_impl(const double x, double * result)
{
  const double xmin    = 2.0 * GSL_DBL_MIN;
  const double x_small = ROOT_EIGHT * GSL_SQRT_DBL_EPSILON;
  const double y = fabs(x);

  if(y == 0.0) {
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
    *result = x * (0.875 + gsl_sf_cheb_eval(&bi1_cs, y*y/4.5-1.0));
    return GSL_SUCCESS;
  }
  else if(y < GSL_LOG_DBL_MAX) {
    *result = exp(y) * gsl_sf_bessel_I1_scaled(x);
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0; /* FIXME: should be Inf */
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
