/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* based on SLATEC besi0 */

/* chebyshev expansions

 series for bi0        on the interval  0.	    to  9.00000d+00
					with weighted error   2.46e-18
					 log weighted error  17.61
			       significant figures required  17.90
				    decimal places required  18.15

 series for ai0        on the interval  1.25000d-01 to  3.33333d-01
					with weighted error   7.87e-17
					 log weighted error  16.10
			       significant figures required  14.69
				    decimal places required  16.76


 series for ai02       on the interval  0.	    to  1.25000d-01
					with weighted error   3.79e-17
					 log weighted error  16.42
			       significant figures required  14.86
				    decimal places required  17.09
*/

static double bi0_data[12] = {
  -.07660547252839144951,
  1.927337953993808270,
   .2282644586920301339, 
   .01304891466707290428,
   .00043442709008164874,
   .00000942265768600193,
   .00000014340062895106,
   .00000000161384906966,
   .00000000001396650044,
   .00000000000009579451,
   .00000000000000053339,
   .00000000000000000245
};
static struct gsl_sf_ChebSeries bi0_cs = {
  bi0_data,
  11,
  -1, 1,
  (double *)0,
  (double *)0
};

static double ai0_data[21] = {
   .07575994494023796, 
   .00759138081082334,
   .00041531313389237,
   .00001070076463439,
  -.00000790117997921,
  -.00000078261435014,
   .00000027838499429,
   .00000000825247260,
  -.00000001204463945,
   .00000000155964859,
   .00000000022925563,
  -.00000000011916228,
   .00000000001757854,
   .00000000000112822,
  -.00000000000114684,
   .00000000000027155,
  -.00000000000002415,
  -.00000000000000608,
   .00000000000000314,
  -.00000000000000071,
   .00000000000000007
};
static struct gsl_sf_ChebSeries ai0_cs = {
  ai0_data,
  20,
  -1, 1,
  (double *)0,
  (double *)0
};

static double ai02_data[22] = {
   .05449041101410882,
   .00336911647825569,
   .00006889758346918,
   .00000289137052082,
   .00000020489185893,
   .00000002266668991,
   .00000000339623203,
   .00000000049406022,
   .00000000001188914,
  -.00000000003149915,
  -.00000000001321580,
  -.00000000000179419,
   .00000000000071801,
   .00000000000038529,
   .00000000000001539,
  -.00000000000004151,
  -.00000000000000954,
   .00000000000000382,
   .00000000000000176,
  -.00000000000000034,
  -.00000000000000027,
   .00000000000000003
};
static struct gsl_sf_ChebSeries ai02_cs = {
  ai02_data,
  21,
  -1, 1,
  (double *)0,
  (double *)0
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_I0_scaled_impl(const double x, double * result)
{
  double y = fabs(x);

  if(y < 2.0 * GSL_SQRT_MACH_EPS) {
    *result = 1.;
  }
  else if(y <= 3.0) {
    *result = exp(-y) * (2.75 + gsl_sf_cheb_eval(y*y/4.5-1., &bi0_cs));
  }
  else if(y <= 8.0) {
    *result = (.375 + gsl_sf_cheb_eval((48./y-11.)/5., &ai0_cs)) / sqrt(y);
  }
  else {
    *result = (.375 + gsl_sf_cheb_eval(16./y-1., &ai02_cs)) / sqrt(y);
  }
  return GSL_SUCCESS;
}

int gsl_sf_bessel_I0_impl(const double x, double * result)
{
  const double x_small = 2. * GSL_SQRT_MACH_EPS;
  const double xmax    = GSL_LOG_DBL_MAX - 1.;   /* alog (r1mach(2)) */
  double y = fabs(x);

  if(y < x_small) {
    *result = 1.;
    return GSL_SUCCESS;
  }
  else if(y <= 3.0) {
    *result = 2.75 + gsl_sf_cheb_eval(y*y/4.5-1.0, &bi0_cs);
    return GSL_SUCCESS;
  }
  else if(y < xmax) {
    double b_scaled;
    gsl_sf_bessel_I0_scaled_impl(x, &b_scaled);
    *result = exp(y) * b_scaled;
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_I0_scaled_e(const double x, double * result)
{
  int status = gsl_sf_bessel_I0_scaled_impl(x, result);  
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_I0_e", status);
  }
  return status;
}

int gsl_sf_bessel_I0_e(const double x, double * result)
{
  int status = gsl_sf_bessel_I0_impl(x, result);  
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_I0_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_I0_scaled(const double x)
{
  double y;
  int status = gsl_sf_bessel_I0_scaled_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_I0_scaled", status);
  }
  return y;
}


double gsl_sf_bessel_I0(const double x)
{
  double y;
  int status = gsl_sf_bessel_I0_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_I0", status);
  }
  return y;
}
