/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC bk0(), bk0e() */

/* chebyshev expansions 

 series for bk0        on the interval  0.	    to  4.00000d+00
					with weighted error   3.57e-19
					 log weighted error  18.45
			       significant figures required  17.99
				    decimal places required  18.97

 series for ak0        on the interval  1.25000d-01 to  5.00000d-01
					with weighted error   5.34e-17
					 log weighted error  16.27
			       significant figures required  14.92
				    decimal places required  16.89

 series for ak02       on the interval  0.	    to  1.25000d-01
					with weighted error   2.34e-17
					 log weighted error  16.63
			       significant figures required  14.67
				    decimal places required  17.20
*/

static double bk0_data[11] = {
  -.03532739323390276872,
   .3442898999246284869, 
   .03597993651536150163,
   .00126461541144692592,
   .00002286212103119451,
   .00000025347910790261,
   .00000000190451637722,
   .00000000001034969525,
   .00000000000004259816,
   .00000000000000013744,
   .00000000000000000035
};

static struct gsl_sf_ChebSeries bk0_cs = {
  bk0_data,
  10,
  -1, 1
};

static double ak0_data[17] = {
  -.07643947903327941,
  -.02235652605699819,
   .00077341811546938,
  -.00004281006688886,
   .00000308170017386,
  -.00000026393672220,
   .00000002563713036,
  -.00000000274270554,
   .00000000031694296,
  -.00000000003902353,
   .00000000000506804,
  -.00000000000068895,
   .00000000000009744,
  -.00000000000001427,
   .00000000000000215,
  -.00000000000000033,
   .00000000000000005
};

static struct gsl_sf_ChebSeries ak0_cs = {
  ak0_data,
  16,
  -1, 1
};

static double ak02_data[14] = {
  -.01201869826307592,
  -.00917485269102569,
   .00014445509317750,
  -.00000401361417543,
   .00000015678318108,
  -.00000000777011043,
   .00000000046111825,
  -.00000000003158592,
   .00000000000243501,
  -.00000000000020743,
   .00000000000001925,
  -.00000000000000192,
   .00000000000000020,
  -.00000000000000002
};

static struct gsl_sf_ChebSeries ak02_cs = {
  ak02_data,
  13,
  -1, 1
};

int gsl_sf_bessel_K0_scaled_impl(const double x, double * result)
{
  const double x_small = 2.0 * GSL_SQRT_MACH_EPS;

  if(x <= 0.0) {
    return GSL_EDOM;
  }
  else if(x < x_small) {
    *result = exp(x) * (-log(0.5*x)*gsl_sf_bessel_I0(x) - .25
      	      	        + gsl_sf_cheb_eval(-1., &bk0_cs)
      	      	       );
    return GSL_SUCCESS;
  }
  else if(x <= 2.) {
    double y = x*x;
    *result = exp(x) * (-log(0.5*x)*gsl_sf_bessel_I0(x) - .25
      	      	        + gsl_sf_cheb_eval(0.5*y-1., &bk0_cs)
      	      	       );
    return GSL_SUCCESS;
  }
  else if(x <= 8.) {
    *result = (1.25 + gsl_sf_cheb_eval((16./x-5.)/3., &ak0_cs)) / sqrt(x);
    return GSL_SUCCESS;
  }
  else {
    *result = (1.25 + gsl_sf_cheb_eval(16./x-1., &ak02_cs)) / sqrt(x);
    return GSL_SUCCESS;
  } 
}


int gsl_sf_bessel_K0_impl(const double x, double * result)
{
  const double x_small = 2.*GSL_SQRT_MACH_EPS;
  const double xmax = GSL_LOG_DBL_MAX - 0.5 * 6.9 /* FIXME: ?? */  - 0.01;
    /*
      xmax = -alog(r1mach(1))
      xmax = xmax - 0.5*xmax*alog(xmax)/(xmax+0.5) - 0.01
    */

  if(x <= 0.) {
    return GSL_EDOM;
  }
  else if(x < x_small) {
    *result = -log(0.5*x)*gsl_sf_bessel_I0(x) - .25 + gsl_sf_cheb_eval(-1., &bk0_cs);
    return GSL_SUCCESS;
  }
  else if(x <= 2.) {
    double y = x*x;
    *result = -log(0.5*x)*gsl_sf_bessel_I0(x) - .25 + gsl_sf_cheb_eval(0.5*y-1., &bk0_cs);
    return GSL_SUCCESS;
  }
  else if(x < xmax) {
    *result = exp(-x) * gsl_sf_bessel_K0_scaled(x);
    return GSL_SUCCESS;
  }
  else {
    return GSL_EUNDRFLW;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_K0_scaled_e(const double x, double * result)
{
  int status = gsl_sf_bessel_K0_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_K0_scaled_e", status);
  }
  return status;
}


int gsl_sf_bessel_K0_e(const double x, double * result)
{
  int status = gsl_sf_bessel_K0_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_K0_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_K0_scaled(const double x)
{
  double y;
  int status = gsl_sf_bessel_K0_scaled_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_K0_scaled", status);
  }
  return y;
}


double gsl_sf_bessel_K0(const double x)
{
  double y;
  int status = gsl_sf_bessel_K0_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_K0", status);
  }
  return y;
}
