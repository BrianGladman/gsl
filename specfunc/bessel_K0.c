/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

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
  -0.03532739323390276872,
   0.3442898999246284869, 
   0.03597993651536150163,
   0.00126461541144692592,
   0.00002286212103119451,
   0.00000025347910790261,
   0.00000000190451637722,
   0.00000000001034969525,
   0.00000000000004259816,
   0.00000000000000013744,
   0.00000000000000000035
};
static gsl_sf_cheb_series bk0_cs = {
  bk0_data,
  10,
  -1, 1,
  (double *)0,
  (double *)0,
  10
};

static double ak0_data[17] = {
  -0.07643947903327941,
  -0.02235652605699819,
   0.00077341811546938,
  -0.00004281006688886,
   0.00000308170017386,
  -0.00000026393672220,
   0.00000002563713036,
  -0.00000000274270554,
   0.00000000031694296,
  -0.00000000003902353,
   0.00000000000506804,
  -0.00000000000068895,
   0.00000000000009744,
  -0.00000000000001427,
   0.00000000000000215,
  -0.00000000000000033,
   0.00000000000000005
};
static gsl_sf_cheb_series ak0_cs = {
  ak0_data,
  16,
  -1, 1,
  (double *)0,
  (double *)0,
  10
};

static double ak02_data[14] = {
  -0.01201869826307592,
  -0.00917485269102569,
   0.00014445509317750,
  -0.00000401361417543,
   0.00000015678318108,
  -0.00000000777011043,
   0.00000000046111825,
  -0.00000000003158592,
   0.00000000000243501,
  -0.00000000000020743,
   0.00000000000001925,
  -0.00000000000000192,
   0.00000000000000020,
  -0.00000000000000002
};
static gsl_sf_cheb_series ak02_cs = {
  ak02_data,
  13,
  -1, 1,
  (double *)0,
  (double *)0,
  8
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_K0_scaled_impl(const double x, double * result)
{
  if(x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x <= 2.0) {
    double c = gsl_sf_cheb_eval(&bk0_cs, 0.5*x*x-1.0);
    double I0;
    int stat_I0 = gsl_sf_bessel_I0_impl(x, &I0);
    *result = exp(x) * ((-log(x)+M_LN2)*I0 - 0.25 + c);
    return stat_I0;
  }
  else if(x <= 8.0) {
    *result = (1.25 + gsl_sf_cheb_eval(&ak0_cs, (16.0/x-5.0)/3.0)) / sqrt(x);
    return GSL_SUCCESS;
  }
  else {
    *result = (1.25 + gsl_sf_cheb_eval(&ak02_cs, 16.0/x-1.0)) / sqrt(x);
    return GSL_SUCCESS;
  } 
}


int gsl_sf_bessel_K0_impl(const double x, double * result)
{
  if(x <= 0.0) {
    return GSL_EDOM;
  }
  else if(x <= 2.0) {
    double c = gsl_sf_cheb_eval(&bk0_cs, 0.5*x*x-1.0);
    double I0;
    int stat_I0 = gsl_sf_bessel_I0_impl(x, &I0);
    *result = (-log(x)+M_LN2)*I0 - 0.25 + c;
    return stat_I0;
  }
  else {
    double K0_scaled;
    int stat_K0 = gsl_sf_bessel_K0_scaled_impl(x, &K0_scaled);
    int stat_e  = gsl_sf_exp_mult_impl(-x, K0_scaled, result);
    return GSL_ERROR_SELECT_2(stat_e, stat_K0);
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
