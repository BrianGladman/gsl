/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besk1(), besk1e() */

/* chebyshev expansions 

 series for bk1        on the interval  0.	    to  4.00000d+00
					with weighted error   7.02e-18
					 log weighted error  17.15
			       significant figures required  16.73
				    decimal places required  17.67

 series for ak1        on the interval  1.25000d-01 to  5.00000d-01
					with weighted error   6.06e-17
					 log weighted error  16.22
			       significant figures required  15.41
				    decimal places required  16.83

 series for ak12       on the interval  0.	    to  1.25000d-01
					with weighted error   2.58e-17
					 log weighted error  16.59
			       significant figures required  15.22
				    decimal places required  17.16
*/

static double bk1_data[11] = {
   0.0253002273389477705,
  -0.3531559607765448760, 
  -0.1226111808226571480, 
  -0.0069757238596398643,
  -0.0001730288957513052,
  -0.0000024334061415659,
  -0.0000000221338763073,
  -0.0000000001411488392,
  -0.0000000000006666901,
  -0.0000000000000024274,
  -0.0000000000000000070
};

static gsl_sf_cheb_series bk1_cs = {
  bk1_data,
  10,
  -1, 1,
  (double *)0,
  (double *)0,
  8
};

static double ak1_data[17] = {
   0.27443134069738830, 
   0.07571989953199368,
  -0.00144105155647540,
   0.00006650116955125,
  -0.00000436998470952,
   0.00000035402774997,
  -0.00000003311163779,
   0.00000000344597758,
  -0.00000000038989323,
   0.00000000004720819,
  -0.00000000000604783,
   0.00000000000081284,
  -0.00000000000011386,
   0.00000000000001654,
  -0.00000000000000248,
   0.00000000000000038,
  -0.00000000000000006
};
static gsl_sf_cheb_series ak1_cs = {
  ak1_data,
  16,
  -1, 1,
  (double *)0,
  (double *)0,
  9
};

static double ak12_data[14] = {
   0.06379308343739001,
   0.02832887813049721,
  -0.00024753706739052,
   0.00000577197245160,
  -0.00000020689392195,
   0.00000000973998344,
  -0.00000000055853361,
   0.00000000003732996,
  -0.00000000000282505,
   0.00000000000023720,
  -0.00000000000002176,
   0.00000000000000215,
  -0.00000000000000022,
   0.00000000000000002
};
static gsl_sf_cheb_series ak12_cs = {
  ak12_data,
  13,
  -1, 1,
  (double *)0,
  (double *)0,
  7
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_K1_scaled_impl(const double x, double * result)
{
  if(x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 2.0*GSL_DBL_MIN) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(x <= 2.0) {
    double c = gsl_sf_cheb_eval(&bk1_cs, 0.5*x*x-1.0);
    double I1;
    int stat_I1 = gsl_sf_bessel_I1_impl(x, &I1);
    *result = exp(x) * ((log(x)-M_LN2)*I1 + (0.75 + c)/x);
    return stat_I1;
  }
  else if(x <= 8.0) {
    *result = (1.25 + gsl_sf_cheb_eval(&ak1_cs, (16.0/x-5.0)/3.0)) / sqrt(x);
    return GSL_SUCCESS;
  }
  else {
    *result = (1.25 + gsl_sf_cheb_eval(&ak12_cs, 16.0/x-1.0)) / sqrt(x);
    return GSL_SUCCESS;
  }
}

int gsl_sf_bessel_K1_impl(const double x, double * result)
{
  if(x <= 0.0) {
    return GSL_EDOM;
  }
  else if(x < 2.0*GSL_DBL_MIN) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(x <= 2.0) {
    double c = gsl_sf_cheb_eval(&bk1_cs, 0.5*x*x-1.0);
    double I1;
    int stat_I1 = gsl_sf_bessel_I1_impl(x, &I1);
    *result = (log(x)-M_LN2)*I1 + (0.75 + c)/x;
    return stat_I1;
  }
  else {
    double K1_scaled;
    int stat_K1 = gsl_sf_bessel_K1_scaled_impl(x, &K1_scaled);
    int stat_e  = gsl_sf_exp_mult_impl(-x, K1_scaled, result);
    return GSL_ERROR_SELECT_2(stat_e, stat_K1);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_K1_scaled_e(const double x, double * result)
{
  int status = gsl_sf_bessel_K1_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_K1_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_K1_e(const double x, double * result)
{
  int status = gsl_sf_bessel_K1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_K1_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_K1_scaled(const double x)
{
  double y;
  int status = gsl_sf_bessel_K1_scaled_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_K1_scaled", status);
  }
  return y;
}

double gsl_sf_bessel_K1(const double x)
{
  double y;
  int status = gsl_sf_bessel_K1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_K1", status);
  }
  return y;
}
