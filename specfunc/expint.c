/* Author: G. Jungman
 * RCS: $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_expint.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/*
 Chebyshev expansions: based on SLATEC e1.f, W. Fullerton
 
 Series for AE11       on the interval -1.00000D-01 to  0.
					with weighted error   1.76E-17
					 log weighted error  16.75
			       significant figures required  15.70
				    decimal places required  17.55


 Series for AE12       on the interval -2.50000D-01 to -1.00000D-01
					with weighted error   5.83E-17
					 log weighted error  16.23
			       significant figures required  15.76
				    decimal places required  16.93


 Series for E11        on the interval -4.00000D+00 to -1.00000D+00
					with weighted error   1.08E-18
					 log weighted error  17.97
			       significant figures required  19.02
				    decimal places required  18.61


 Series for E12        on the interval -1.00000D+00 to  1.00000D+00
					with weighted error   3.15E-18
					 log weighted error  17.50
			approx significant figures required  15.8
				    decimal places required  18.10


 Series for AE13       on the interval  2.50000D-01 to  1.00000D+00
					with weighted error   2.34E-17
					 log weighted error  16.63
			       significant figures required  16.14
				    decimal places required  17.33


 Series for AE14       on the interval  0.	    to  2.50000D-01
					with weighted error   5.41E-17
					 log weighted error  16.27
			       significant figures required  15.38
				    decimal places required  16.97

   770701  DATE WRITTEN
   890531  Changed all specific intrinsics to generic.  (WRB)
   891115  Modified prologue description.  (WRB)
   891115  REVISION DATE from Version 3.2
   891214  Prologue converted to Version 4.0 format.  (BAB)
   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
   920618  Removed space from variable names.  (RWC, WRB)
*/

static double AE11_data[39] = {
   .121503239716065790,
  -.065088778513550150,
   .004897651357459670,
  -.000649237843027216,
   .000093840434587471,
   .000000420236380882,
  -.000008113374735904,
   .000002804247688663,
   .000000056487164441,
  -.000000344809174450,
   .000000058209273578,
   .000000038711426349,
  -.000000012453235014,
  -.000000005118504888,
   .000000002148771527,
   .000000000868459898,
  -.000000000343650105,
  -.000000000179796603,
   .000000000047442060,
   .000000000040423282,
  -.000000000003543928,
  -.000000000008853444,
  -.000000000000960151,
   .000000000001692921,
   .000000000000607990,
  -.000000000000224338,
  -.000000000000200327,
  -.000000000000006246,
   .000000000000045571,
   .000000000000016383,
  -.000000000000005561,
  -.000000000000006074,
  -.000000000000000862,
   .000000000000001223,
   .000000000000000716,
  -.000000000000000024,
  -.000000000000000201,
  -.000000000000000082,
   .000000000000000017
};

static double AE12_data[25] = {
   .582417495134726740,
  -.158348850905782750,
  -.006764275590323141,
   .005125843950185725,
   .000435232492169391,
  -.000143613366305483,
  -.000041801320556301,
  -.000002713395758640,
   .000001151381913647,
   .000000420650022012,
   .000000066581901391,
   .000000000662143777,
  -.000000002844104870,
  -.000000000940724197,
  -.000000000177476602,
  -.000000000015830222,
   .000000000002905732,
   .000000000001769356,
   .000000000000492735,
   .000000000000093709,
   .000000000000010707,
  -.000000000000000537,
  -.000000000000000716,
  -.000000000000000244,
  -.000000000000000058
};

static double E11_data[19] = {
  -16.11346165557149402600,
    7.79407277874268027690,
   -1.95540581886314195070,
     .37337293866277945612,
    -.05692503191092901938,
     .00721107776966009185,
    -.00078104901449841593,
     .00007388093356262168,
    -.00000620286187580820,
     .00000046816002303176,
    -.00000003209288853329,
     .00000000201519974874,
    -.00000000011673686816,
     .00000000000627627066,
    -.00000000000031481541,
     .00000000000001479904,
    -.00000000000000065457,
     .00000000000000002733,
    -.00000000000000000108
};

static double E12_data[16] = {
  -0.03739021479220279500,
   0.04272398606220957700,
   -.13031820798497005440,
    .01441912402469889073,
   -.00134617078051068022,
    .00010731029253063780,
   -.00000742999951611943,
    .00000045377325690753,
   -.00000002476417211390,
    .00000000122076581374,
   -.00000000005485141480,
    .00000000000226362142,
   -.00000000000008635897,
    .00000000000000306291,
   -.00000000000000010148,
    .00000000000000000315
};

static double AE13_data[25] = {
  -.605773246640603460,
  -.112535243483660900,
   .013432266247902779,
  -.001926845187381145,
   .000309118337720603,
  -.000053564132129618,
   .000009827812880247,
  -.000001885368984916,
   .000000374943193568,
  -.000000076823455870,
   .000000016143270567,
  -.000000003466802211,
   .000000000758754209,
  -.000000000168864333,
   .000000000038145706,
  -.000000000008733026,
   .000000000002023672,
  -.000000000000474132,
   .000000000000112211,
  -.000000000000026804,
   .000000000000006457,
  -.000000000000001568,
   .000000000000000383,
  -.000000000000000094,
   .000000000000000023
};

static double AE14_data[26] = {
  -.18929180007530170,
  -.08648117855259871,
   .00722410154374659,
  -.00080975594575573,
   .00010999134432661,
  -.00001717332998937,
   .00000298562751447,
  -.00000056596491457,
   .00000011526808397,
  -.00000002495030440,
   .00000000569232420,
  -.00000000135995766,
   .00000000033846628,
  -.00000000008737853,
   .00000000002331588,
  -.00000000000641148,
   .00000000000181224,
  -.00000000000052538,
   .00000000000015592,
  -.00000000000004729,
   .00000000000001463,
  -.00000000000000461,
   .00000000000000148,
  -.00000000000000048,
   .00000000000000016,
  -.00000000000000005
};

static struct gsl_sf_cheb_series AE11_cs = {
  AE11_data,
  38,
  -1, 1,
  (double *)0,
  (double *)0
};

static struct gsl_sf_cheb_series AE12_cs = {
  AE12_data,
  24,
  -1, 1,
  (double *)0,
  (double *)0
};

static struct gsl_sf_cheb_series E11_cs = {
  E11_data,
  18,
  -1, 1,
  (double *)0,
  (double *)0
};

static struct gsl_sf_cheb_series E12_cs = {
  E12_data,
  15,
  -1, 1,
  (double *)0,
  (double *)0
};

static struct gsl_sf_cheb_series AE13_cs = {
  AE13_data,
  24,
  -1, 1,
  (double *)0,
  (double *)0
};

static struct gsl_sf_cheb_series AE14_cs = {
  AE14_data,
  25,
  -1, 1,
  (double *)0,
  (double *)0
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* checked OK [GJ] */
int gsl_sf_expint_E1_impl(const double x, double * result)
{
  const double xmaxt = -GSL_LOG_DBL_MIN;      /* XMAXT = -LOG (R1MACH(1)) */
  const double xmax  = xmaxt - log(xmaxt);    /* XMAX = XMAXT - LOG(XMAXT) */

  if(x < -xmax) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(x <= -10.) {
    *result = exp(-x)/x * (1. + gsl_sf_cheb_eval(20./x+1., &AE11_cs));
    return GSL_SUCCESS;
  }
  else if(x <= -4.) {
    *result = exp(-x)/x * (1. + gsl_sf_cheb_eval((40./x+7.)/3., &AE12_cs));
    return GSL_SUCCESS;
  }
  else if(x <= -1.) {
    *result = -log(fabs(x)) + gsl_sf_cheb_eval((2.*x+5.)/3., &E11_cs);
    return GSL_SUCCESS;
  }
  else if(x == 0.) {
    return GSL_EDOM;
  }
  else if(x <= 1.) {
    *result = -log(fabs(x)) - 0.6875 + x + gsl_sf_cheb_eval(x, &E12_cs);
    return GSL_SUCCESS;
  }
  else if(x <= 4.) {
    *result = exp(-x)/x * (1. + gsl_sf_cheb_eval((8./x-5.)/3., &AE13_cs));
    return GSL_SUCCESS;
  }
  else if(x <= xmax) {
    *result = exp(-x)/x * (1. +  gsl_sf_cheb_eval(8./x-1., &AE14_cs));
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
}

int gsl_sf_expint_E2_impl(const double x, double * result)
{
  const double xmaxt = -GSL_LOG_DBL_MIN;
  const double xmax  = xmaxt - log(xmaxt);

  if(x < -xmax) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(x <= xmax) {
    double E1;
    int stat_E1 = gsl_sf_expint_E1_impl(x, &E1);
    *result = 0.5 * (exp(-x) - x*E1);
    return stat_E1;
  }
  else {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
}


/* checked OK [GJ] */
int gsl_sf_expint_Ei_impl(const double x, double * result)
{
  int status = gsl_sf_expint_E1_impl(-x, result);
  if(status == GSL_SUCCESS) {
    *result = - *result;
  }
  return status;
}

#if 0
static double recurse_En(int n, double x, double E1)
{
  int i;
  double En = E1;
  double ex = exp(-x);
  for(i=2; i<=n; i++) {
    En = 1./(i-1) * (ex - x * En);
  }
  return En;
}
#endif

/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_expint_E1_e(const double x, double * result)
{
  int status = gsl_sf_expint_E1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_expint_E1_e", status);
  }
  return status;
}

int gsl_sf_expint_E2_e(const double x, double * result)
{
  int status = gsl_sf_expint_E2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_expint_E2_e", status);
  }
  return status;
}

int gsl_sf_expint_Ei_e(const double x, double * result)
{
  int status = gsl_sf_expint_Ei_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_expint_Ei_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_expint_E1(double x)
{
  double y;
  int status = gsl_sf_expint_E1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_expint_E1", status);
  }
  return y;
}

double gsl_sf_expint_E2(double x)
{
  double y;
  int status = gsl_sf_expint_E2_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_expint_E2", status);
  }
  return y;
}

double gsl_sf_expint_Ei(double x)
{
  double y;
  int status = gsl_sf_expint_Ei_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_expint_Ei", status);
  }
  return y;
}
