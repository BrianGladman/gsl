/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_debye.h"

/* based on MISCFUN implementations by Allan J. MacLeod  [TOMS 757] */

/* debye_n(x) = n / x^n  [Integral {0 to x} t^n/(exp(t)-1) dt] */

static double adeb1_data[19] = {
   2.40065971903814101941,
   0.19372130421893600885,
  -0.623291245548957703e-02,
   0.35111747702064800e-03,
  -0.2282224667012310e-04,
   0.158054678750300e-05,
  -0.11353781970719e-06,
   0.835833611875e-08,
  -0.62644247872e-09,
   0.4760334890e-10,
  -0.365741540e-11,
   0.28354310e-12,
  -0.2214729e-13,
   0.174092e-14,
  -0.13759e-15,
   0.1093e-16,
  -0.87e-18,
   0.7e-19,
  -0.1e-19
};
static struct gsl_sf_cheb_series adeb1_cs = {
  adeb1_data,
  18,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};

static double adeb2_data[19] = {
   2.59438102325707702826,
   0.28633572045307198337,
  -0.1020626561580467129e-01,
   0.60491097753468435e-03,
  -0.4052576589502104e-04,
   0.286338263288107e-05,
  -0.20863943030651e-06,
   0.1552378758264e-07,
  -0.117312800866e-08,
   0.8973585888e-10,
  -0.693176137e-11,
   0.53980568e-12,
  -0.4232405e-13,
   0.333778e-14,
  -0.26455e-15,
   0.2106e-16,
  -0.168e-17,
   0.13e-18,
  -0.1e-19
};
static struct gsl_sf_cheb_series adeb2_cs = {
  adeb2_data,
  18,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};

static double adeb3_data[19] = {
   2.70773706832744094526,
   0.34006813521109175100,
  -0.1294515018444086863e-01,
   0.79637553801738164e-03,
  -0.5463600095908238e-04,
   0.392430195988049e-05,
  -0.28940328235386e-06,
   0.2173176139625e-07,
  -0.165420999498e-08,
   0.12727961892e-09,
  -0.987963459e-11,
   0.77250740e-12,
  -0.6077972e-13,
   0.480759e-14,
  -0.38204e-15,
   0.3048e-16,
  -0.244e-17,
   0.20e-18,
  -0.2e-19
};
static struct gsl_sf_cheb_series adeb3_cs = {
  adeb3_data,
  18,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};

static double adeb4_data[19] = {
   2.78186941502052346008,
   0.37497678352689286364,
  -0.1494090739903158326e-01,
   0.94567981143704274e-03,
  -0.6613291613893255e-04,
   0.481563298214449e-05,
  -0.35880839587593e-06,
   0.2716011874160e-07,
  -0.208070991223e-08,
   0.16093838692e-09,
  -0.1254709791e-10,
   0.98472647e-12,
  -0.7772369e-13,
   0.616483e-14,
  -0.49107e-15,
   0.3927e-16,
  -0.315e-17,
   0.25e-18,
  -0.2e-19
};
static struct gsl_sf_cheb_series adeb4_cs = {
  adeb4_data,
  18,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};


int gsl_sf_debye_1_impl(const double x, double * result)
{
  const double DEBINF = 0.60792710185402662866;
  const double xlim   = -GSL_LOG_DBL_MIN;
  const double xup    = -(M_LN2 + GSL_LOG_MACH_EPS);

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 2.0*GSL_SQRT_MACH_EPS) {
    *result = ((x-9.0) * x + 36.0) / 36.0;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double t = x*x/8.0 - 1.0;
    *result = gsl_sf_cheb_eval(t, &adeb1_cs) - 0.25 * x;
    return GSL_SUCCESS;
  }
  else if(x < xup) {
    int i;
    int nexp = floor(xlim/x);
    double expmx = exp(-x);
    double sum   = 0.0;
    double xk    = nexp * x;
    double rk    = nexp;
    double t;
    for(i=nexp; i>=1; i--) {
      t = (1.0 + 1.0/xk)/rk;
      sum = sum*expmx + t;
      rk -= 1.0;
      xk -= x;
    }
    *result = 1.0/(x*DEBINF) - sum*expmx;
  }
  else if(x < xlim) {
    *result = 1.0/(x*DEBINF) - exp(-x)*(1.0 + 1.0/x);
    return GSL_SUCCESS;
  }
  else {
    *result = 1.0/(x*DEBINF);
    return GSL_SUCCESS;
  }
} 
     
int gsl_sf_debye_2_impl(const double x, double * result)
{
  const double DEBINF = 4.80822761263837714160;
  const double xlo = 2.0*M_SQRT2*GSL_SQRT_MACH_EPS;
  const double xup = -(M_LN2 + GSL_LOG_MACH_EPS);
  const double xlim1  = -GSL_LOG_DBL_MIN;
  const double xlim2  = sqrt(DEBINF) / GSL_SQRT_DBL_MIN;

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < xlo) {
    *result = ((x-8.0)*x + 24.0) / 24.0;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double t = x*x/8.0 - 1.0;
    *result = gsl_sf_cheb_eval(t, &adeb2_cs) - x/3.0;
    return GSL_SUCCESS;
  }
  else if(x < xup) {
    int i;
    int nexp = floor(xlim1/x);
    double expmx = exp(-x);
    double xk    = nexp * x;
    double rk    = nexp;
    double sum   = 0.0;
    double t;
    for(i=nexp; i>=1; i--) {
      t = (1.0 + 2.0/xk + 2.0/(xk*xk)) / rk;
      sum = sum*expmx + t;
      rk -= 1.0;
      xk -= x;
    }
    *result = DEBINF /(x*x) - 2.0 * sum * expmx;
    return GSL_SUCCESS;
  }
  else if(x < xlim1) {
    double sum = ((x+2.0)*x + 2.0) / (x*x);
    *result = DEBINF /(x*x) - 2.0 * sum * exp(-x);
    return GSL_SUCCESS;
  }
  else if(x < xlim2) {
    *result = DEBINF /(x*x);
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
}

int gsl_sf_debye_3_impl(const double x, double * result)
{
  const double DEBINF = 0.51329911273421675946e-01;
  const double xlo = 2.0*M_SQRT2*GSL_SQRT_MACH_EPS;
  const double xup = -(M_LN2 + GSL_LOG_MACH_EPS);
  const double xlim1 = -GSL_LOG_DBL_MIN;
  const double xlim2 = pow(1.0/DEBINF,1.0/3.0) / GSL_ROOT3_DBL_MIN;

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < xlo) {
    *result = ((x-7.5)*x + 20.0) / 20.0;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double t = x*x/8.0 - 1.0;
    *result = gsl_sf_cheb_eval(t, &adeb3_cs) - 0.375*x;
    return GSL_SUCCESS;
  }
  else if(x < xup) {
    int i;
    int nexp = floor(xlim1/x);
    double expmx = exp(-x);
    double xk    = nexp * x;
    double rk    = nexp;
    double sum   = 0.0;
    double t, xk_inv;
    for(i=nexp; i>=1; i--) {
      xk_inv = 1.0/xk;
      t = (((6.0*xk_inv + 6.0)*xk_inv + 3.0)*xk_inv + 1.0) / rk;
      sum = sum*expmx + t;
      rk -= 1.0;
      xk -= x;
    }
    *result = 1.0/(DEBINF * x*x*x) - 3.0 * sum * expmx;
    return GSL_SUCCESS;
  }
  else if(x < xlim1) {
    double sum = (((x + 3.0)*x + 6.0)*x + 6.0) / (x*x*x);
    *result = 1.0/(DEBINF * x*x*x) - 3.0 * sum * exp(-x);
    return GSL_SUCCESS;
  }
  else if(x < xlim2) {
    *result = 1.0/(DEBINF * x*x*x);
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
}

int gsl_sf_debye_4_impl(const double x, double * result)
{
  const double DEBINF = 99.54506449376351292781;
  const double xlo = 2.0*M_SQRT2*GSL_SQRT_MACH_EPS;
  const double xup = -(M_LN2 + GSL_LOG_MACH_EPS);
  const double xlim1 = -GSL_LOG_DBL_MIN;
  const double xlim2 = pow(DEBINF, 0.25) / GSL_ROOT4_DBL_MIN;
  
  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < xlo) {
    *result = ((2.5*x - 18.0)*x + 45.0) / 45.0;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double t = x*x/8.0 - 1.0;
    *result = gsl_sf_cheb_eval(t, &adeb4_cs) - 2.0*x/5.0;
    return GSL_SUCCESS;
  }
  else if(x < xup) {
    int i;
    int nexp = floor(xlim1/x);
    double expmx = exp(-x);
    double xk    = nexp * x;
    double rk    = nexp;
    double sum   = 0.0;
    double t, xk_inv;
    for(i=nexp; i>=1; i--) {
      xk_inv = 1.0/xk;
      t = ((((24.0*xk_inv + 24.0)*xk_inv + 12.0)*xk_inv + 4.0)*xk_inv + 1.0) / rk;
      sum = sum*expmx + t;
      rk -= 1.0;
      xk -= x;
    }
    *result = DEBINF/(x*x*x*x) - 4.0 * sum * expmx;
    return GSL_SUCCESS;
  }
  else if(x < xlim1) {
    double sum = ((((x + 4.0)*x + 12.0)*x + 24.0)*x + 24.0) / (x*x*x*x);
    *result = DEBINF/(x*x*x*x) - 4.0 * sum * exp(-x);
    return GSL_SUCCESS;
  }
  else if(x < xlim2) {
    double x2 = x*x;
    *result = (DEBINF/x2)/x2;
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
}
