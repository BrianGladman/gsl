/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_debye.h"


static double adeb1_data[17] = {
   2.4006597190381410194,
   0.1937213042189360089,
  -0.62329124554895770e-02,
   0.3511174770206480e-03,
  -0.228222466701231e-04,
   0.15805467875030e-05,
  -0.1135378197072e-06,
   0.83583361188e-08,
  -0.6264424787e-09,
   0.476033489e-10,
  -0.36574154e-11,
   0.2835431e-12,
  -0.221473e-13,
   0.17409e-14,
  -0.1376e-15,
   0.109e-16,
  -0.9e-18
};
static gsl_sf_cheb_series adeb1_cs = {
  adeb1_data,
  16,
  -1.0, 1.0,
  (double *)0,
  (double *)0,
  9
};

static double adeb2_data[18] = {
   2.5943810232570770282,
   0.2863357204530719834,
  -0.102062656158046713e-01,
   0.6049109775346844e-03,
  -0.405257658950210e-04,
   0.28633826328811e-05,
  -0.2086394303065e-06,
   0.155237875826e-07,
  -0.11731280087e-08,
   0.897358589e-10,
  -0.69317614e-11,
   0.5398057e-12,
  -0.423241e-13,
   0.33378e-14,
  -0.2645e-15,
   0.211e-16,
  -0.17e-17,
   0.1e-18
};
static gsl_sf_cheb_series adeb2_cs = {
  adeb2_data,
  17,
  -1.0, 1.0,
  (double *)0,
  (double *)0,
  10
};

static double adeb3_data[17] = {
   2.707737068327440945,
   0.340068135211091751,
  -0.12945150184440869e-01,
   0.7963755380173816e-03,
  -0.546360009590824e-04,
   0.39243019598805e-05,
  -0.2894032823539e-06,
   0.217317613962e-07,
  -0.16542099950e-08,
   0.1272796189e-09,
  -0.987963460e-11,
   0.7725074e-12,
  -0.607797e-13,
   0.48076e-14,
  -0.3820e-15,
   0.305e-16,
  -0.24e-17
};
static gsl_sf_cheb_series adeb3_cs = {
  adeb3_data,
  16,
  -1.0, 1.0,
  (double *)0,
  (double *)0,
  10
};

static double adeb4_data[17] = {
   2.781869415020523460,
   0.374976783526892863,
  -0.14940907399031583e-01,
   0.945679811437042e-03,
  -0.66132916138933e-04,
   0.4815632982144e-05,
  -0.3588083958759e-06,
   0.271601187416e-07,
  -0.20807099122e-08,
   0.1609383869e-09,
  -0.125470979e-10,
   0.9847265e-12,
  -0.777237e-13,
   0.61648e-14,
  -0.4911e-15,
   0.393e-16,
  -0.32e-17
};
static gsl_sf_cheb_series adeb4_cs = {
  adeb4_data,
  16,
  -1.0, 1.0,
  (double *)0,
  (double *)0,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_debye_1_impl(const double x, gsl_sf_result * result)
{
  const double val_infinity = 1.64493406684822644;
  const double xcut = -GSL_LOG_DBL_MIN;

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x < 2.0*GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0 - 0.25*x + x*x/36.0;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    const double t = x*x/8.0 - 1.0;
    gsl_sf_result c;
    gsl_sf_cheb_eval_impl(&adeb1_cs, t, &c);
    result->val = c.val - 0.25 * x;
    result->err = c.err + 0.25 * x * GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < -(M_LN2 + GSL_LOG_DBL_EPSILON)) {
    const int nexp = floor(xcut/x);
    const double ex  = exp(-x);
    double sum = 0.0;
    double xk  = nexp * x;
    double rk  = nexp;
    int i;
    for(i=nexp; i>=1; i--) {
      sum *= ex;
      sum += (1.0 + 1.0/xk)/rk;
      rk -= 1.0;
      xk -= x;
    }
    result->val = val_infinity/x - sum*ex;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < xcut) {
    result->val = (val_infinity - exp(-x)*(x+1.0)) / x;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = val_infinity/x;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

    
int gsl_sf_debye_2_impl(const double x, gsl_sf_result * result)
{
  const double val_infinity = 4.80822761263837714;
  const double xcut = -GSL_LOG_DBL_MIN;

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x < 2.0*M_SQRT2*GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0 - x/3.0 + x*x/24.0;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    const double t = x*x/8.0 - 1.0;
    gsl_sf_result c;
    gsl_sf_cheb_eval_impl(&adeb2_cs, t, &c);
    result->val = c.val - x/3.0;
    result->err = c.err + GSL_DBL_EPSILON * x/3.0;
    return GSL_SUCCESS;
  }
  else if(x < -(M_LN2 + GSL_LOG_DBL_EPSILON)) {
    const int nexp = floor(xcut/x);
    const double ex  = exp(-x);
    double xk  = nexp * x;
    double rk  = nexp;
    double sum = 0.0;
    int i;
    for(i=nexp; i>=1; i--) {
      sum *= ex;
      sum += (1.0 + 2.0/xk + 2.0/(xk*xk)) / rk;
      rk -= 1.0;
      xk -= x;
    }
    result->val = val_infinity/(x*x) - 2.0 * sum * ex;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < xcut) {
    const double x2  = x*x;
    const double sum = 2.0 + 2.0*x + x2;
    result->val = (val_infinity - 2.0 * sum * exp(-x)) / x2;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = (val_infinity/x)/x;
    result->err = GSL_DBL_EPSILON * result->val;
    if(result->val == 0.0)
      return GSL_EUNDRFLW;
    else
      return GSL_SUCCESS;
  }
}


int gsl_sf_debye_3_impl(const double x, gsl_sf_result * result)
{
  const double val_infinity = 19.4818182068004875;
  const double xcut = -GSL_LOG_DBL_MIN;

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x < 2.0*M_SQRT2*GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0 - 3.0*x/8.0 + x*x/20.0;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    const double t = x*x/8.0 - 1.0;
    gsl_sf_result c;
    gsl_sf_cheb_eval_impl(&adeb3_cs, t, &c);
    result->val = c.val - 0.375*x;
    result->err = c.err + GSL_DBL_EPSILON * 0.375*x;
    return GSL_SUCCESS;
  }
  else if(x < -(M_LN2 + GSL_LOG_DBL_EPSILON)) {
    const int nexp = floor(xcut/x);
    const double ex  = exp(-x);
    double xk  = nexp * x;
    double rk  = nexp;
    double sum = 0.0;
    int i;
    for(i=nexp; i>=1; i--) {
      double xk_inv = 1.0/xk;
      sum *= ex;
      sum += (((6.0*xk_inv + 6.0)*xk_inv + 3.0)*xk_inv + 1.0) / rk;
      rk -= 1.0;
      xk -= x;
    }
    result->val = val_infinity/(x*x*x) - 3.0 * sum * ex;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x < xcut) {
    const double x3 = x*x*x;
    const double sum = 6.0 + 6.0*x + 3.0*x*x + x3;
    result->val = (val_infinity - 3.0 * sum * exp(-x)) / x3;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    result->val = ((val_infinity/x)/x)/x;
    result->err = GSL_DBL_EPSILON * result->val;
    if(result->val == 0.0)
      return GSL_EUNDRFLW;
    else
      return GSL_SUCCESS;
  }
}


int gsl_sf_debye_4_impl(const double x, gsl_sf_result * result)
{
  const double val_infinity = 99.5450644937635129;
  const double xcut = -GSL_LOG_DBL_MIN;

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x < 2.0*M_SQRT2*GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0 - 2.0*x/5.0 + x*x/18.0;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    const double t = x*x/8.0 - 1.0;
    gsl_sf_result c;
    gsl_sf_cheb_eval_impl(&adeb4_cs, t, &c);
    result->val = c.val - 2.0*x/5.0;
    result->err = c.err + GSL_DBL_EPSILON * 2.0*x/5.0;
    return GSL_SUCCESS;
  }
  else if(x < -(M_LN2 + GSL_LOG_DBL_EPSILON)) {
    const int nexp = floor(xcut/x);
    const double ex  = exp(-x);
    double xk  = nexp * x;
    double rk  = nexp;
    double sum = 0.0;
    int i;
    for(i=nexp; i>=1; i--) {
      double xk_inv = 1.0/xk;
      sum *= ex;
      sum += ((((24.0*xk_inv + 24.0)*xk_inv + 12.0)*xk_inv + 4.0)*xk_inv + 1.0) / rk;
      rk -= 1.0;
      xk -= x;
    }
    result->val = val_infinity/(x*x*x*x) - 4.0 * sum * ex;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x < xcut) {
    const double x2 = x*x;
    const double x4 = x2*x2;
    const double sum = 24.0 + 24.0*x + 12.0*x2 + 4.0*x2*x + x4;
    result->val = (val_infinity - 4.0 * sum * exp(-x)) / x4;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    result->val = (((val_infinity/x)/x)/x)/x;
    result->err = GSL_DBL_EPSILON * result->val;
    if(result->val == 0.0)
      return GSL_EUNDRFLW;
    else
      return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_debye_1_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_debye_1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_debye_1_e", status);
  }
  return status;
}

int
gsl_sf_debye_2_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_debye_2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_debye_2_e", status);
  }
  return status;
}

int
gsl_sf_debye_3_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_debye_3_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_debye_3_e", status);
  }
  return status;
}

int
gsl_sf_debye_4_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_debye_4_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_debye_4_e", status);
  }
  return status;
}
