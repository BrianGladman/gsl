/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_transport.h"


static double transport2_data[18] = {
   1.671760446434538503,
  -0.147735359946794490,
   0.148213819946936338e-01,
  -0.14195330326305613e-02,
   0.1306541324415708e-03,
  -0.117155795867579e-04,
   0.10333498445756e-05,
  -0.901911304223e-07,
   0.78177169833e-08,
  -0.6744565684e-09,
   0.579946394e-10,
  -0.49747619e-11,
   0.425961e-12,
  -0.36422e-13,
   0.3111e-14,
  -0.265e-15,
   0.23e-16,
  -0.19e-17
};
static gsl_sf_cheb_series transport2_cs = {
  transport2_data,
  17,
  -1, 1,
  (double *)0,
  (double *)0,
  9
};

static double transport3_data[18] = {
   0.762012543243872007,
  -0.105674387705058533,
   0.119778084819657810e-01,
  -0.12144015203698307e-02,
   0.1155099769392855e-03,
  -0.105815992124423e-04,
   0.9474663385302e-06,
  -0.836221212858e-07,
   0.73109099278e-08,
  -0.6350594779e-09,
   0.549118282e-10,
  -0.47321395e-11,
   0.4067695e-12,
  -0.348971e-13,
   0.29892e-14,
  -0.256e-15,
   0.219e-16,
  -0.19e-17
};
static gsl_sf_cheb_series transport3_cs = {
  transport3_data,
  17,
  -1, 1,
  (double *)0,
  (double *)0,
  9
};


static double transport4_data[18] = {
  0.4807570994615110579,
 -0.8175378810321083956e-01,
  0.1002700665975162973e-01,
 -0.10599339359820151e-02,
  0.1034506245030405e-03,
 -0.96442705485899e-05,
  0.8745544408515e-06,
 -0.779321207981e-07,
  0.68649886141e-08,
 -0.5999571076e-09,
  0.521366241e-10,
 -0.45118382e-11,
  0.3892159e-12,
 -0.334936e-13,
  0.28767e-14,
 -0.2467e-15,
  0.211e-16,
 -0.18e-17
};
static gsl_sf_cheb_series transport4_cs = {
  transport4_data,
  17,
  -1, 1,
  (double *)0,
  (double *)0,
  9
};


static double transport5_data[18] = {
   0.347777777133910789,
  -0.66456988976050428e-01,
   0.8611072656883309e-02,
  -0.9396682223755538e-03,
   0.936324806081513e-04,
  -0.88571319340833e-05,
   0.811914989145e-06,
  -0.72957654233e-07,
   0.646971455e-08,
  -0.568490283e-09,
   0.49625598e-10,
  -0.4310940e-11,
   0.373100e-12,
  -0.32198e-13,
   0.2772e-14,
  -0.238e-15,
   0.21e-16,
  -0.18e-17
};
static gsl_sf_cheb_series transport5_cs = {
  transport5_data,
  17,
  -1, 1,
  (double *)0,
  (double *)0,
  9
};


static
double
transport_sumexp(const int numexp, const int order, const double t, double x)
{
  double rk = (double)numexp;
  double sumexp = 0.0;
  int k;
  for(k=1; k<=numexp; k++) {
    double sum2 = 1.0;
    double xk  = 1.0/(rk*x);
    double xk1 = 1.0;
    int j;
    for(j=1; j<=order; j++) {
      sum2 = sum2*xk1*xk + 1.0;
      xk1 += 1.0;
    }
    sumexp *= t;
    sumexp += sum2;
    rk -= 1.0;
  }
  return sumexp;
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_transport_2_impl(const double x, double * result)
{
  const double val_infinity = 3.289868133696452873;

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 3.0*GSL_SQRT_MACH_EPS) {
    *result = x;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double t = (x*x/8.0 - 0.5) - 0.5;
    *result = x * gsl_sf_cheb_eval(&transport2_cs, t);
    return GSL_SUCCESS;
  }
  else if(x < -GSL_LOG_MACH_EPS) {
    int    numexp = (int)((-GSL_LOG_MACH_EPS)/x) + 1;
    double sumexp = transport_sumexp(numexp, 2, exp(-x), x);
    double t = 2.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
  else if(x < 2.0/GSL_MACH_EPS) {
    int    numexp = 1;
    double sumexp = transport_sumexp(numexp, 2, 1.0, x);
    double t = 2.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
  else {
    double t = 2.0 * log(x) - x;
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
}


int
gsl_sf_transport_3_impl(const double x, double * result)
{ 
  const double val_infinity = 7.212341418957565712;

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 3.0*GSL_SQRT_MACH_EPS) {
    *result = 0.5*x*x;
    if(*result == 0.0)
      return GSL_EUNDRFLW;
    else
      return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double x2 = x*x;
    double t = (x2/8.0 - 0.5) - 0.5;
    *result = x2 * gsl_sf_cheb_eval(&transport3_cs, t);
    return GSL_SUCCESS;
  }
  else if(x < -GSL_LOG_MACH_EPS) {
    int    numexp = (int)((-GSL_LOG_MACH_EPS)/x) + 1;
    double sumexp = transport_sumexp(numexp, 3, exp(-x), x);
    double t = 3.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
  else if(x < 3.0/GSL_MACH_EPS) {
    int    numexp = 1;
    double sumexp = transport_sumexp(numexp, 3, 1.0, x);
    double t = 3.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
  else {
    double t = 3.0 * log(x) - x;
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
}


int
gsl_sf_transport_4_impl(const double x, double * result)
{
  const double val_infinity = 25.97575760906731660;

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 3.0*GSL_SQRT_MACH_EPS) {
    *result = x*x*x/3.0;
    if(*result == 0.0)
      return GSL_EUNDRFLW;
    else
      return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double x2 = x*x;
    double t = (x2/8.0 - 0.5) - 0.5;
    *result = x2*x * gsl_sf_cheb_eval(&transport4_cs, t);
    return GSL_SUCCESS;
  }
  else if(x < -GSL_LOG_MACH_EPS) {
    int    numexp = (int)((-GSL_LOG_MACH_EPS)/x) + 1;
    double sumexp = transport_sumexp(numexp, 4, exp(-x), x);
    double t = 4.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
  else if(x < 3.0/GSL_MACH_EPS) {
    int    numexp = 1;
    double sumexp = transport_sumexp(numexp, 4, 1.0, x);
    double t = 4.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
  else {
    double t = 4.0 * log(x) - x;
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
}


int
gsl_sf_transport_5_impl(const double x, double * result)
{
  const double val_infinity = 124.4313306172043912;

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 3.0*GSL_SQRT_MACH_EPS) {
    *result = x*x*x*x/4.0;
    if(*result == 0.0)
      return GSL_EUNDRFLW;
    else
      return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double x2 = x*x;
    double t = (x2/8.0 - 0.5) - 0.5;
    *result = x2*x2 * gsl_sf_cheb_eval(&transport5_cs, t);
    return GSL_SUCCESS;
  }
  else if(x < -GSL_LOG_MACH_EPS) {
    int    numexp = (int)((-GSL_LOG_MACH_EPS)/x) + 1;
    double sumexp = transport_sumexp(numexp, 5, exp(-x), x);
    double t = 5.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
  else if(x < 3.0/GSL_MACH_EPS) {
    int    numexp = 1;
    double sumexp = transport_sumexp(numexp, 5, 1.0, x);
    double t = 5.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
  else {
    double t = 5.0 * log(x) - x;
    if(t < GSL_LOG_MACH_EPS) {
      *result = val_infinity;
    }
    else {
      *result = val_infinity - exp(t);
    }
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_transport_2_e(const double x, double * result)
{
  int status = gsl_sf_transport_2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_transport_2_e", status);
  }
  return status;
}

int gsl_sf_transport_3_e(const double x, double * result)
{
  int status = gsl_sf_transport_3_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_transport_3_e", status);
  }
  return status;
}

int gsl_sf_transport_4_e(const double x, double * result)
{
  int status = gsl_sf_transport_4_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_transport_4_e", status);
  }
  return status;
}

int gsl_sf_transport_5_e(const double x, double * result)
{
  int status = gsl_sf_transport_5_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_transport_5_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*/

double gsl_sf_transport_2(const double x)
{
  double y;
  int status = gsl_sf_transport_2_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_transport_2", status);
  }
  return y;
}

double gsl_sf_transport_3(const double x)
{
  double y;
  int status = gsl_sf_transport_3_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_transport_3", status);
  }
  return y;
}

double gsl_sf_transport_4(const double x)
{
  double y;
  int status = gsl_sf_transport_4_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_transport_4", status);
  }
  return y;
}

double gsl_sf_transport_5(const double x)
{
  double y;
  int status = gsl_sf_transport_5_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_transport_5", status);
  }
  return y;
}
