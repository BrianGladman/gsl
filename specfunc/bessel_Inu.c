/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_temme.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#define locMin(a,b)  ((a) < (b) ? (a) : (b))
#define locMax(a,b)  ((a) > (b) ? (a) : (b))


#if 0
/* Evaluate continued fraction CF1 to obtain I'/I.
 */
static
int
bessel_I_CF1(const double nu, const double x, double * result)
{
  const int max_iter = 5000;
  int i = 0;
  double x_inv = 1.0/x;
  double r = nu * x_inv;
  double b = 2.0*nu*x_inv;
  double d = 0.0;
  double c;
  r = locMax(r, 1.0e-100);
  c = r;
  while(i < max_iter) {
    double del;
    b  += 2.0*x_inv;
    d   = 1.0/(b+d);
    c   = b + 1.0/c;
    del = c*d;
    r  *= del;
    if (fabs(del-1.0) < GSL_MACH_EPS) break;
    ++i;
  }
  
  *result = r;
  if(i == max_iter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}
#endif /* 0 */


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Inu_scaled_impl(double nu, double x, double * result)
{
  if(x < 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else if(x*x < 10.0*(nu+1.0)*GSL_ROOT5_MACH_EPS) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 4, &b);
    *result = ex * b;
    return stat;
  }
  else if(x*x < 10.0*(nu+1.0)) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 100, &b);
    *result = ex * b;
    return stat;
  }
  else if(0.5/(nu*nu + x*x) < GSL_ROOT3_MACH_EPS) {
    return gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nu, x, result);
  }
  else {
    int N = floor(nu);
    int M = ceil(sqrt(0.5/GSL_ROOT3_MACH_EPS - x*x + 1.0)) + 1;
    double nu_frac = nu - N;
    double mu = M + nu_frac;
    double I_mu, Ip_mu, I_mup1;
    double I_nu, Ip_nu;
    gsl_sf_bessel_Inu_scaled_asymp_unif_impl(mu,     x, &I_mu);
    gsl_sf_bessel_Inu_scaled_asymp_unif_impl(mu+1.0, x, &I_mup1);
    Ip_mu = mu/x * I_mu + I_mup1;
    gsl_sf_bessel_I_recur(nu, x, M-N, I_mu, Ip_mu, &I_nu, &Ip_nu, (double*)0, (double *)0);
    *result = I_nu;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_bessel_Inu_impl(double nu, double x, double * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    *result = 0.0;
    return GSL_EOVRFLW;
  }
  else {
    double ex = exp(x);
    double b = 0.0;
    int stat_I = gsl_sf_bessel_Inu_scaled_impl(nu, x, &b);
    *result = ex * b;
    return stat_I;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Inu_scaled_e(double nu, double x, double * result)
{
  int status = gsl_sf_bessel_Inu_scaled_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Inu_scaled_e", status);
  }
  return status;
}


int
gsl_sf_bessel_Inu_e(double nu, double x, double * result)
{
  int status = gsl_sf_bessel_Inu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Inu_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*/

double
gsl_sf_bessel_Inu_scaled(double nu, double x)
{
  double y;
  int status = gsl_sf_bessel_Inu_scaled_impl(nu, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Inu_scaled", status);
  }
  return y;
}


double
gsl_sf_bessel_Inu(double nu, double x)
{
  double y;
  int status = gsl_sf_bessel_Inu_impl(nu, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Inu", status);
  }
  return y;
}
