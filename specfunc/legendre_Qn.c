/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"
#include "gsl_sf_elementary.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"


/* Evaluate f_{ell+1}/f_ell
 * f_ell := Q^{b}_{a+ell}(x)
 * x > 1
 */
static
int
legendreQ_CF1_xgt1(int ell, double a, double b, double x, double * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int n = 1;
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = ell + 1.0 + a + b;
  double b1 = (2.0*(ell+1.0+a) + 1.0) * x;
  double An = b1*Anm1 + a1*Anm2;
  double Bn = b1*Bnm1 + a1*Bnm2;
  double an, bn;
  double fn = An/Bn;

  while(n < maxiter) {
    double old_fn;
    double del;
    double lna;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    lna = ell + n + a;
    an = b*b - lna*lna;
    bn = (2.0*lna + 1.0) * x;
    An = bn*Anm1 + an*Anm2;
    Bn = bn*Bnm1 + an*Bnm2;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An /= RECUR_BIG;
      Bn /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
    }

    old_fn = fn;
    fn = An/Bn;
    del = old_fn/fn;

    if(fabs(del - 1.0) < 10.0*GSL_MACH_EPS) break;
  }

  *result = fn;

  if(n == maxiter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS; 
}


/* Uniform asymptotic for Q_l(x).
 * Assumes x > -1.0 and x != 1.0.
 * Discards second order and higher terms.
 */
static
int
legendre_Ql_asymp_unif(const double ell, const double x, double * result)
{
  if(x < 1.0) {
    double u   = ell + 0.5;
    double th  = acos(x);
    double Y0, Y1;
    int stat_Y0, stat_Y1;
    int stat_m;
    double pre;
    double B00;
    double sum;

    /* B00 = 1/8 (1 - th cot(th) / th^2
     * pre = sqrt(th/sin(th))
     */
    if(th < GSL_ROOT4_MACH_EPS) {
      B00 = (1.0 + th*th/15.0)/24.0;
      pre = 1.0 + th*th/12.0;
    }
    else {
      double sin_th = sqrt(1.0 - x*x);
      double cot_th = x / sin_th;
      B00 = 1.0/8.0 * (1.0 - th * cot_th) / (th*th);
      pre = sqrt(th/sin_th);
    }

    stat_Y0 = gsl_sf_bessel_Y0_impl(u*th, &Y0);
    stat_Y1 = gsl_sf_bessel_Y1_impl(u*th, &Y1);

    sum = -0.5*M_PI * (Y0 + th/u * Y1 * B00);

    stat_m = gsl_sf_multiply_impl(pre, sum, result);

    return GSL_ERROR_SELECT_3(stat_m, stat_Y0, stat_Y1);
  }
  else {
    double u   = ell + 0.5;
    double xi  = acosh(x);
    double K0_scaled, K1_scaled;
    int stat_K0, stat_K1;
    int stat_e;
    double pre;
    double B00;
    double sum;

    /* B00 = -1/8 (1 - xi coth(xi) / xi^2
     * pre = sqrt(xi/sinh(xi))
     */
    if(xi < GSL_ROOT4_MACH_EPS) {
      B00 = (1.0-xi*xi/15.0)/24.0;
      pre = 1.0 - xi*xi/12.0;
    }
    else {
      double sinh_xi = sqrt(x*x - 1.0);
      double coth_xi = x / sinh_xi;
      B00 = -1.0/8.0 * (1.0 - xi * coth_xi) / (xi*xi);
      pre = sqrt(xi/sinh_xi);
    }

    stat_K0 = gsl_sf_bessel_K0_scaled_impl(u*xi, &K0_scaled);
    stat_K1 = gsl_sf_bessel_K1_scaled_impl(u*xi, &K1_scaled);

    sum = K0_scaled - xi/u * K1_scaled * B00;

    stat_e = gsl_sf_exp_mult_impl(-u*xi, pre * sum, result);

    return GSL_ERROR_SELECT_3(stat_e, stat_K0, stat_K1);
  }
}



/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_legendre_Q0_impl(const double x, double * result)
{
  if(x <= -1.0 || x == 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0){
    *result = 0.5 * log((1.0+x)/(1.0-x));
    return GSL_SUCCESS;
  }
  else {
    *result = 0.5 * log((x+1.0)/(x-1.0));
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Q1_impl(const double x, double * result)
{
  if(x <= -1.0 || x == 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0){
    *result = 0.5 * x * log((1.0+x)/(1.0-x)) - 1.0;
    return GSL_SUCCESS;
  }
  else {
    *result = 0.5 * x * log((x+1.0)/(x-1.0)) - 1.0;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Q2_impl(const double x, double * result)
{
  if(x <= -1.0 || x == 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0){
    *result = 0.25 * (3.0*x*x-1.0) * log((1.0+x)/(1.0-x)) - 1.5*x;
    return GSL_SUCCESS;
  }
  else {
    *result = 0.25 * (3.0*x*x-1.0) * log((x+1.0)/(x-1.0)) - 1.5*x;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Ql_impl(const int l, const double x, double * result)
{
  if(x <= -1.0 || x == 1.0 || l < 0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(l == 0) {
    return gsl_sf_legendre_Q0_impl(x, result);
  }
  else if(l == 1) {
    return gsl_sf_legendre_Q1_impl(x, result);
  }
  else if(l > 100000) {
    return legendre_Ql_asymp_unif(l, x, result);
  }
  else if(x < 1.0){
    /* Forward recurrence.
     */
    double Q0, Q1;
    int stat_Q0 = gsl_sf_legendre_Q0_impl(x, &Q0);
    int stat_Q1 = gsl_sf_legendre_Q1_impl(x, &Q1);
    double Qellm1 = Q0;
    double Qell   = Q1;
    double Qellp1;
    int ell;
    for(ell=1; ell<l; ell++) {
      Qellp1 = (x*(2.0*ell + 1.0) * Qell - ell * Qellm1) / (ell + 1.0);
      Qellm1 = Qell;
      Qell   = Qellp1;
    }
    *result = Qell;
    return GSL_ERROR_SELECT_2(stat_Q0, stat_Q1);
  }
  else {
    /* x > 1.0 */

    double rat;
    int stat_CF1  = legendreQ_CF1_xgt1(l, 0.0, 0.0, x, &rat);
    int stat_Q;
    double Qellp1 = rat * GSL_SQRT_DBL_MIN;
    double Qell   = GSL_SQRT_DBL_MIN;
    double Qellm1;
    int ell;
    for(ell=l; ell>0; ell--) {
      Qellm1 = (x * (2.0*ell + 1.0) * Qell - (ell+1.0) * Qellp1) / ell;
      Qellp1 = Qell;
      Qell   = Qellm1;
    }

    if(fabs(Qell) > fabs(Qellp1)) {
      double Q0;
      stat_Q = gsl_sf_legendre_Q0_impl(x, &Q0);
      *result = GSL_SQRT_DBL_MIN * Q0 / Qell;
    }
    else {
      double Q1;
      stat_Q = gsl_sf_legendre_Q1_impl(x, &Q1);
      *result = GSL_SQRT_DBL_MIN * Q1 / Qellp1;
    }

    return GSL_ERROR_SELECT_2(stat_Q, stat_CF1);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_legendre_Q0_e(const double x, double * result)
{
  int status = gsl_sf_legendre_Q0_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Q0_e", status);
  }
  return status;
}


int
gsl_sf_legendre_Q1_e(const double x, double * result)
{
  int status = gsl_sf_legendre_Q1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Q1_e", status);
  }
  return status;
}


int
gsl_sf_legendre_Q2_e(const double x, double * result)
{
  int status = gsl_sf_legendre_Q2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Q2_e", status);
  }
  return status;
}


int
gsl_sf_legendre_Ql_e(const int l, const double x, double * result)
{
  int status = gsl_sf_legendre_Ql_impl(l, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Ql_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double
gsl_sf_legendre_Q0(const double x)
{
  double y;
  int status = gsl_sf_legendre_Q0_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Q0", status);
  }
  return y;
}


double
gsl_sf_legendre_Q1(const double x)
{
  double y;
  int status = gsl_sf_legendre_Q1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Q1", status);
  }
  return y;
}


double
gsl_sf_legendre_Q2(const double x)
{
  double y;
  int status = gsl_sf_legendre_Q2_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Q2", status);
  }
  return y;
}


double
gsl_sf_legendre_Ql(const int l, const double x)
{
  double y;
  int status = gsl_sf_legendre_Ql_impl(l, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Ql", status);
  }
  return y;
}
