/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_olver.h"
#include "bessel_temme.h"
#include "gsl_sf_bessel.h"


#if 0
static
int
bessel_J_CF1(const double nu, const double x, double * result, double * sign)
{
  const int max_iter = 5000;
  const double SMALL = 1.0e-100;

  int i = 1;
  int isign = 1;
  double x_inv = 1.0/x;
  double r = nu * x_inv;
  double b = 2.0*nu*x_inv;
  double d = 0.0;
  double c;
  r = GSL_MAX(r, SMALL);
  c = r;
  while(i < max_iter) {
    double del;
    b  += 2.0*x_inv;
    d   = b - d;
    if(fabs(d) < SMALL) d = SMALL;
    c   = b - 1.0/c;
    if(fabs(c) < SMALL) c = SMALL;
    d = 1.0/d;
    del = c*d;
    r *= del;
    if(d < 0.0) isign = -isign;
    if(fabs(del-1.0) < 0.5*GSL_DBL_EPSILON) break;
  }
  
  *result = r;
  *sign = isign;
  if(i == max_iter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}
#endif /* 0 */


/* Evaluate at large enough nu to apply asymptotic
 * results and apply backward recurrence.
 */
static
int
bessel_J_recur_asymp(const double nu, const double x, double * Jnu, double * Jnup1)
{
  const double nu_cut = 25.0;
  int n;
  int steps = ceil(nu_cut - nu) + 1;

  gsl_sf_result r_Jnp1;
  gsl_sf_result r_Jn;
  int stat_O1 = gsl_sf_bessel_Jnu_asymp_Olver_impl(nu + steps + 1.0, x, &r_Jnp1);
  int stat_O2 = gsl_sf_bessel_Jnu_asymp_Olver_impl(nu + steps,       x, &r_Jn);
  double Jnp1 = r_Jnp1.val;
  double Jn   = r_Jn.val;
  double Jnm1;
  double Jnp1_save;

  for(n=steps; n>0; n--) {
    Jnm1 = 2.0*(nu+n)/x * Jn - Jnp1;
    Jnp1 = Jn;
    Jnp1_save = Jn;
    Jn   = Jnm1;
  }

  *Jnu   = Jn;
  *Jnup1 = Jnp1_save;
  return GSL_ERROR_SELECT_2(stat_O1, stat_O2);
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Jnu_impl(const double nu, const double x, gsl_sf_result * result)
{
  const double nu_cut = 25.0;

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x < 0.0 || nu < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    double b;
    int stat = gsl_sf_bessel_JnuYnu_zero(nu, &b, (double *)0, (double *)0, (double *)0);
    result->val = b;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
  else if(x*x < 10.0*(nu+1.0)*GSL_ROOT5_DBL_EPSILON) {
    double b;
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, -1, 50, GSL_DBL_EPSILON, &b);
    result->val = b;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
  else if(x*x < 10.0*(nu+1.0)) {
    double b;
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, -1, 100, GSL_DBL_EPSILON, &b);
    result->val = b;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
  else if(nu > nu_cut) {
    return gsl_sf_bessel_Jnu_asymp_Olver_impl(nu, x, result);
  }
  else {
    double Jnu, Jnup1;
    int status = bessel_J_recur_asymp(nu, x, &Jnu, &Jnup1);
    result->val = Jnu;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(0.5 * GSL_MAX_DBL(nu,1.0) * Jnu);
    return status;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Jnu_e(const double nu, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_Jnu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Jnu_e", status);
  }
  return status;
}
