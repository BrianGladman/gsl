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



/* Evaluate at large enough nu to apply asymptotic
 * results and apply backward recurrence.
 */
static
int
bessel_J_recur_asymp(const double nu, const double x,
                     gsl_sf_result * Jnu, gsl_sf_result * Jnup1)
{
  const double nu_cut = 25.0;
  int n;
  int steps = ceil(nu_cut - nu) + 1;

  gsl_sf_result r_Jnp1;
  gsl_sf_result r_Jn;
  int stat_O1 = gsl_sf_bessel_Jnu_asymp_Olver_impl(nu + steps + 1.0, x, &r_Jnp1);
  int stat_O2 = gsl_sf_bessel_Jnu_asymp_Olver_impl(nu + steps,       x, &r_Jn);
  double r_fe = fabs(r_Jnp1.err/r_Jnp1.val) + fabs(r_Jn.err/r_Jn.val);
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

  Jnu->val = Jn;
  Jnu->err = (r_fe + GSL_DBL_EPSILON * (steps + 1.0)) * fabs(Jn);
  Jnup1->val = Jnp1_save;
  Jnup1->err = (r_fe + GSL_DBL_EPSILON * (steps + 1.0)) * fabs(Jnp1_save);

  return GSL_ERROR_SELECT_2(stat_O1, stat_O2);
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Jnu_impl(const double nu, const double x, gsl_sf_result * result)
{
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
  else if(nu > 25.0) {
    return gsl_sf_bessel_Jnu_asymp_Olver_impl(nu, x, result);
  }
  else {
    gsl_sf_result Jnup1;
    return bessel_J_recur_asymp(nu, x, result, &Jnup1);
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
