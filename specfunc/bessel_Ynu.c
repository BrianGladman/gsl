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


/* Perform forward recurrence for Y_nu(x) and Y'_nu(x)
 *
 *        Y_{nu+1} =  nu/x Y_nu - Y'_nu
 *       Y'_{nu+1} = -(nu+1)/x Y_{nu+1} + Y_nu
 */
static
int
bessel_Y_recur(const double nu_min, const double x, const int kmax,
               const double Y_start, const double Yp_start,
	       double * Y_end, double * Yp_end)
{
  double x_inv = 1.0/x;
  double nu = nu_min;
  double Y_nu  = Y_start;
  double Yp_nu = Yp_start;
  int k;

  for(k=1; k<=kmax; k++) {
    double nuox = nu*x_inv;
    double Y_nu_save = Y_nu;
    Y_nu  = -Yp_nu + nuox * Y_nu;
    Yp_nu = Y_nu_save - (nuox+x_inv) * Y_nu;
    nu += 1.0;
  }
  *Y_end  = Y_nu;
  *Yp_end = Yp_nu;
  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Ynu_impl(double nu, double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= 0.0 || nu < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(nu > 50.0) {
    return gsl_sf_bessel_Ynu_asymp_Olver_impl(nu, x, result);
  }
  else {
    /* -1/2 <= mu <= 1/2 */
    int N = (int)(nu + 0.5);
    double mu = nu - N;
 
    double Ymu, Ymup1;
    int stat_mu;
    double Ynm1;
    double Yn;
    double Ynp1;
    int n;

    if(x < 2.0) {
      /* Determine Ymu, Ymup1 directly. This is really
       * an optimization since this case could as well
       * be handled by a call to gsl_sf_bessel_JY_mu_restricted(),
       * as below.
       */
      gsl_sf_result Y_mu, Y_mup1;
      stat_mu = gsl_sf_bessel_Y_temme(mu, x, &Y_mu, &Y_mup1);
      Ymu   = Y_mu.val;
      Ymup1 = Y_mup1.val;
    }
    else {
      /* Determine Ymu, Ymup1 and Jmu, Jmup1.
       */
      gsl_sf_result Y_mu, Y_mup1;
      gsl_sf_result J_mu, J_mup1;
      stat_mu = gsl_sf_bessel_JY_mu_restricted(mu, x, &J_mu, &J_mup1, &Y_mu, &Y_mup1);
      Ymu   = Y_mu.val;
      Ymup1 = Y_mup1.val;
    }

    /* Forward recursion to get Ynu, Ynup1.
     */
    Ynm1 = Ymu;
    Yn   = Ymup1;
    for(n=1; n<=N; n++) {
      Ynp1 = 2.0*(mu+n)/x * Yn - Ynm1;
      Ynm1 = Yn;
      Yn   = Ynp1;
    }

    result->val = Ynm1; /* Y_nu */
    result->err = GSL_DBL_EPSILON * (0.5*N + 2.0) * fabs(Ynm1);

    return stat_mu;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Ynu_e(const double nu, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_Ynu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Ynu_e", status);
  }
  return status;
}
