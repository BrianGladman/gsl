/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_legendre.h"

#include "legendre_impl.h"


/* [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991) */

static inline double olver_U1(double beta2, double p)
{
  return (p-1.)/(24.*(1.+beta2)) * (3. + beta2*(2. + 5.*p*(1.+p)));
}
static inline double olver_U2(double beta2, double p)
{
  double beta4 = beta2*beta2;
  double p2    = p*p;
  double poly1 =  4.*beta4 + 84.*beta2 - 63.;
  double poly2 = 16.*beta4 + 90.*beta2 - 81.;
  double poly3 = beta2*p2*(97.*beta2 - 432. + 77.*p*(beta2-6.) - 385.*beta2*p2*(1. + p));
  return (1.-p)/(1152.*(1.+beta2)) * (poly1 + poly2 + poly3);
}

/* P^{-mu}_{-1/2 + I tau}, mu -> Inf */
int gsl_sf_conical_P_xlt1_large_mu_impl(double mu, double tau, double x, double * result)
{
  double beta  = tau/mu;
  double beta2 = beta*beta;
  double S     = beta * acos((1.-beta2)/(1.+beta2));
  double p     = x/sqrt(beta2*(1.-x*x) + 1.);
  double ln_pre_1 =  0.5*mu*(S - log(1.+beta2) + log((1.-p)/(1.+p))) - gsl_sf_lngamma(mu+1.);
  double ln_pre_2 = -0.25 * log(1. + beta2*(1.-x));
  double ln_pre_3 = -tau * atan(p*beta);
  double ln_pre = ln_pre_1 + ln_pre_2 + ln_pre_3;
  
  if(ln_pre > GSL_LOG_DBL_MAX) {
    *result = 0.; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(ln_pre < GSL_LOG_DBL_MIN) {
    *result = 0.;
    return GSL_EUNDRFLW;
  }
  else {
    double sum = 1. - olver_U1(beta2, p)/mu + olver_U2(beta2, p)/(mu*mu);
    *result = exp(ln_pre) * sum;
    return GSL_SUCCESS;
  }
}


/* P^{-mu}_{-1/2 + I tau}, mu -> Inf */
int gsl_sf_conical_P_xgt1_large_mu_impl(double mu, double tau, double x, double * result)
{
}


/* A_n^{-mu}, B_n^{-mu}   [Olver, p.465, 469] */

static double olver_B0_xi(double mu, double xi)
{
  return (1. - 4.*mu*mu)/(8.*xi) * (coth(xi) - 1./xi);
}
static double olver_A1_xi(double mu, double xi, double x)
{
  double B = olver_B0_xi(mu, xi);
  double psi = (4.*mu*mu - 1.)/16. * (1./(x*x-1.) - 1./(xi*xi));
  return 0.5*xi*xi*B*B + (mu+0.5)*B - psi + mu/6.*(0.25 - mu*mu);
}


/* P^{-mu}_{-1/2 + I tau}, tau -> Inf   [Olver, p. 473]
 * 
 */
int gsl_sf_conical_P_xlt1_large_tau_impl(double mu, double tau, double x, double * result)
{
  
}
