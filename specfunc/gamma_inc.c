/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_errno.h>
#include <gsl_math.h>
#include "gsl_sf_erf.h"
#include "gsl_sf_log.h"
#include "gsl_sf_gamma.h"


/* The dominant part,
 * D(a,x) := x^a e^(-x) / Gamma(a+1)
 */
static
int
gamma_inc_D(const double a, const double x, double * result)
{
  if(a < 10.0) {
    double lnr;
    double lg;
    gsl_sf_lngamma_impl(a+1.0, &lg);
    lnr = a * log(x) - x - lg;
    *result = exp(lnr);
    return GSL_SUCCESS;
  }
  else {
    double mu = (x-a)/a;
    double ln_term;
    double gstar;
    gsl_sf_log_1plusx_mx_impl(mu, &ln_term);  /* log(1+mu) - mu */
    gsl_sf_gammastar_impl(a, &gstar);
    *result = exp(a*ln_term)/sqrt(2.0*M_PI*a)/gstar;
    return GSL_SUCCESS;
  }
}


/* P series representation.
 */
static
int
gamma_inc_P_series(const double a, const double x, double * result)
{
  const int nmax = 5000;

  double D;
  int stat_D = gamma_inc_D(a, x, &D);
  
  double sum  = 1.0;
  double term = 1.0;
  int n;
  for(n=1; n<nmax; n++) {
    term *= x/(a+n);
    sum  += term;
    if(fabs(term/sum) < GSL_MACH_EPS) break;
  }

  *result = D * sum;

  if(n == nmax)
    return GSL_EMAXITER;
  else
    return stat_D;
}


/* Q large x asymptotic
 */
static
int
gamma_inc_Q_large_x(const double a, const double x, double * result)
{
  const int nmax = 5000;

  double D;
  const int stat_D = gamma_inc_D(a, x, &D);
  const double pre = D * (a/x);

  double sum  = 1.0;
  double term = 1.0;
  double last = 1.0;
  int n;
  for(n=1; n<nmax; n++) {
    term *= (a-n)/x;
    if(fabs(last/term) > 1.0) break;
    if(fabs(term/sum) < GSL_MACH_EPS) break;
    sum  += term;
    last  = term;
  }

  *result = pre * sum;

  if(n == nmax)
    return GSL_EMAXITER;
  else
    return stat_D;
}


/* Uniform asymptotic for x near a, a and x large.
 * See [Temme, p. 285]
 */
static
int
gamma_inc_Q_asymp_unif(const double a, const double x, double * result)
{
  const double rta = sqrt(a);
  const double eps = (x-a)/a;

  double ln_term;
  const int stat_ln = gsl_sf_log_1plusx_mx_impl(eps, &ln_term);  /* log(1+eps) - eps */
  const double eta  = eps * sqrt(-2.0*ln_term/(eps*eps));

  double erfc = gsl_sf_erfc(eta*M_SQRT2*rta);

  double R;
  double c0, c1;
  
  if(fabs(eps) < GSL_ROOT5_MACH_EPS) {
    c0 = -1.0/3.0 + eps*(1.0/12.0 - eps*(23.0/540.0 - eps*(353.0/12960.0 - eps*589.0/30240.0)));
    c1 = 0.0;
  }
  else {
    double rt_term;
    rt_term = sqrt(-2.0 * ln_term/(eps*eps));
    c0 = (1.0 - 1.0/rt_term)/eps;
    c1 = 0.0;
  }

  R = exp(-0.5*a*eta*eta)/(M_SQRT2*M_SQRTPI*rta) * (c0 + c1/a);

  *result = 0.5 * erfc + R;
  return stat_ln;
}


int
gsl_sf_gamma_inc_Q_impl(const double a, const double x, double * result)
{
  if(a <= 0.0 || x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(a > 1.0e+05 && (x-a)*(x-a) < a) {
    return gamma_inc_Q_asymp_unif(a, x, result);
  }
  else if(x < a) {
    double P;
    int stat_P = gamma_inc_P_series(a, x, &P);
    *result = 1.0 - P;
    return stat_P;
  }
  else {
    return gamma_inc_Q_large_x(a, x, result);
  }
}


int
gsl_sf_gamma_inc_P_impl(const double a, const double x, double * result)
{
  if(a <= 0.0 || x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(a > 1.0e+05 && (x-a)*(x-a) < a) {
    double Q;
    int stat_Q = gamma_inc_Q_asymp_unif(a, x, &Q);
    *result = 1.0 - Q;
    return stat_Q;
  }
  else if(x < a) {
    return gamma_inc_P_series(a, x, result);
  }
  else {
    double Q;
    int stat_Q = gamma_inc_Q_large_x(a, x, &Q);
    *result = 1.0 - Q;
    return stat_Q;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/


int
gsl_sf_gamma_inc_P_e(const double a, const double x, double * result)
{
  int status = gsl_sf_gamma_inc_P_impl(a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_gamma_inc_P_e", status);
  }
  return status;
}

int
gsl_sf_gamma_inc_Q_e(const double a, const double x, double * result)
{
  int status = gsl_sf_gamma_inc_Q_impl(a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_gamma_inc_Q_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double
gsl_sf_gamma_inc_P(const double a, const double x)
{
  double y;
  int status = gsl_sf_gamma_inc_P_impl(a, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_gamma_inc_P", status);
  }
  return y;
}

double
gsl_sf_gamma_inc_Q(const double a, const double x)
{
  double y;
  int status = gsl_sf_gamma_inc_Q_impl(a, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_gamma_inc_Q", status);
  }
  return y;
}
