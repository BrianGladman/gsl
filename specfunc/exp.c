/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
#include "gsl_sf_exp.h"


/* Evaluate the continued fraction for exprel.
 * [Abramowitz+Stegun, 4.2.41]
 */
static
int
exprel_n_CF(const int N, const double x, double * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int n = 1;
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = 1.0;
  double b1 = 1.0;
  double a2 = -x;
  double b2 = N+1;
  double an, bn;

  double fn;

  double An = b1*Anm1 + a1*Anm2;   /* A1 */
  double Bn = b1*Bnm1 + a1*Bnm2;   /* B1 */
  
  /* One explicit step, before we get to the main pattern. */
  n++;
  Anm2 = Anm1;
  Bnm2 = Bnm1;
  Anm1 = An;
  Bnm1 = Bn;
  An = b2*Anm1 + a2*Anm2;   /* A2 */
  Bn = b2*Bnm1 + a2*Bnm2;   /* B2 */

  fn = An/Bn;

  while(n < maxiter) {
    double old_fn;
    double del;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    an = ( GSL_IS_ODD(n) ? ((n-1)/2)*x : -(N+(n/2)-1)*x );
    bn = N + n - 1;
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


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_exp_impl(const double x, double * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    *result = 0.0;
    return GSL_EOVRFLW;
  }
  else if(x < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    *result = exp(x);
    return GSL_SUCCESS;
  }
}

int gsl_sf_exp_sgn_impl(const double x, const double sgn, double * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    *result = 0.0;
    return GSL_EOVRFLW;
  }
  else if(x < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    *result = GSL_SIGN(sgn) * exp(x);
    return GSL_SUCCESS;
  }
}

int gsl_sf_exp_mult_impl(const double x, const double y, gsl_sf_result * result)
{
  const double ay  = fabs(y);

  if(y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
          && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
    ) {
    double ex = exp(x);
    result->val = y * ex;
    result->err = fabs(result->val * GSL_DBL_EPSILON);
    return GSL_SUCCESS;
  }
  else {
    const double ly  = log(ay);
    const double lnr = x + ly;

    if(lnr > GSL_LOG_DBL_MAX - 0.01) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EOVRFLW;
    }
    else if(lnr < GSL_LOG_DBL_MIN + 0.01) {
      *result = 0.0;
      return GSL_EUNDRFLW;
    }
    else {
      const double sy  = GSL_SIGN(y);
      const double M   = floor(x);
      const double N   = floor(ly);
      const double a   = x  - M;
      const double b   = ly - N;
      result->val = sy * exp(M+N) * exp(a+b);
      result->err = fabs(result->val * GSL_DBL_EPSILON);
      return GSL_SUCCESS;
    }
  }
}

int gsl_sf_expm1_impl(const double x, double * result)
{
  const double cut = 0.002;

  if(x < GSL_LOG_DBL_MIN) {
    *result = -1.0;
    return GSL_SUCCESS;
  }
  else if(x < -cut) {
    *result = exp(x) - 1.0;
    return GSL_SUCCESS;
  }
  else if(x < cut) {
    *result = x * (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
    return GSL_SUCCESS;
  } 
  else if(x < GSL_LOG_DBL_MAX) {
    *result = exp(x) - 1.0;
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
}


int gsl_sf_exprel_impl(const double x, double * result)
{
  const double cut = 0.002;

  if(x < GSL_LOG_DBL_MIN) {
    *result = -1.0/x;
    return GSL_SUCCESS;
  }
  else if(x < -cut) {
    *result = (exp(x) - 1.0)/x;
    return GSL_SUCCESS;
  }
  else if(x < cut) {
    *result = (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
    return GSL_SUCCESS;
  } 
  else if(x < GSL_LOG_DBL_MAX) {
    *result = (exp(x) - 1.0)/x;
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
}


int gsl_sf_exprel_2_impl(double x, double * result)
{
  const double cut = 0.002;

  if(x < GSL_LOG_DBL_MIN) {
    *result = -2.0/x*(1.0 + 1.0/x);
    return GSL_SUCCESS;
  }
  else if(x < -cut) {
    *result = 2.0*(exp(x) - 1.0 - x)/(x*x);
    return GSL_SUCCESS;
  }
  else if(x < cut) {
    *result = (1.0 + 1.0/3.0*x*(1.0 + 0.25*x*(1.0 + 0.2*x*(1.0 + 1.0/6.0*x))));
    return GSL_SUCCESS;
  } 
  else if(x < GSL_LOG_DBL_MAX) {
    *result = 2.0*(exp(x) - 1.0 - x)/(x*x);
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
}


int
gsl_sf_exprel_n_impl(const int N, const double x, double * result)
{
  if(N < 0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(fabs(x) < GSL_ROOT3_MACH_EPS * N) {
    *result = 1.0 + x/(N+1) * (1.0 + x/(N+2));
    return GSL_SUCCESS;
  }
  else if(N == 0) {
    return gsl_sf_exp_impl(x, result);
  }
  else if(N == 1) {
    return gsl_sf_exprel_impl(x, result);
  }
  else if(N == 2) {
    return gsl_sf_exprel_2_impl(x, result);
  }
  else {
    if(x > N && (-x + N*(1.0 + log(x/N)) < GSL_LOG_MACH_EPS)) {
      /* x is much larger than n.
       * Ignore polynomial part, so
       * exprel_N(x) ~= e^x N!/x^N
       */
      double lnf_N;
      double lnr;
      gsl_sf_lnfact_impl(N, &lnf_N);
      lnr = x + lnf_N - N*log(x);
      return gsl_sf_exp_impl(lnr, result);
    }
    else if(x > N) {
      /* Write the identity
       *   exprel_n(x) = e^x n! / x^n (1 - Gamma[n,x]/Gamma[n])
       * then use the asymptotic expansion
       * Gamma[n,x] ~ x^(n-1) e^(-x) (1 + (n-1)/x + (n-1)(n-2)/x^2 + ...)
       */
      double ln_x = log(x);
      double lnf_N, lg_N;
      double lnpre;
      gsl_sf_lnfact_impl(N, &lnf_N);    /* log(N!)       */
      lg_N  = lnf_N - log(N);           /* log(Gamma(N)) */
      lnpre = x + lnf_N - N*ln_x;
      if(lnpre < GSL_LOG_DBL_MAX - 5.0) {
        int stat_eG;
	double pre = exp(lnpre);
	double bigG_ratio;
        double ln_bigG_ratio_pre = -x + (N-1)*ln_x - lg_N;
	double bigGsum = 1.0;
	double term = 1.0;
	int k;
	for(k=1; k<N; k++) {
	  term *= (N-k)/x;
	  bigGsum += term;
	}
	stat_eG = gsl_sf_exp_mult_impl(ln_bigG_ratio_pre, bigGsum, &bigG_ratio);
	if(stat_eG == GSL_SUCCESS) {
          *result = pre * (1.0 - bigG_ratio);
	  return GSL_SUCCESS;
	}
	else {
	  *result = 0.0;
	  return stat_eG;
	}
      }
      else {
        *result = 0.0;
	return GSL_EOVRFLW;
      }
    }
    else if(x > -10.0*N) {
      return exprel_n_CF(N, x, result);
    }
    else {
      /* x -> -Inf asymptotic:
       * exprel_n(x) ~ e^x n!/x^n - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
       *             ~ - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
       */
      double sum  = 1.0;
      double term = 1.0;
      int k;
      for(k=1; k<N; k++) {
        term *= (N-k)/x;
	sum  += term;
      }
      *result = -N/x * sum;
      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_exp_err_impl(const double x, const double dx, gsl_sf_result * result)
{
  const double adx = fabs(dx);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x + adx > GSL_LOG_DBL_MAX) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EOVRFLW;
  }
  else if(x - adx < GSL_LOG_DBL_MIN) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    const double ex  = exp(x);
    const double edx = exp(adx);
    result->val = ex;
    result->err = ex * GSL_MAX_DBL(GSL_DBL_EPSILON, edx - 1.0/edx);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_exp_e(const double x, double * result)
{
  int status = gsl_sf_exp_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_exp_e", status);
  }
  return status;
}


int gsl_sf_exp_sgn_e(const double x, const double sgn, double * result)
{
  int status = gsl_sf_exp_sgn_impl(x, sgn, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_exp_sgn_e", status);
  }
  return status;
}


int gsl_sf_exp_mult_e(const double x, const double y, double * result)
{
  int status = gsl_sf_exp_mult_impl(x, y, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_exp_mult_e", status);
  }
  return status;
}


int gsl_sf_expm1_e(const double x, double * result)
{
  int status = gsl_sf_expm1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_expm1_e", status);
  }
  return status;
}


int gsl_sf_exprel_e(const double x, double * result)
{
  int status = gsl_sf_exprel_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_exprel_e", status);
  }
  return status;
}


int gsl_sf_exprel_2_e(const double x, double * result)
{
  int status = gsl_sf_exprel_2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_exprel_2_e", status);
  }
  return status;
}


int gsl_sf_exprel_n_e(const int n, const double x, double * result)
{
  int status = gsl_sf_exprel_n_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_exprel_n_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_exp(const double x)
{
  double y;
  int status = gsl_sf_exp_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_exp", status);
  }
  return y;
}


double gsl_sf_exp_sgn(const double x, const double sgn)
{
  double y;
  int status = gsl_sf_exp_sgn_impl(x, sgn, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_exp_sgn", status);
  }
  return y;
}


double gsl_sf_exp_mult(const double x, const double y)
{
  double r;
  int status = gsl_sf_exp_mult_impl(x, y, &r);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_exp_mult", status);
  }
  return r;
}


double gsl_sf_expm1(const double x)
{
  double y;
  int status = gsl_sf_expm1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_expm1", status);
  }
  return y;
}


double gsl_sf_exprel(const double x)
{
  double y;
  int status = gsl_sf_exprel_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_exprel", status);
  }
  return y;
}


double gsl_sf_exprel_2(const double x)
{
  double y;
  int status = gsl_sf_exprel_2_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_exprel_2", status);
  }
  return y;
}


double gsl_sf_exprel_n(const int n, const double x)
{
  double y;
  int status = gsl_sf_exprel_n_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_exprel_n", status);
  }
  return y;
}
