/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "legendre.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_trig.h"
#include "gsl_sf_legendre.h"


/* See [Abbott+Schaefer, Ap.J. 308, 546 (1986)] for
 * enough details to follow what is happening here.
 */


/* Logarithm of normalization factor, Log[N(ell,lambda)].
 * N(ell,lambda) = Product[ lambda^2 + n^2, {n,0,ell} ]
 *               = |Gamma(ell + 1 + I lambda)|^2  lambda sinh(Pi lambda) / Pi
 * Assumes ell >= 0.
 */
static
int
legendre_H3d_lnnorm(const int ell, const double lambda, double * result)
{
  double abs_lam = fabs(lambda);

  if(abs_lam == 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(lambda > (ell + 1.0)/GSL_ROOT3_DBL_EPSILON) {
    /* There is a cancellation between the sinh(Pi lambda)
     * term and the log(gamma(ell + 1 + i lambda) in the
     * result below, so we show some care and save some digits.
     * Note that the above guarantees that lambda is large,
     * since ell >= 0. We use Stirling and a simple expansion
     * of sinh.
     */
    double rat = (ell+1.0)/lambda;
    double ln_lam2ell2  = 2.0*log(lambda) + log(1.0 + rat*rat);
    double lg_corrected = -2.0*(ell+1.0) + M_LNPI + (ell+0.5)*ln_lam2ell2 + 1.0/(288.0*lambda*lambda);
    double angle_terms  = lambda * 2.0 * rat * (1.0 - rat*rat/3.0);
    *result = log(abs_lam) + lg_corrected + angle_terms - M_LNPI;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result lg_r;
    gsl_sf_result lg_theta;
    gsl_sf_result ln_sinh;
    gsl_sf_lngamma_complex_impl(ell+1.0, lambda, &lg_r, &lg_theta);
    gsl_sf_lnsinh_impl(M_PI * abs_lam, &ln_sinh);
    *result = log(abs_lam) + ln_sinh.val + 2.0*lg_r.val - M_LNPI;
    return GSL_SUCCESS;
  }
}


/* Calculate series for small eta*lambda.
 * Assumes eta > 0, lambda != 0.
 *
 * This is just the defining hypergeometric for the Legendre function.
 *
 * P^{mu}_{-1/2 + I lam}(z) = 1/Gamma(l+3/2) ((z+1)/(z-1)^(mu/2)
 *                            2F1(1/2 - I lam, 1/2 + I lam; l+3/2; (1-z)/2)
 * We use
 *       z = cosh(eta)
 * (z-1)/2 = sinh^2(eta/2)
 *
 * And recall
 * H3d = sqrt(Pi Norm /(2 lam^2 sinh(eta))) P^{-l-1/2}_{-1/2 + I lam}(cosh(eta))
 */
static
int
legendre_H3d_series(const int ell, const double lambda, const double eta,
                    gsl_sf_result * result)
{
  const int nmax = 5000;
  const double shheta = sinh(0.5*eta);
  const double ln_zp1 = M_LN2 + log(1.0 + shheta*shheta);
  const double ln_zm1 = M_LN2 + 2.0*log(shheta);
  const double zeta = -shheta*shheta;
  gsl_sf_result lg_lp32;
  double term = 1.0;
  double sum  = 1.0;
  gsl_sf_result lnsheta;
  double lnN;
  double lnpre_val, lnpre_err, lnprepow;
  int stat_e;
  int n;

  gsl_sf_lngamma_impl(ell + 3.0/2.0, &lg_lp32);
  gsl_sf_lnsinh_impl(eta, &lnsheta);
  legendre_H3d_lnnorm(ell, lambda, &lnN);
  lnprepow = 0.5*(ell + 0.5) * (ln_zm1 - ln_zp1);
  lnpre_val = lnprepow + 0.5*(lnN + M_LNPI - M_LN2 - lnsheta.val) - lg_lp32.val - log(fabs(lambda));
  lnpre_err = lnsheta.err + lg_lp32.err + GSL_DBL_EPSILON * fabs(lnpre_val);
  for(n=1; n<nmax; n++) {
    double aR = n - 0.5;
    term *= (aR*aR + lambda*lambda)*zeta/(ell + n + 0.5)/n;
    sum  += term;
    if(fabs(term/sum) < 4.0 * GSL_DBL_EPSILON) break;
  }

  stat_e = gsl_sf_exp_mult_err_impl(lnpre_val, lnpre_err, sum, fabs(term), result);
  return GSL_ERROR_SELECT_2(stat_e, (n==nmax ? GSL_EMAXITER : GSL_SUCCESS));
}


/* Evaluate legendre_H3d(ell+1)/legendre_H3d(ell)
 * by continued fraction.
 */
#if 0
static
int
legendre_H3d_CF1(const int ell, const double lambda, const double coth_eta, double * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int n = 1;
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = sqrt(lambda*lambda + (ell+1.0)*(ell+1.0));
  double b1 = (2.0*ell + 3.0) * coth_eta;
  double An = b1*Anm1 + a1*Anm2;
  double Bn = b1*Bnm1 + a1*Bnm2;
  double an, bn;
  double fn = An/Bn;

  while(n < maxiter) {
    double old_fn;
    double del;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    an = -(lambda*lambda + ((double)ell + n)*((double)ell + n));
    bn = (2.0*ell + 2.0*n + 1.0) * coth_eta;
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
    
    if(fabs(del - 1.0) < 4.0*GSL_DBL_EPSILON) break;
  }

  *result = fn;

  if(n == maxiter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}
#endif /* 0 */


/* Evaluate legendre_H3d(ell+1)/legendre_H3d(ell)
 * by continued fraction. Use the Gautschi (Euler)
 * equivalent series.
 */
static
int
legendre_H3d_CF1_ser(const int ell, const double lambda, const double coth_eta, double * result)
{
  const int maxk = 20000;
  double tk   = 1.0;
  double sum  = 1.0;
  double rhok = 0.0;
  int k;
 
  for(k=1; k<maxk; k++) {
    double tlk = (2.0*ell + 1.0 + 2.0*k);
    double l1k = (ell + 1.0 + k);
    double ak = -(lambda*lambda + l1k*l1k)/(tlk*(tlk+2.0)*coth_eta*coth_eta);
    rhok = -ak*(1.0 + rhok)/(1.0 + ak*(1.0 + rhok));
    tk  *= rhok;
    sum += tk;
    if(fabs(tk/sum) < GSL_DBL_EPSILON) break;
  }

  *result = sqrt(lambda*lambda+(ell+1.0)*(ell+1.0))/((2.0*ell+3)*coth_eta)*sum;

  if(k == maxk)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}



/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_legendre_H3d_0_impl(const double lambda, const double eta, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(eta < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(eta == 0.0 || lambda == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    double lam_eta = lambda * eta;
    gsl_sf_result s;
    gsl_sf_sin_err_impl(lam_eta, GSL_DBL_EPSILON * fabs(lam_eta), &s);
    if(eta > -0.5*GSL_LOG_DBL_EPSILON) {
      double f = 2.0 / lambda * exp(-eta);
      result->val  = f * s.val;
      result->err  = fabs(f) * (fabs(s.val * GSL_DBL_EPSILON) + s.err);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
    else {
      double f = 1.0/(lambda*sinh(eta));
      result->val  = f * s.val;
      result->err  = fabs(f) * (fabs(s.val * GSL_DBL_EPSILON) + s.err);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_H3d_1_impl(const double lambda, const double eta, gsl_sf_result * result)
{
  const double xi    = fabs(eta*lambda);
  const double lsq   = lambda*lambda;
  const double lsqp1 = lsq + 1.0;

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(eta < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(eta == 0.0 || lambda == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(xi < GSL_ROOT5_DBL_EPSILON && eta < GSL_ROOT5_DBL_EPSILON) {
    double etasq = eta*eta;
    double xisq  = xi*xi;
    double term1 = (etasq + xisq)/3.0;
    double term2 = -(2.0*etasq*etasq + 5.0*etasq*xisq + 3.0*xisq*xisq)/90.0;
    double sinh_term = 1.0 - eta*eta/6.0 * (1.0 - 7.0/60.0*eta*eta);
    double pre = sinh_term/sqrt(lsqp1) / eta;
    result->val  = pre * (term1 + term2);
    result->err  = pre * GSL_DBL_EPSILON * (fabs(term1) + fabs(term2));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double sin_term;     /*  sin(xi)/xi     */
    double cos_term;     /*  cos(xi)        */
    double coth_term;    /*  eta/tanh(eta)  */
    double sinh_term;    /*  eta/sinh(eta)  */
    double sin_term_err;
    double cos_term_err;
    double pre = sinh_term/sqrt(lsqp1) / eta;
    double term1;
    double term2;
    if(xi < GSL_ROOT5_DBL_EPSILON) {
      sin_term = 1.0 - xi*xi/6.0 * (1.0 - xi*xi/20.0);
      cos_term = 1.0 - 0.5*xi*xi * (1.0 - xi*xi/12.0);
      sin_term_err = GSL_DBL_EPSILON;
      cos_term_err = GSL_DBL_EPSILON;
    }
    else {
      sin_term = sin(xi)/xi;
      cos_term = cos(xi);
      sin_term_err = 2.0 * fabs(cos_term * GSL_DBL_EPSILON);
      cos_term_err = 2.0 * fabs((sin_term * xi) * (xi * GSL_DBL_EPSILON));
    }
    if(eta < GSL_ROOT5_DBL_EPSILON) {
      coth_term = 1.0 + eta*eta/3.0 * (1.0 - eta*eta/15.0);
      sinh_term = 1.0 - eta*eta/6.0 * (1.0 - 7.0/60.0*eta*eta);
    }
    else {
      coth_term = eta/tanh(eta);
      sinh_term = eta/sinh(eta);
    }
    term1 = sin_term*coth_term;
    term1 = cos_term;
    result->val  = pre * (term1 - term2);
    result->err  = pre * (sin_term_err * coth_term + cos_term_err);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_H3d_impl(const int ell, const double lambda, const double eta,
                         gsl_sf_result * result)
{
  const double abs_lam = fabs(lambda);
  const double lsq     = abs_lam*abs_lam;
  const double xi      = abs_lam * eta;
  const double cosh_eta = cosh(eta);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(eta < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(eta > GSL_LOG_DBL_MAX) {
    /* cosh(eta) is too big. */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EOVRFLW;
  }
  else if(ell == 0) {
    return gsl_sf_legendre_H3d_0_impl(lambda, eta, result);
  }
  else if(ell == 1) {
    return gsl_sf_legendre_H3d_1_impl(lambda, eta, result);
  }
  else if(eta == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(xi < 1.0) {
    return legendre_H3d_series(ell, lambda, eta, result);
  }
  else if((ell*ell+lsq)/sqrt(1.0+lsq)/(cosh_eta*cosh_eta) < 5.0*GSL_ROOT3_DBL_EPSILON) {
    /* Large argument.
     */
    gsl_sf_result P;
    double lm;
    int stat_P = gsl_sf_conicalP_large_x_impl(-ell-0.5, lambda, cosh_eta, &P, &lm);
    if(P.val == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return stat_P;
    }
    else {
      double lnN;
      gsl_sf_result lnsh;
      double lnpre_val, lnpre_err;
      int stat_e;
      gsl_sf_lnsinh_impl(eta, &lnsh);
      legendre_H3d_lnnorm(ell, lambda, &lnN);
      lnpre_val  = 0.5*(M_LNPI + lnN - M_LN2 - lnsh.val) - log(abs_lam);
      lnpre_err  = lnsh.err;
      lnpre_err += 2.0 * GSL_DBL_EPSILON * fabs(lnpre_val);
      stat_e = gsl_sf_exp_mult_err_impl(lnpre_val + lm, lnpre_err, P.val, 0.0, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_P);
    }
  }
  else if(abs_lam > 1000.0*ell*ell) {
    /* Large degree.
     */
    gsl_sf_result P;
    double lm;
    int stat_P = gsl_sf_conicalP_xgt1_neg_mu_largetau_impl(ell+0.5,
                                                           lambda,
                                                           cosh_eta, eta,
                                                           &P, &lm);
    if(P.val == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return stat_P;
    }
    else {
      double lnN;
      gsl_sf_result lnsh;
      double lnpre_val, lnpre_err;
      int stat_e;
      gsl_sf_lnsinh_impl(eta, &lnsh);
      legendre_H3d_lnnorm(ell, lambda, &lnN);
      lnpre_val  = 0.5*(M_LNPI + lnN - M_LN2 - lnsh.val) - log(abs_lam);
      lnpre_err  = lnsh.err;
      lnpre_err += 2.0 * GSL_DBL_EPSILON * fabs(lnpre_val);
      stat_e = gsl_sf_exp_mult_err_impl(lnpre_val + lm, lnpre_err, P.val, 0.0, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_P);
    }
  }
  else {
    /* Backward recurrence.
     */
    const double coth_eta = 1.0/tanh(eta);
    double rH;
    int stat_CF1 = legendre_H3d_CF1_ser(ell, lambda, coth_eta, &rH);
    double Hlm1;
    double Hl    = GSL_SQRT_DBL_MIN;
    double Hlp1  = rH * Hl;
    int lp;
    for(lp=ell; lp>0; lp--) {
      double root_term_0 = sqrt(lambda*lambda + (double)lp*lp);
      double root_term_1 = sqrt(lambda*lambda + (lp+1.0)*(lp+1.0));
      Hlm1 = ((2.0*lp + 1.0)*coth_eta*Hl - root_term_1 * Hlp1)/root_term_0;
      Hlp1 = Hl;
      Hl   = Hlm1;
    }

    if(fabs(Hl) > fabs(Hlp1)) {
      gsl_sf_result H0;
      int stat_H0 = gsl_sf_legendre_H3d_0_impl(lambda, eta, &H0);
      result->val  = GSL_SQRT_DBL_MIN/Hl * H0.val;
      result->err  = GSL_SQRT_DBL_MIN/Hl * H0.err;
      result->err += GSL_DBL_EPSILON * ell * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_H0, stat_CF1);
    }
    else {
      gsl_sf_result H1;
      int stat_H1 = gsl_sf_legendre_H3d_1_impl(lambda, eta, &H1);
      result->val = GSL_SQRT_DBL_MIN/Hlp1 * H1.val;
      result->err = GSL_SQRT_DBL_MIN/Hlp1 * H1.err;
      result->err += GSL_DBL_EPSILON * ell * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_H1, stat_CF1);
    }
  }
}


int
gsl_sf_legendre_H3d_array_impl(const int lmax, const double lambda, const double eta, double * result_array)
{
  if(result_array == 0) {
    return GSL_EFAULT;
  }
  else if(eta < 0.0 || lmax < 0) {
    int ell;
    for(ell=0; ell<=lmax; ell++) result_array[ell] = 0.0;
    return GSL_EDOM;
  }
  else if(eta > GSL_LOG_DBL_MAX) {
    /* cosh(eta) is too big. */
    int ell;
    for(ell=0; ell<=lmax; ell++) result_array[ell] = 0.0;
    return GSL_EOVRFLW;
  }
  else if(lmax == 0) {
    gsl_sf_result H0;
    int stat = gsl_sf_legendre_H3d_impl(0, lambda, eta, &H0);
    result_array[0] = H0.val;
    return stat;
  }
  else {
    /* Not the most efficient method. But what the hell... it's simple.
     */
    gsl_sf_result r_Hlp1;
    gsl_sf_result r_Hl;
    int stat_lmax   = gsl_sf_legendre_H3d_impl(lmax,   lambda, eta, &r_Hlp1);
    int stat_lmaxm1 = gsl_sf_legendre_H3d_impl(lmax-1, lambda, eta, &r_Hl);
    int stat_max = GSL_ERROR_SELECT_2(stat_lmax, stat_lmaxm1);

    const double coth_eta = 1.0/tanh(eta);
    int stat_recursion = GSL_SUCCESS;
    double Hlp1 = r_Hlp1.val;
    double Hl   = r_Hl.val;
    double Hlm1;
    int ell;

    result_array[lmax]   = Hlp1;
    result_array[lmax-1] = Hl;

    for(ell=lmax-1; ell>0; ell--) {
      double root_term_0 = sqrt(lambda*lambda + (double)ell*ell);
      double root_term_1 = sqrt(lambda*lambda + (ell+1.0)*(ell+1.0));
      Hlm1 = ((2.0*ell + 1.0)*coth_eta*Hl - root_term_1 * Hlp1)/root_term_0;
      result_array[ell-1] = Hlm1;
      if(!(Hlm1 < DBL_MAX)) stat_recursion = GSL_EOVRFLW;
      Hlp1 = Hl;
      Hl   = Hlm1;
    }

    return GSL_ERROR_SELECT_2(stat_recursion, stat_max);
  }
}
  

/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_legendre_H3d_0_e(const double lambda, const double eta, gsl_sf_result * result)
{
  int status = gsl_sf_legendre_H3d_0_impl(lambda, eta, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_H3d_0_e", status);
  }
  return status;
}

int
gsl_sf_legendre_H3d_1_e(const double lambda, const double eta, gsl_sf_result * result)
{
  int status = gsl_sf_legendre_H3d_1_impl(lambda, eta, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_H3d_1_e", status);
  }
  return status;
}

int
gsl_sf_legendre_H3d_e(const int l, const double lambda, const double eta, gsl_sf_result * result)
{
  int status = gsl_sf_legendre_H3d_impl(l, lambda, eta, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_H3d_e", status);
  }
  return status;
}

int
gsl_sf_legendre_H3d_array_e(const int lmax, const double lambda, const double eta, double * result_array)
{
  int status = gsl_sf_legendre_H3d_array_impl(lmax, lambda, eta, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_H3d_array_e", status);
  }
  return status;
}
