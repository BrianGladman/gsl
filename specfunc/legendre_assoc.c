/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"

#include "recurse.h"

#define Root_2OverPi_  0.797884560802865355879892


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* P_{-1/2 + I lambda}^{1/2} (x)
 * [Abramowitz+Stegun 8.6.8, 8.6.12]
 * checked OK [GJ]
 */
int gsl_sf_conical_sph_irr_1_impl(const double lambda,
                                  const double one_minus_x,
                                  const double one_plus_x,
                                  double * result
                                  )
{
  double x = 0.5 * (one_plus_x - one_minus_x);
  if(fabs(x) < 1.) {
    double ac  = acos(x);
    double den = sqrt(sqrt(one_minus_x*one_plus_x));
    *result = Root_2OverPi_ / den * cosh(ac * lambda);
    return GSL_SUCCESS;
  }
  else if(one_minus_x < 0.) { /* x > 1. */
    double sq_term = sqrt(-one_minus_x*one_plus_x);
    double ln_term = log(x + sq_term);
    double den = sqrt(sq_term);
    *result = Root_2OverPi_ / den * cos(lambda * ln_term);
    return GSL_SUCCESS;
  }
  else {
    /* |x| == 1 */
    return GSL_EDOM;
  }
}


/* P_{-1/2 + I lambda}^{-1/2} (x)
 * [Abramowitz+Stegun 8.6.9, 8.6.14]
 * checked OK [GJ] 
 */
int gsl_sf_conical_sph_reg_0_impl(const double lambda,
                                  const double one_minus_x,
                                  const double one_plus_x,
                                  double * result
                                  )
{
  double x = 0.5 * (one_plus_x - one_minus_x);
  if(fabs(x) < 1.) {
    double ac  = acos(x);
    double den = sqrt(sqrt(one_minus_x*one_plus_x));
    double arg = ac * lambda;
    if(fabs(arg) < GSL_SQRT_MACH_EPS) {
      *result = Root_2OverPi_ / den * ac;
    }
    else {
      *result = Root_2OverPi_ / (den*lambda) * sinh(arg);
    }
    return GSL_SUCCESS;
  }
  else if(one_minus_x < 0.) { /* x > 1. */
    double sq_term = sqrt(-one_minus_x*one_plus_x);
    double ln_term = log(x + sq_term);
    double den = sqrt(sq_term);
    double arg = lambda * ln_term;
    if(arg < GSL_SQRT_MACH_EPS) {
      *result = Root_2OverPi_ / den * ln_term;
    }
    else {
      *result = Root_2OverPi_ / (den*lambda) * sin(arg);
    }
    return GSL_SUCCESS;
  }
  else if(one_minus_x == 0.) { /* x == 1. */
    *result = 0.;
    return GSL_SUCCESS;
  }  
  else {
    /* x <= -1 */
    return GSL_EDOM;
  }
}

/* P_{-1/2 + I lambda}^mu (x)     mu <= 1/2, lambda >= 0, |1-x| < 2, x > 1
 * [Abramowitz+Stegun, 8.1.2]
 */
static double conical_reg_series(const int N, const double mu, const double lambda, const double x)
{
  int i;
  const double mu_term = 1. - mu;
  const double lam2    = lambda*lambda;
  const double arg     = 0.5 * (1-x);
  double term      = (0.25 + lam2)/mu_term * arg;
  double F_series  = 1. + term;
  for(i=2; i<=N; i++) {
    int odd = 2*i+1;
    term     *= (0.25*odd*odd + lam2)/(mu_term + 1.) * arg;
    F_series += term;
  }
  return exp(-gsl_sf_lngamma(mu_term)) * pow((x-1)/(x+1),-0.5*mu) * F_series;
}

/* P_{-1/2 + I lambda}^{-1/2 - n} (x)    x > 1
 *
 * p[0] = lambda^2
 * p[1] = x/Sqrt(x^2 - 1)
 */
#define REC_COEFF_A(n,p) (-(2.*(n)+1.)*p[1]/(((n)+1)*((n)+1)+p[0]))
#define REC_COEFF_B(n,p) (1./(((n)+1)*((n)+1)+p[0]))
GEN_RECURSE_BACKWARD_MINIMAL_SIMPLE(conical_sph_reg_xgt1)
#undef REC_COEFF_A
#undef REC_COEFF_B

/* P_{-1/2 + I lambda}^{-1/2 - n} (x)    -1 < x < 1
 *
 * p[0] = lambda^2
 * p[1] = x/Sqrt(1 - x^2)
 */
#define REC_COEFF_A(n,p) ((2.*(n)+1.)*p[1]/(((n)+1)*((n)+1)+p[0]))
#define REC_COEFF_B(n,p) (-1./(((n)+1)*((n)+1)+p[0]))
GEN_RECURSE_FORWARD_SIMPLE(conical_sph_reg)
#undef REC_COEFF_A
#undef REC_COEFF_B

/* P_{-1/2 + I lambda}^{-1/2 - n} (x)    -1 < x < 1
 *
 * p[0] = lambda^2
 * p[1] = x/Sqrt(1 - x^2)
 */
#define REC_COEFF_A(n,p) ((2.*(n)+1.)*p[1]/(((n)+1)*((n)+1)+p[0]))
#define REC_COEFF_B(n,p) (-1./(((n)+1)*((n)+1)+p[0]))
GEN_RECURSE_BACKWARD_MINIMAL_SIMPLE(conical_sph_reg_xlt1)
#undef REC_COEFF_A
#undef REC_COEFF_B


int gsl_sf_conical_sph_reg_impl(const int lmax, const double lambda,
                                const double one_m_x, const double one_p_x,
                                double * result, double * harvest
				)
{
  double x = 0.5 * (one_p_x - one_m_x);

  if(fabs(x) < 1.) {
    double f0;
    double p[2];
    p[0] = lambda*lambda;
    p[1] = x/sqrt(one_m_x*one_p_x);
    gsl_sf_conical_sph_reg_0_impl(lambda, one_m_x, one_p_x, &f0);  /* l =  0  */ 
    recurse_backward_minimal_simple_conical_sph_reg_xlt1(lmax+30, lmax, 0, p, f0, harvest, result);
    
  }
  else if(x > 1.) {
    double f0;
    double p[2];
    p[0] = lambda*lambda;
    p[1] = x/sqrt(-one_m_x*one_p_x);
    gsl_sf_conical_sph_reg_0_impl(lambda, one_m_x, one_p_x, &f0);
    recurse_backward_minimal_simple_conical_sph_reg_xgt1(lmax+30, lmax, 0, p, f0, harvest, result);
  }
  else {
    return GSL_EDOM;
  }
}





int gsl_sf_hyper_0_impl(const double lambda, const double x, double * result)
{
  *result = sin(lambda*x)/(lambda*sinh(x));
}

int gsl_sf_hyper_1_impl(const double lambda, const double x, double * result)
{
  *result = sin(lambda*x)/(lambda*sinh(x)) 
	  /sqrt(lambda*lambda+1.) * (1./tanh(x) - lambda/tan(lambda*x));
}

int gsl_sf_hyper_array_impl(int lmax, double lambda, double x, double * result, double * harvest)
{
  double X = 1./tanh(x);
  double y2, y1, y0;
  int ell;

  gsl_sf_hyper_0_impl(lambda, x, &y2);
  gsl_sf_hyper_1_impl(lambda, x, &y1);

  harvest[0] = y2;
  harvest[1] = y1;

  for(ell=2; ell<=lmax; ell++) {
    double a = sqrt(lambda*lambda + ell*ell);
    double b = sqrt(lambda*lambda + (ell-1)*(ell-1));
    y0 = ((2*ell-1)*X*y1 - b*y2) / a;
    y2 = y1;
    y1 = y0;
    harvest[ell] = y0;
  }
  
  *result = y0;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_conical_sph_irr_1_e(const double lambda, const double x, double * result)
{
  int status = gsl_sf_conical_sph_irr_1_impl(lambda, 1.-x, 1.+x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_sph_irr_1_e", status);
  }
  return status;
}

int gsl_sf_conical_sph_reg_0_e(const double lambda, const double x, double * result)
{
  int status = gsl_sf_conical_sph_reg_0_impl(lambda, 1.-x, 1.+x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_sph_reg_0_e", status);
  }
  return status;
}

int gsl_sf_conical_sph_reg_e(const int l, const double lambda, const double x, double * result)
{
  int status = gsl_sf_conical_sph_reg_impl(l, lambda, 1.-x, 1.+x, result, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_sph_reg_e", status);
  }
  return status;
}

int gsl_sf_conical_sph_reg_array_e(const int l, const double lambda, const double x, double * result_array)
{
  double y;
  int status = gsl_sf_conical_sph_reg_impl(l, lambda, 1.-x, 1.+x, &y, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_sph_reg_array_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_conical_sph_irr_1(const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conical_sph_irr_1_impl(lambda, 1.-x, 1.+x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conical_sph_irr_1_e");
  }
  return y;
}

double gsl_sf_conical_sph_reg_0(const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conical_sph_reg_0_impl(lambda, 1.-x, 1.+x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conical_sph_reg_0_e");
  }
  return y;
}

double gsl_sf_conical_sph_reg(const int l, const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conical_sph_reg_impl(l, lambda, 1.-x, 1.+x, &y, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conical_sph_reg");
  }
  return y;
}
