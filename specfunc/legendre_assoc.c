/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"

#include "recurse.h"

#define Root_2OverPi_  0.797884560802865355879892


/* P_{-1/2 + I lambda}^{1/2}
 *
 * checked OK [GJ]
 */
int gsl_sf_conical_sph_irr_1_impl(double lambda,
                                  double one_minus_x,
                                  double one_plus_x,
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


/* P_{-1/2 + I lambda}{-1/2} (x)
 *
 * checked OK [GJ] 
 */
int gsl_sf_conical_sph_reg_0_impl(double lambda,
                                  double one_minus_x,
                                  double one_plus_x,
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
  else if(x == 1.) {
    *result = 0.;
    return GSL_SUCCESS;
  }  
  else {
    /* x <= -1 */
    return GSL_EDOM;
  }
}


/* P_{-1/2 + I lambda}^{-1/2-n} (x)
 *
 * p[0] = lambda^2
 * p[1] = x/Sqrt(x^2 - 1)
 */
#define REC_COEFF_A(n,p) (-(2.*n+1.)*p[1]/((n+1)*(n+1)+p[0]))
#define REC_COEFF_B(n,p) (1./((n+1)*(n+1)+p[0]))
GEN_RECURSE_BACKWARD_MINIMAL_SIMPLE(conical_sph_reg)
#undef REC_COEFF_A
#undef REC_COEFF_B


int gsl_sf_conical_sph_reg_array_impl(int lmax, double lambda, double x, double * result, double * harvest)
{
  double f0;
  double p[2];
  p[0] = lambda*lambda;
  p[1] = x/sqrt(x*x-1.);
  gsl_sf_conical_sph_reg_0_impl(lambda, 1-x, 1+x, &f0);
  recurse_backward_minimal_simple_conical_sph_reg(lmax+30, lmax, 0, p, f0, harvest, result);
}




int gsl_sf_hyper_0_impl(double lambda, double x, double * result)
{
  *result = sin(lambda*x)/(lambda*sinh(x));
}

int gsl_sf_hyper_1_impl(double lambda, double x, double * result)
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
