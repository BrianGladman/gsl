/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"


#define Root_2OverPi_  0.797884560802865355879892


/* checked OK [GJ] */
int gsl_sf_legendre_sph_con_reg_m1_impl(double lambda,
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


/* checked OK [GJ] */
int gsl_sf_legendre_sph_con_reg_0_impl(double lambda,
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

int gsl_sf_legendre_sph_con_reg_array_impl(int lmax, double lambda, double x, double * result, double * harvest)
{
  double nu_nup1 = -0.25 - lambda*lambda;
  double X = x/sqrt(x*x-1.);
  double y2, y1, y0;
  int ell;
  double ell_fact = 1.;

  gsl_sf_legendre_sph_con_reg_m1_impl(lambda, 1-x, 1+x, &y2);
  gsl_sf_legendre_sph_con_reg_0_impl(lambda, 1-x, 1+x, &y1);
  harvest[0] = y1;

  for(ell=1; ell<=lmax; ell++) {
    y0 = (-y2 - 2*(0.5-ell)*X*y1) / (ell*ell - 0.25 - nu_nup1);
printf("%22.17g   %22.17g    %22.17g\n", y2, -2.*(0.5-ell)*X*y1, -y2 - 2*(0.5-ell)*X*y1);
    y2 = y1;
    y1 = y0;
    harvest[ell] = y0;
  }
  
  *result = y0;
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
