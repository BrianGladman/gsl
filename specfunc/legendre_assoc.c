/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"


#define Root_2OverPi_  0.797884560803


int gsl_sf_legendre_sph_con_reg_m1_impl(double lambda, double x, double * result)
{
  if(fabs(x) < 1.) {
    double ac  = acos(x);
    double den = sqrt(sqrt(1.-x*x));
    *result = Root_2OverPi_ / den * cosh(ac * lambda);
    return GSL_SUCCESS;
  }
  else if(fabs(x) > 1.) {
    double ln_term = log(x + sqrt(x*x-1.));
    double den = sqrt(sqrt(x*x-1.));
    *result = Root_2OverPi_ / den * cos(lambda * ln_term);
    return GSL_SUCCESS;
  }
  else {
    /* |x| == 1 */
    return GSL_EDOM;
  }
}

int gsl_sf_legendre_sph_con_reg_0_impl(double lambda, double x, double * result)
{
  if(fabs(x) < 1.) {
    double ac  = acos(x);
    double den = sqrt(sqrt(1.-x*x));
    double arg = ac * lambda;
    if(fabs(arg) < GSL_SQRT_MACH_EPS) {
      *result = Root_2OverPi_ / den * ac;
    }
    else {
      *result = Root_2OverPi_ / (den*lambda) * sinh(arg);
    }
    return GSL_SUCCESS;
  }
  else if(x > 1.) {
    doubel sq_term = sqrt(x*x-1.);
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
  else {
    /* x == 1  or  x <= -1 */
    return GSL_EDOM;
  }
}
