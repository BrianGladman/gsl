/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_bessel.h"

#define ACC GSL_SQRT_MACH_EPS
#define RootPiOver2_  0.886226925453
#define Gamma1pt5_    RootPiOver2_


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_j_steed_array_impl(double x, int lmax, double * jl_x)
{
  if(lmax < 0 || x < 0.) {
    return GSL_EDOM;
  }

  if(x < 2.*GSL_ROOT4_MACH_EPS) {
    /* first two terms of Taylor series */
    double inv_fact = 1.;  /* 1/(1 3 5 ... (2l+1)) */
    double x_l      = 1.;  /* x^l */
    int l;
    for(l=0; l<=lmax; l++) {
      jl_x[l]  = x_l * inv_fact;
      jl_x[l] *= 1. - 0.5*x*x/(2.*l+3.);
      inv_fact /= 2.*l+3.;
      x_l      *= x;
    }
  }
  else {
    /* Steed/Barnett algorithm */
    double x_inv = 1./x;
    double W = 2.*x_inv;
    double F = 1.;
    double FP = (lmax+1.) * x_inv;
    double B = 2.*FP + x_inv;
    double end = B + 20000.*W;
    double D = 1./B;
    double del = -D;
    
    FP += del;
    
    /* continued fraction */
    do {
      B += W;
      D = 1./(B-D);
      del *= (B*D - 1.);
      FP += del;
      if(D < 0.) F = -F;
      if(B > end) {
	/*
	GSL_ERROR_RETURN(
		"gsl_sf_bessel_j_steed: continued fraction not converging",
		GSL_EFAILED,
		0.
		);
		*/
	return GSL_ETIMEOUT;
      }
    }
    while(fabs(del) >= fabs(FP) * ACC);
    
    FP *= F;
    
    if(lmax > 0) {
      /* downward recursion */
      double XP2 = FP;
      double PL = lmax * x_inv;
      int L  = lmax;
      int LP;
      jl_x[lmax] = F;
      for(LP = 1; LP<=lmax; LP++) {
	jl_x[L-1] = PL * jl_x[L] + XP2;
	FP = PL*jl_x[L-1] - jl_x[L];
	XP2 = FP;
	PL -= x_inv;
	--L;
      }
      F = jl_x[0];
    }
    
    /* normalization */
    W = x_inv / sqrt(FP*FP + F*F);
    jl_x[0] = W*F;
    if(lmax > 0) {
      int L;
      for(L=1; L<=lmax; L++) {
	jl_x[L] *= W;
      }
    }
  }
  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_j_steed_array_e(double x, int lmax, double * jl_x)
{
  int status = gsl_sf_bessel_j_steed_array_impl(x, lmax, jl_x);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_j_steed_array_e", status);
  }
  return status;
}
