/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_j0_impl(const double x, double * result)
{
  if(fabs(x) < GSL_ROOT4_MACH_EPS) {
    *result = 1. - x*x/6.;
  }
  else {
    *result = sin(x)/x;
  }
  return GSL_SUCCESS;
}

int gsl_sf_bessel_j1_impl(const double x, double * result)
{
  if(fabs(x) < 3.*DBL_MIN) {
    *result = 0.;
    return GSL_EUNDRFLW;
  }
  else if(fabs(x) < 2.*GSL_ROOT4_MACH_EPS) {
    *result = x/3. * (1. - x*x/10.);
    return GSL_SUCCESS;
  }
  else {
    double cos_x = cos(x);
    double sin_x = sin(x);
    *result = sin_x/(x*x) - cos_x/x;
    return GSL_SUCCESS;
  }
}

int gsl_sf_bessel_j2_impl(const double x, double * result)
{
  if(fabs(x) < GSL_SQRT_DBL_MIN) {
    *result = 0.;
    return GSL_EUNDRFLW;
  }
  else if(fabs(x) < 2.*GSL_ROOT4_MACH_EPS) {
    *result = x*x/15. * (1. - x*x/14.);
    return GSL_SUCCESS;
  }
  else {
    double x2 = x*x;
    double cos_x = cos(x);
    double sin_x = sin(x);
    *result =  (3./x2 - 1.) * sin_x/x - 3.*cos_x/x2;
    return GSL_SUCCESS;
  }
}  

int gsl_sf_bessel_jl_impl(const int l, const double x, double * result)
{
  if(l < 0 || x < 0.0) {
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(x*x < 10.*(l+1.5)*GSL_ROOT5_MACH_EPS) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_Jnu_taylor_impl(l+0.5, x, -1, 4, &b);
    *result = sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(GSL_ROOT3_MACH_EPS * x > (l*l + l + 1)) {
    double b = 0.0;
    int status = gsl_sf_bessel_Jnu_asympx_impl(l + 0.5, x, &b);
    *result = sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(l > 30) {
    double b = 0.0;
    int status = gsl_sf_bessel_Jnu_asymp_Olver_impl(l + 0.5, x, &b);
    *result = sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(l == 0) {
    return gsl_sf_bessel_j0_impl(x, result);
  }
  else if(l == 1) {
    return gsl_sf_bessel_j1_impl(x, result);
  }
  else if(l == 2) {
    return gsl_sf_bessel_j2_impl(x, result);
  }
  else {
    /* recurse down from safe values */
    double rt_term = sqrt(M_PI/(2.0*x));
    double jellp1, jell, jellm1;
    const int LMAX = 31;
    int ell;
    gsl_sf_bessel_asymp_Jnu_Olver_impl(LMAX + 1 + 0.5, x, &jellp1);
    gsl_sf_bessel_asymp_Jnu_Olver_impl(LMAX     + 0.5, x, &jell);
    jellp1 *= rt_term;
    jell   *= rt_term;
    for(ell = LMAX; ell >= l+1; ell--) {
      jellm1 = -jellp1 + (2*ell + 1)/x * jell;
      jellp1 = jell;
      jell   = jellm1;
    }
    *result = jellm1;
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_jl_array_impl(const int lmax, const double x, double * result_array)
{
  int ell;
  double jellp1, jell, jellm1;
  gsl_sf_bessel_jl_impl(lmax+1, x, &jellp1);
  gsl_sf_bessel_jl_impl(lmax,   x, &jell);
  result_array[lmax] = jell;
  for(ell = lmax; ell >= 1; ell--) {
    jellm1 = -jellp1 + (2*ell + 1)/x * jell;
    jellp1 = jell;
    jell   = jellm1;
    result_array[ell-1] = jellm1;
  }
  return GSL_SUCCESS;
}

int gsl_sf_bessel_j_steed_array_impl(const int lmax, const double x, double * jl_x)
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
    /* Steed/Barnett algorithm [Comp. Phys. Comm. 21, 297 (1981)] */
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
	return GSL_EMAXITER;
      }
    }
    while(fabs(del) >= fabs(FP) * GSL_MACH_EPS);
    
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

int gsl_sf_bessel_j0_e(const double x, double * result)
{
  int status = gsl_sf_bessel_j0_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_j0_e", status);
  }
  return status;
}

int gsl_sf_bessel_j1_e(const double x, double * result)
{
  int status = gsl_sf_bessel_j1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_j1_e", status);
  }
  return status;
}

int gsl_sf_bessel_j2_e(const double x, double * result)
{
  int status = gsl_sf_bessel_j2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_j2_e", status);
  }
  return status;
}

int gsl_sf_bessel_jl_e(const int l, const double x, double * result)
{
  int status = gsl_sf_bessel_jl_impl(l, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_jl_e", status);
  }
  return status;
}

int gsl_sf_bessel_jl_array_e(const int lmax, const double x, double * jl_array)
{
  int status = gsl_sf_bessel_jl_array_impl(lmax, x, jl_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_jl_array_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_j0(const double x)
{
  double y;
  int status = gsl_sf_bessel_j0_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_j0", status);
  }
  return y;
}

double gsl_sf_bessel_j1(const double x)
{
  double y;
  int status = gsl_sf_bessel_j1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_j1", status);
  }
  return y;
}

double gsl_sf_bessel_j2(const double x)
{
  double y;
  int status = gsl_sf_bessel_j2_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_j2", status);
  }
  return y;
}

double gsl_sf_bessel_jl(const int l, const double x)
{
  double y;
  int status = gsl_sf_bessel_jl_impl(l, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_jl", status);
  }
  return y;
}
