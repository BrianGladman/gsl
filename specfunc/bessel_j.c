/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_olver.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_trig.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_j0_impl(const double x, double * result)
{
  double ax = fabs(x);

  if(ax < 0.1) {
    double x2 = x*x;
    *result = 1.0 - x2/6.0 * (1.0 - x2/20.0 * (1.0 - x2/42.0 * (1.0 - x2/72.0)));
    return GSL_SUCCESS;
  }
  else {
    double arg = x;
    int stat = gsl_sf_angle_restrict_pos_impl(&arg);
    *result = sin(arg)/x;
    return stat;
  }
  return GSL_SUCCESS;
}

int gsl_sf_bessel_j1_impl(const double x, double * result)
{
  double ax = fabs(x);

  if(ax < 3.1 * DBL_MIN && x != 0.0) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(ax < 0.05) {
    double x2 = x*x;
    *result = x/3.0 * (1.0 - x2/10.0 * (1.0 - x2/28.0 * (1.0 - x2/54.0)));
    return GSL_SUCCESS;
  }
  else {
    double arg = x;
    int stat = gsl_sf_angle_restrict_pos_impl(&arg);
    double cos_x = cos(arg);
    double sin_x = sin(arg);
    *result = (sin_x/x - cos_x)/x;
    return stat;
  }
}

int gsl_sf_bessel_j2_impl(const double x, double * result)
{
  double ax = fabs(x);

  if(ax < 4.0*GSL_SQRT_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(ax < 1.3) {
    const double y  = x*x;
    const double c1 = -1.0/14.0;
    const double c2 =  1.0/504.0;
    const double c3 = -1.0/33264.0;
    const double c4 =  1.0/3459456.0;
    const double c5 = -1.0/518918400;
    const double c6 =  1.0/105859353600.0;
    const double c7 = -1.0/28158588057600.0;
    const double c8 =  1.0/9461285587353600.0;
    const double c9 = -1.0/3916972233164390400.0;
    *result = y/15.0*(1.0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*(c8+y*c9)))))))));
    return GSL_SUCCESS;
  }
  else {
    double arg = x;
    int stat = gsl_sf_angle_restrict_pos_impl(&arg);
    double cos_x = cos(arg);
    double sin_x = sin(arg);
    *result =  ((3.0/(x*x) - 1.0) * sin_x - 3.0*cos_x/x)/x;
    return stat;
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
  else if(x*x < 10.0*(l+1.5)*GSL_ROOT5_MACH_EPS) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_Jnu_taylor_impl(l+0.5, x, -1, 50, GSL_MACH_EPS, &b);
    *result = sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(x*x < 10.0*(l+0.5)/M_E) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_Jnu_taylor_impl(l+0.5, x, -1, 50, GSL_MACH_EPS, &b);
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
    gsl_sf_bessel_Jnu_asymp_Olver_impl(LMAX + 1 + 0.5, x, &jellp1);
    gsl_sf_bessel_Jnu_asymp_Olver_impl(LMAX     + 0.5, x, &jell);
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
  if(lmax < 0 || x < 0.0) {
    int j;
    for(j=0; j<=lmax; j++) result_array[j] = 0.0;
    return GSL_EDOM;
  }
  else {
    int ell;
    double jellp1, jell, jellm1;
    int stat_0 = gsl_sf_bessel_jl_impl(lmax+1, x, &jellp1);
    int stat_1 = gsl_sf_bessel_jl_impl(lmax,   x, &jell);
    result_array[lmax] = jell;
    for(ell = lmax; ell >= 1; ell--) {
      jellm1 = -jellp1 + (2*ell + 1)/x * jell;
      jellp1 = jell;
      jell   = jellm1;
      result_array[ell-1] = jellm1;
    }
    return GSL_ERROR_SELECT_2(stat_0, stat_1);
  }
}

int gsl_sf_bessel_jl_steed_array_impl(const int lmax, const double x, double * jl_x)
{
  if(lmax < 0 || x < 0.0) {
    int j;
    for(j=0; j<=lmax; j++) jl_x[j] = 0.0;
    return GSL_EDOM;
  }

  if(x < 2.0*GSL_ROOT4_MACH_EPS) {
    /* first two terms of Taylor series */
    double inv_fact = 1.0;  /* 1/(1 3 5 ... (2l+1)) */
    double x_l      = 1.0;  /* x^l */
    int l;
    for(l=0; l<=lmax; l++) {
      jl_x[l]  = x_l * inv_fact;
      jl_x[l] *= 1.0 - 0.5*x*x/(2.0*l+3.0);
      inv_fact /= 2.0*l+3.;
      x_l      *= x;
    }
  }
  else {
    /* Steed/Barnett algorithm [Comp. Phys. Comm. 21, 297 (1981)] */
    double x_inv = 1.0/x;
    double W = 2.0*x_inv;
    double F = 1.0;
    double FP = (lmax+1.0) * x_inv;
    double B = 2.0*FP + x_inv;
    double end = B + 20000.0*W;
    double D = 1.0/B;
    double del = -D;
    
    FP += del;
    
    /* continued fraction */
    do {
      B += W;
      D = 1.0/(B-D);
      del *= (B*D - 1.);
      FP += del;
      if(D < 0.0) F = -F;
      if(B > end) {
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
