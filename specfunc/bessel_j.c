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

int gsl_sf_bessel_j0_impl(const double x, gsl_sf_result * result)
{
  double ax = fabs(x);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(ax < 0.5) {
    const double y = x*x;
    const double c1 = -1.0/6.0;
    const double c2 =  1.0/120.0;
    const double c3 = -1.0/5040.0;
    const double c4 =  1.0/362880.0;
    const double c5 = -1.0/39916800.0;
    const double c6 =  1.0/6227020800.0;
    result->val = 1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*c6)))));
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    result->val = sin(x)/x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_j1_impl(const double x, gsl_sf_result * result)
{
  double ax = fabs(x);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(ax < 3.1 * DBL_MIN) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(ax < 0.25) {
    const double y = x*x;
    const double c1 = -1.0/10.0;
    const double c2 =  1.0/280.0;
    const double c3 = -1.0/15120.0;
    const double c4 =  1.0/1330560.0;
    const double c5 = -1.0/172972800.0;
    const double sum = 1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*c5))));
    result->val = x/3.0 * sum;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double cos_x = cos(x);
    double sin_x = sin(x);
    result->val = (sin_x/x - cos_x)/x;
    result->err = GSL_DBL_EPSILON*(fabs(sin_x/(x*x)) + fabs(cos_x/x));
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_j2_impl(const double x, gsl_sf_result * result)
{
  double ax = fabs(x);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(ax < 4.0*GSL_SQRT_DBL_MIN) {
    result->val = 0.0;
    result->err = 0.0;
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
    const double sum = 1.0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*(c8+y*c9))))))));
    result->val = y/15.0 * sum;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double cos_x = cos(x);
    double sin_x = sin(x);
    result->val =  ((3.0/(x*x) - 1.0) * sin_x - 3.0*cos_x/x)/x;
    result->err = GSL_DBL_EPSILON * (fabs(sin_x/x) + 3.0*fabs(cos_x/(x*x)));
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_jl_impl(const int l, const double x, gsl_sf_result * result)
{
  if(l < 0 || x < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
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
  else if(x*x < 10.0*(l+0.5)/M_E) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_Jnu_taylor_impl(l+0.5, x, -1, 50, GSL_DBL_EPSILON, &b);
    result->val = sqrt((0.5*M_PI)/x) * b;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return status;
  }
  else if(GSL_ROOT3_DBL_EPSILON * x > (l*l + l + 1.0)) {
    gsl_sf_result b;
    int status = gsl_sf_bessel_Jnu_asympx_impl(l + 0.5, x, &b);
    double pre = sqrt((0.5*M_PI)/x);
    result->val = pre * b.val;
    result->err = GSL_DBL_EPSILON * fabs(result->val) + pre * b.err;
    return status;
  }
  else if(l > 1.0/GSL_ROOT6_DBL_EPSILON) {
    gsl_sf_result b;
    int status = gsl_sf_bessel_Jnu_asymp_Olver_impl(l + 0.5, x, &b);
    double pre = sqrt((0.5*M_PI)/x);
    result->val = pre * b.val;
    result->err = GSL_DBL_EPSILON * fabs(result->val) + pre * b.err;
    return status;
  }
  else {
    /* recurse down from safe values */
    double rt_term = sqrt((0.5*M_PI)/x);
    const int LMAX = 2 + (int)(1.0/GSL_ROOT6_DBL_EPSILON);

    gsl_sf_result r_jellp1;
    gsl_sf_result r_jell;
    int stat_0 = gsl_sf_bessel_Jnu_asymp_Olver_impl(LMAX + 1 + 0.5, x, &r_jellp1);
    int stat_1 = gsl_sf_bessel_Jnu_asymp_Olver_impl(LMAX     + 0.5, x, &r_jell);
    double jellp1 = r_jellp1.val;
    double jell   = r_jell.val;
    double jellm1;
    int ell;

    jellp1 *= rt_term;
    jell   *= rt_term;
    for(ell = LMAX; ell >= l+1; ell--) {
      jellm1 = -jellp1 + (2*ell + 1)/x * jell;
      jellp1 = jell;
      jell   = jellm1;
    }

    result->val = jellm1;
    result->err = fabs(result->val)*(fabs(r_jellp1.err/r_jellp1.val) + fabs(r_jell.err/r_jell.val));

    return GSL_ERROR_SELECT_2(stat_0, stat_1);
  }
}


int gsl_sf_bessel_jl_array_impl(const int lmax, const double x, double * result_array)
{
  if(result_array == 0) {
    return GSL_EFAULT;
  }
  else if(lmax < 0 || x < 0.0) {
    int j;
    for(j=0; j<=lmax; j++) result_array[j] = 0.0;
    return GSL_EDOM;
  }
  else {
    gsl_sf_result r_jellp1;
    gsl_sf_result r_jell;
    int stat_0 = gsl_sf_bessel_jl_impl(lmax+1, x, &r_jellp1);
    int stat_1 = gsl_sf_bessel_jl_impl(lmax,   x, &r_jell);
    double jellp1 = r_jellp1.val;
    double jell   = r_jell.val;
    double jellm1;
    int ell;

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
  if(jl_x == 0) {
    return GSL_EFAULT;
  }
  else if(lmax < 0 || x < 0.0) {
    int j;
    for(j=0; j<=lmax; j++) jl_x[j] = 0.0;
    return GSL_EDOM;
  }
  else if(x < 2.0*GSL_ROOT4_DBL_EPSILON) {
    /* first two terms of Taylor series */
    double inv_fact = 1.0;  /* 1/(1 3 5 ... (2l+1)) */
    double x_l      = 1.0;  /* x^l */
    int l;
    for(l=0; l<=lmax; l++) {
      jl_x[l]  = x_l * inv_fact;
      jl_x[l] *= 1.0 - 0.5*x*x/(2.0*l+3.0);
      inv_fact /= 2.0*l+3.0;
      x_l      *= x;
    }
    return GSL_SUCCESS;
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
    while(fabs(del) >= fabs(FP) * GSL_DBL_EPSILON);
    
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

    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_j0_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_j0_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_j0_e", status);
  }
  return status;
}

int gsl_sf_bessel_j1_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_j1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_j1_e", status);
  }
  return status;
}

int gsl_sf_bessel_j2_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_j2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_j2_e", status);
  }
  return status;
}

int gsl_sf_bessel_jl_e(const int l, const double x, gsl_sf_result * result)
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
