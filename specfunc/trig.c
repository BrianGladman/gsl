/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_log.h"
#include "gsl_sf_trig.h"


/* sinh(x) series
 * double-precision for |x| < 1.0
 */
#ifdef HAVE_INLINE
inline
#endif
static
int
sinh_series(const double x, double * result)
{
  const double y = x*x;
  const double c0 = 1.0/6.0;
  const double c1 = 1.0/120.0;
  const double c2 = 1.0/5040.0;
  const double c3 = 1.0/362880.0;
  const double c4 = 1.0/39916800.0;
  const double c5 = 1.0/6227020800.0;
  const double c6 = 1.0/1307674368000.0;
  const double c7 = 1.0/355687428096000.0;
  *result = x*(1.0 + y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))))));
  return GSL_SUCCESS;
}


/* cosh(x)-1 series
 * double-precision for |x| < 1.0
 */
#ifdef HAVE_INLINE
inline
#endif
static
int
cosh_m1_series(const double x, double * result)
{
  const double y = x*x;
  const double c0 = 0.5;
  const double c1 = 1.0/24.0;
  const double c2 = 1.0/720.0;
  const double c3 = 1.0/40320.0;
  const double c4 = 1.0/3628800.0;
  const double c5 = 1.0/479001600.0;
  const double c6 = 1.0/87178291200.0;
  const double c7 = 1.0/20922789888000.0;
  const double c8 = 1.0/6402373705728000.0;
  *result = y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*c8))))))));
  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_complex_sin_impl(const double zr, const double zi,
                        gsl_sf_result * szr, gsl_sf_result * szi)
{
  if(szr == 0 || szi == 0) {
    return GSL_EFAULT;
  }
  else if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    szr->val = sin(zr)*(ch_m1 + 1.0);
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    szr->val = sin(zr)*ch;
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else {
    szr->val = 0.0;
    szi->val = 0.0;
    szr->err = 0.0;
    szi->err = 0.0;
    return GSL_EOVRFLW;
  }
}


int
gsl_sf_complex_cos_impl(const double zr, const double zi,
                        gsl_sf_result * czr, gsl_sf_result * czi)
{
  if(czr == 0 || czi == 0) {
    return GSL_EFAULT;
  }
  else if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    czr->val =  cos(zr)*(ch_m1 + 1.0);
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    czr->val =  cos(zr)*ch;
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else {
    czr->val = 0.0;
    czi->val = 0.0;
    czr->err = 0.0;
    czi->err = 0.0;
    return GSL_EOVRFLW;
  }
}


int
gsl_sf_complex_logsin_impl(const double zr, const double zi,
                           gsl_sf_result * lszr, gsl_sf_result * lszi)
{
  if(lszr == 0 || lszi == 0) {
    return GSL_EFAULT;
  }
  else if(zi > 60.0) {
    lszr->val = -M_LN2 + zi;
    lszi->val =  0.5*M_PI - zr;
    lszr->err = GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else if(zi < -60.0) {
    lszr->val = -M_LN2 - zi;
    lszi->val = -0.5*M_PI + zr;
    lszr->err = GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else {
    gsl_sf_result sin_r, sin_i;
    int status;
    gsl_sf_complex_sin_impl(zr, zi, &sin_r, &sin_i); /* ok by construction */
    status = gsl_sf_complex_log_impl(sin_r.val, sin_i.val, lszr, lszi);
    if(status == GSL_EDOM) {
      lszr->val = 0.0;
      lszi->val = 0.0;
      lszr->err = 0.0;
      lszi->err = 0.0;
      return GSL_EDOM;
    }
  }
  return gsl_sf_angle_restrict_symm_impl(&(lszi->val));
}


int
gsl_sf_lnsinh_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(fabs(x) < 1.0) {
    double eps;
    sinh_series(x, &eps);
    result->val = log(eps);
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = x + log(0.5*(1.0 - exp(-2.0*x)));
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_lncosh_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(fabs(x) < 1.0) {
    double eps;
    cosh_m1_series(x, &eps);
    return gsl_sf_log_1plusx_impl(eps, result);
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = x + log(0.5*(1.0 + exp(-2.0*x)));
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/*
inline int gsl_sf_sincos_impl(const double theta, double * s, double * c)
{
  double tan_half = tan(0.5 * theta);
  double den = 1. + tan_half*tan_half;
  double cos_theta = (1.0 - tan_half*tan_half) / den;
  double sin_theta = 2.0 * tan_half / den;
}
*/

int
gsl_sf_polar_to_rect_impl(const double r, const double theta,
                          gsl_sf_result * x, gsl_sf_result * y)
{
  double t   = theta;
  int status = gsl_sf_angle_restrict_symm_impl(&t);
  x->val = r * cos(t);
  y->val = r * sin(t);
  x->err = GSL_DBL_EPSILON * fabs(x->val);
  y->err = GSL_DBL_EPSILON * fabs(y->val);
  return status;
}


int
gsl_sf_rect_to_polar_impl(const double x, const double y,
                          gsl_sf_result * r, gsl_sf_result * theta)
{
  r->val = hypot(x, y);
  r->err = GSL_DBL_EPSILON * fabs(r->val);

  if(r->val > 0.0) {
    theta->val = atan2(y, x);
    theta->err = GSL_DBL_EPSILON * fabs(theta->val);
    return GSL_SUCCESS;
  }
  else {
    theta->val = 0.0;
    theta->err = 0.0;
    return GSL_EDOM;
  }
}


int gsl_sf_angle_restrict_symm_impl(double * theta)
{
  /* synthetic extended precision constants */
  const double P1 = 4 * 7.8539812564849853515625e-1;
  const double P2 = 4 * 3.7748947079307981766760e-8;
  const double P3 = 4 * 2.6951514290790594840552e-15;
  const double TwoPi = 2*(P1 + P2 + P3);

  const double t = *theta;
  const double y = 2*floor(t/TwoPi);
  double r = ((t - y*P1) - y*P2) - y*P3;

  if(r >  M_PI) r -= TwoPi;
  *theta = r;

  if(t > 0.0625/GSL_DBL_EPSILON)
    return GSL_ELOSS;
  else
    return GSL_SUCCESS;
}


int gsl_sf_angle_restrict_pos_impl(double * theta)
{
  /* synthetic extended precision constants */
  const double P1 = 4 * 7.85398125648498535156e-1;
  const double P2 = 4 * 3.77489470793079817668e-8;
  const double P3 = 4 * 2.69515142907905952645e-15;
  const double TwoPi = 2*(P1 + P2 + P3);

  const double t = *theta;
  const double y = 2*floor(t/TwoPi);

  *theta = ((t - y*P1) - y*P2) - y*P3;

  if(t > 0.0625/GSL_DBL_EPSILON)
    return GSL_ELOSS;
  else
    return GSL_SUCCESS;
}


int gsl_sf_sin_err_impl(const double x, const double dx, gsl_sf_result * result)
{
  double s = sin(x);
  result->val = s;
  result->err = fabs((1.0 - s*s) * dx) + GSL_DBL_EPSILON * fabs(s);
  return GSL_SUCCESS;
}


int gsl_sf_cos_err_impl(const double x, const double dx, gsl_sf_result * result)
{
  double c = cos(x);
  result->val = c;
  result->err = fabs((1.0 - c*c) * dx) + GSL_DBL_EPSILON * fabs(c);
  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_complex_sin_e(const double zr, const double zi, gsl_sf_result * szr, gsl_sf_result * szi)
{
  int status = gsl_sf_complex_sin_impl(zr, zi, szr, szi);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_sin_e", status);
  }
  return status;
}

int gsl_sf_complex_logsin_e(const double zr, const double zi, gsl_sf_result * lszr, gsl_sf_result * lszi)
{
  int status = gsl_sf_complex_logsin_impl(zr, zi, lszr, lszi);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_logsin_e", status);
  }
  return status;
}

int gsl_sf_complex_cos_e(const double zr, const double zi, gsl_sf_result * czr, gsl_sf_result * czi)
{
  int status = gsl_sf_complex_cos_impl(zr, zi, czr, czi);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_cos_e", status);
  }
  return status;
}

int gsl_sf_lnsinh_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_lnsinh_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_lnsinh_e", status);
  }
  return status;
}

int gsl_sf_lncosh_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_lncosh_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_lncosh_e", status);
  }
  return status;
}

int gsl_sf_polar_to_rect_e(const double r, const double theta, gsl_sf_result * x, gsl_sf_result * y)
{
  int status = gsl_sf_polar_to_rect_impl(r, theta, x, y);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_polar_to_rect_e", status);
  }
  return status;
}

int gsl_sf_rect_to_polar_e(const double x, const double y, gsl_sf_result * r, gsl_sf_result * theta)
{
  int status = gsl_sf_rect_to_polar_impl(x, y, r, theta);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_rect_to_polar_e", status);
  }
  return status;
}

int gsl_sf_angle_restrict_symm_e(double * theta)
{
  int status = gsl_sf_angle_restrict_symm_impl(theta);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_angle_restrict_symm_e", status);
  }
  return status;
}

int gsl_sf_angle_restrict_pos_e(double * theta)
{
  int status = gsl_sf_angle_restrict_pos_impl(theta);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_angle_restrict_pos_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes  *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_sin_pi_x(const double x)
{
  const double N = floor(x + 0.5);
  const double f = x - N;
  double result;

  if(N < INT_MAX && N > INT_MIN) {
    /* Make it an integer if we can. Saves another
     * call to floor().
     */
    const int intN    = (int)N;
    const double sign = ( GSL_IS_ODD(intN) ? -1.0 : 1.0 );
    result = sign * sin(M_PI * f);
  }
  else if(N > 2.0/GSL_DBL_EPSILON || N < -2.0/GSL_DBL_EPSILON) {
    /* All integer-valued floating point numbers
     * bigger than 2/eps=2^53 are actually even.
     */
    result = 0.0;
  }
  else {
    const double resN = N - 2.0*floor(0.5*N); /* 0 for even N, 1 for odd N */
    const double sign = ( fabs(resN) > 0.5 ? -1.0 : 1.0 );
    result = sign * sin(M_PI*f);
  }

  return result;
}
