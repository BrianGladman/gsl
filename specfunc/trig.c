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

int gsl_sf_complex_sin_impl(const double zr, const double zi, double * szr, double * szi)
{
  if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    *szr = sin(zr)*(ch_m1 + 1.0);
    *szi = cos(zr)*sh;
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    *szr = sin(zr)*ch;
    *szi = cos(zr)*sh;
    return GSL_SUCCESS;
  }
  else {
    *szr = 0.0;
    *szi = 0.0;
    return GSL_EOVRFLW;
  }
}


int gsl_sf_complex_cos_impl(const double zr, const double zi, double * czr, double * czi)
{
  if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    *czr =  cos(zr)*(ch_m1 + 1.0);
    *czi = -sin(zr)*sh;
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    *czr =  cos(zr)*ch;
    *czi = -sin(zr)*sh;
    return GSL_SUCCESS;
  }
  else {
    *czr = 0.0;
    *czi = 0.0;
    return GSL_EOVRFLW;
  }
}


int gsl_sf_complex_logsin_impl(const double zr, const double zi, double * lszr, double * lszi)
{
  if(zi > 60.0) {
    *lszr = -M_LN2 + zi;
    *lszi =  0.5*M_PI - zr;
  }
  else if(zi < -60.0) {
    *lszr = -M_LN2 - zi;
    *lszi = -0.5*M_PI + zr;
  }
  else {
    double sin_r, sin_i;
    int status;
    gsl_sf_complex_sin_impl(zr, zi, &sin_r, &sin_i); /* ok by construction */
    status = gsl_sf_complex_log_impl(sin_r, sin_i, lszr, lszi);
    if(status == GSL_EDOM) {
      *lszr = 0.0;
      *lszi = 0.0;
      return GSL_EDOM;
    }
  }
  return gsl_sf_angle_restrict_symm_impl(lszi);
}


int gsl_sf_lnsinh_impl(const double x, double * result)
{
  if(x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(fabs(x) < 1.0) {
    double eps;
    sinh_series(x, &eps);
    *result = log(eps);
    return GSL_SUCCESS;
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    *result = x + log(0.5*(1.0 - exp(-2.0*x)));
    return GSL_SUCCESS;
  }
  else {
    *result = -M_LN2 + x;
    return GSL_SUCCESS;
  }
}


int gsl_sf_lncosh_impl(const double x, double * result)
{
  if(fabs(x) < 1.0) {
    double eps;
    cosh_m1_series(x, &eps);
    return gsl_sf_log_1plusx_impl(eps, result);
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    *result = x + log(0.5*(1.0 + exp(-2.0*x)));
    return GSL_SUCCESS;
  }
  else {
    *result = -M_LN2 + x;
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

int gsl_sf_polar_to_rect_impl(const double r, const double theta, double * x, double * y)
{
  double t   = theta;
  int status = gsl_sf_angle_restrict_symm_impl(&t);
  *x = r * cos(t);
  *y = r * sin(t);
  return status;
}


int gsl_sf_rect_to_polar_impl(const double x, const double y, double * r, double * theta)
{
  *r = hypot(x, y);
  if(*r > 0.0) {
    *theta = atan2(y, x);
    return GSL_SUCCESS;
  }
  else {
    *theta = 0.0;
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


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_complex_sin_e(const double zr, const double zi, double * szr, double * szi)
{
  int status = gsl_sf_complex_sin_impl(zr, zi, szr, szi);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_sin_e", status);
  }
  return status;
}

int gsl_sf_complex_logsin_e(const double zr, const double zi, double * lszr, double * lszi)
{
  int status = gsl_sf_complex_logsin_impl(zr, zi, lszr, lszi);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_logsin_e", status);
  }
  return status;
}

int gsl_sf_complex_cos_e(const double zr, const double zi, double * czr, double * czi)
{
  int status = gsl_sf_complex_cos_impl(zr, zi, czr, czi);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_cos_e", status);
  }
  return status;
}

int gsl_sf_lnsinh_e(const double x, double * result)
{
  int status = gsl_sf_lnsinh_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_lnsinh_e", status);
  }
  return status;
}

int gsl_sf_lncosh_e(const double x, double * result)
{
  int status = gsl_sf_lncosh_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_lncosh_e", status);
  }
  return status;
}

int gsl_sf_polar_to_rect_e(const double r, const double theta, double * x, double * y)
{
  int status = gsl_sf_polar_to_rect_impl(r, theta, x, y);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_polar_to_rect_e", status);
  }
  return status;
}

int gsl_sf_rect_to_polar_e(const double x, const double y, double * r, double * theta)
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

double gsl_sf_lnsinh(const double x)
{
  double y;
  int status = gsl_sf_lnsinh_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_lnsinh", status);
  }
  return y;
}

double gsl_sf_lncosh(const double x)
{
  double y;
  int status = gsl_sf_lncosh_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_lncosh", status);
  }
  return y;
}
