/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_log.h"
#include "gsl_sf_trig.h"


/* sinh(x) series
 * double-precision for |x| < 1.0
 */
inline
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
inline
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


/* Chebyshev expansion for f(t) = sinc((t+1)/2), -1 < t < 1
 */
static double sinc_data[17] = {
  1.133648177811747875422,
 -0.532677564732557348781,
 -0.068293048346633177859,
  0.033403684226353715020,
  0.001485679893925747818,
 -0.000734421305768455295,
 -0.000016837282388837229,
  0.000008359950146618018,
  0.000000117382095601192,
 -0.000000058413665922724,
 -0.000000000554763755743,
  0.000000000276434190426,
  0.000000000001895374892,
 -0.000000000000945237101,
 -0.000000000000004900690,
  0.000000000000002445383,
  0.000000000000000009925
};
static gsl_sf_cheb_series sinc_cs = {
  sinc_data,
  16,
  -1, 1,
  (double *)0,
  (double *)0,
  10
};


/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (sin(x)/x - 1)/(x*x)
 */
static double sin_data[12] = {
  -0.3295190160663511504173,
   0.0025374284671667991990,
   0.0006261928782647355874,
  -4.6495547521854042157541e-06,
  -5.6917531549379706526677e-07,
   3.7283335140973803627866e-09,
   3.0267376484747473727186e-10,
  -1.7400875016436622322022e-12,
  -1.0554678305790849834462e-13,
   5.3701981409132410797062e-16,
   2.5984137983099020336115e-17,
  -1.1821555255364833468288e-19
};
static gsl_sf_cheb_series sin_cs = {
  sin_data,
  11,
  -1, 1,
  (double *)0,
  (double *)0,
  11
};

/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (2(cos(x) - 1)/(x^2) + 1) / x^2
 */
static double cos_data[11] = {
  0.165391825637921473505668118136,
 -0.00084852883845000173671196530195,
 -0.000210086507222940730213625768083,
  1.16582269619760204299639757584e-6,
  1.43319375856259870334412701165e-7,
 -7.4770883429007141617951330184e-10,
 -6.0969994944584252706997438007e-11,
  2.90748249201909353949854872638e-13,
  1.77126739876261435667156490461e-14,
 -7.6896421502815579078577263149e-17,
 -3.7363121133079412079201377318e-18
};
static gsl_sf_cheb_series cos_cs = {
  cos_data,
  10,
  -1, 1,
  (double *)0,
  (double *)0,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* I would have prefered just using the library sin() function.
 * But after some experimentation I decided that there was
 * no good way to understand the error; library sin() is just a black box.
 * So we have to roll our own.
 */
int
gsl_sf_sin_impl(double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else {
    const double P1 = 7.85398125648498535156e-1;
    const double P2 = 3.77489470793079817668e-8;
    const double P3 = 2.69515142907905952645e-15;

    const double sgn_x = GSL_SIGN(x);
    const double abs_x = fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const double x2 = x*x;
      result->val = x * (1.0 - x2/6.0);
      result->err = fabs(x*x2*x2 / 100.0);
      return GSL_SUCCESS;
    }
    else {
      double sgn_result = sgn_x;
      double y = floor(abs_x/(0.25*M_PI));
      int octant = y - ldexp(floor(ldexp(y,-3)),3);
      int stat_cs;
      double z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
	octant &= 07;
	y += 1.0;
      }

      if(octant > 3) {
        octant -= 4;
	sgn_result = -sgn_result;
      }
      
      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result sin_cs_result;
	const double t = 8.0*fabs(z)/M_PI - 1.0;
	stat_cs = gsl_sf_cheb_eval_impl(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result cos_cs_result;
	const double t = 8.0*fabs(z)/M_PI - 1.0;
	stat_cs = gsl_sf_cheb_eval_impl(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}


int
gsl_sf_cos_impl(double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else {
    const double P1 = 7.85398125648498535156e-1;
    const double P2 = 3.77489470793079817668e-8;
    const double P3 = 2.69515142907905952645e-15;

    const double abs_x = fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const double x2 = x*x;
      result->val = 1.0 - 0.5*x2;
      result->err = fabs(x2*x2/12.0);
      return GSL_SUCCESS;
    }
    else {
      double sgn_result = 1.0;
      double y = floor(abs_x/(0.25*M_PI));
      int octant = y - ldexp(floor(ldexp(y,-3)),3);
      int stat_cs;
      double z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
	octant &= 07;
	y += 1.0;
      }

      if(octant > 3) {
        octant -= 4;
	sgn_result = -sgn_result;
      }

      if(octant > 1) {
        sgn_result = -sgn_result;
      }

      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result cos_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = gsl_sf_cheb_eval_impl(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result sin_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = gsl_sf_cheb_eval_impl(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}


int
gsl_sf_hypot_impl(const double x, const double y, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x == 0.0 && y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    const double a = fabs(x);
    const double b = fabs(y);
    const double min = GSL_MIN_DBL(a,b);
    const double max = GSL_MAX_DBL(a,b);
    const double rat = min/max;
    const double root_term = sqrt(1.0 + rat*rat);

    if(max < DBL_MAX/root_term) {
      result->val = max * root_term;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      result->val = 0.0; /* FIXME: should be Inf */
      result->err = 0.0;
      return GSL_EOVRFLW;
    }
  }
}


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
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else if(zi < -60.0) {
    lszr->val = -M_LN2 - zi;
    lszi->val = -0.5*M_PI + zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
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
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = x + log(0.5*(1.0 - exp(-2.0*x)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
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
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
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
  double c = cos(t);
  double s = sin(t);
  x->val = r * cos(t);
  y->val = r * sin(t);
  x->err  = r * fabs(s * GSL_DBL_EPSILON * t);
  x->err += 2.0 * GSL_DBL_EPSILON * fabs(x->val);
  y->err  = r * fabs(c * GSL_DBL_EPSILON * t);
  y->err += 2.0 * GSL_DBL_EPSILON * fabs(y->val);
  return status;
}


int
gsl_sf_rect_to_polar_impl(const double x, const double y,
                          gsl_sf_result * r, gsl_sf_result * theta)
{
  int stat_h = gsl_sf_hypot_impl(x, y, r);
  if(r->val > 0.0) {
    theta->val = atan2(y, x);
    theta->err = 2.0 * GSL_DBL_EPSILON * fabs(theta->val);
    return stat_h;
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
  int stat_s = gsl_sf_sin_impl(x, result);
  result->err += fabs(cos(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_s;
}


int gsl_sf_cos_err_impl(const double x, const double dx, gsl_sf_result * result)
{
  int stat_c = gsl_sf_cos_impl(x, result);
  result->err += fabs(sin(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_c;
}


#if 0
int
gsl_sf_sin_pi_x_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(-100.0 < x && x < 100.0) {
    result->val = sin(M_PI * x) / (M_PI * x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double N = floor(x + 0.5);
    const double f = x - N;

    if(N < INT_MAX && N > INT_MIN) {
      /* Make it an integer if we can. Saves another
       * call to floor().
       */
      const int intN    = (int)N;
      const double sign = ( GSL_IS_ODD(intN) ? -1.0 : 1.0 );
      result->val = sign * sin(M_PI * f);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
    }
    else if(N > 2.0/GSL_DBL_EPSILON || N < -2.0/GSL_DBL_EPSILON) {
      /* All integer-valued floating point numbers
       * bigger than 2/eps=2^53 are actually even.
       */
      result->val = 0.0;
      result->err = 0.0;
    }
    else {
      const double resN = N - 2.0*floor(0.5*N); /* 0 for even N, 1 for odd N */
      const double sign = ( fabs(resN) > 0.5 ? -1.0 : 1.0 );
      result->val = sign * sin(M_PI*f);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
    }

    return GSL_SUCCESS;
  }
}
#endif


int gsl_sf_sinc_impl(double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else {
    const double ax = fabs(x);

    if(ax < 0.8) {
      /* Do not go to the limit of the fit since
       * there is a zero there and the Chebyshev
       * accuracy will go to zero.
       */
      return gsl_sf_cheb_eval_impl(&sinc_cs, 2.0*ax-1.0, result);
    }
    else if(ax < 100.0) {
      /* Small arguments are no problem.
       * We trust the library sin() to
       * roughly machine precision.
       */
      result->val = sin(M_PI * ax)/(M_PI * ax);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* Large arguments must be handled separately.
       */
      const double r = M_PI*ax;
      gsl_sf_result s;
      int stat_s = gsl_sf_sin_impl(r, &s);
      result->val = s.val/r;
      result->err = s.err/r + 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_s;
    }
  }
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/


int gsl_sf_sin_e(const double x, gsl_sf_result * r)
{
  int status = gsl_sf_sin_impl(x, r);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_sin_e", status);
  }
  return status;
}

int gsl_sf_cos_e(const double x, gsl_sf_result * r)
{
  int status = gsl_sf_cos_impl(x, r);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_cos_e", status);
  }
  return status;
}

int gsl_sf_hypot_e(const double x, const double y, gsl_sf_result * r)
{
  int status = gsl_sf_hypot_impl(x, y, r);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hypot_e", status);
  }
  return status;
}

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

#if 0
int gsl_sf_sin_pi_x_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_sin_pi_x_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_sin_pi_x_e", status);
  }
  return status;
}
#endif

int gsl_sf_sinc_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_sinc_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_sinc_e", status);
  }
  return status;
}
