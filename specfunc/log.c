/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_log.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* Chebyshev expansion for log(1 + x(t))/x(t)
 *
 * x(t) = (4t-1)/(2(4-t))
 * t(x) = (8x+1)/(2(x+2))
 * -1/2 < x < 1/2
 * -1 < t < 1
 */
static double lopx_data[21] = {
  2.16647910664395270521272590407,
 -0.28565398551049742084877469679,
  0.01517767255690553732382488171,
 -0.00200215904941415466274422081,
  0.00019211375164056698287947962,
 -0.00002553258886105542567601400,
  2.9004512660400621301999384544e-06,
 -3.8873813517057343800270917900e-07,
  4.7743678729400456026672697926e-08,
 -6.4501969776090319441714445454e-09,
  8.2751976628812389601561347296e-10,
 -1.1260499376492049411710290413e-10,
  1.4844576692270934446023686322e-11,
 -2.0328515972462118942821556033e-12,
  2.7291231220549214896095654769e-13,
 -3.7581977830387938294437434651e-14,
  5.1107345870861673561462339876e-15,
 -7.0722150011433276578323272272e-16,
  9.7089758328248469219003866867e-17,
 -1.3492637457521938883731579510e-17,
  1.8657327910677296608121390705e-18
};
static gsl_sf_cheb_series lopx_cs = {
  lopx_data,
  20,
  -1, 1,
  (double *)0,
  (double *)0,
  10
};

/* Chebyshev expansion for (log(1 + x(t)) - x(t))/x(t)^2
 *
 * x(t) = (4t-1)/(2(4-t))
 * t(x) = (8x+1)/(2(x+2))
 * -1/2 < x < 1/2
 * -1 < t < 1
 */
static double lopxmx_data[20] = {
 -1.12100231323744103373737274541,
  0.19553462773379386241549597019,
 -0.01467470453808083971825344956,
  0.00166678250474365477643629067,
 -0.00018543356147700369785746902,
  0.00002280154021771635036301071,
 -2.8031253116633521699214134172e-06,
  3.5936568872522162983669541401e-07,
 -4.6241857041062060284381167925e-08,
  6.0822637459403991012451054971e-09,
 -8.0339824424815790302621320732e-10,
  1.0751718277499375044851551587e-10,
 -1.4445310914224613448759230882e-11,
  1.9573912180610336168921438426e-12,
 -2.6614436796793061741564104510e-13,
  3.6402634315269586532158344584e-14,
 -4.9937495922755006545809120531e-15,
  6.8802890218846809524646902703e-16,
 -9.5034129794804273611403251480e-17,
  1.3170135013050997157326965813e-17
};
static gsl_sf_cheb_series lopxmx_cs = {
  lopxmx_data,
  19,
  -1, 1,
  (double *)0,
  (double *)0,
  9
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_log_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else {
    result->val = log(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_log_abs_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else {
    result->val = log(fabs(x));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_complex_log_impl(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta)
{
  if(lnr == 0 || theta == 0) {
    return GSL_EFAULT;
  }
  else if(zr != 0.0 || zi != 0.0) {
    const double ax = fabs(zr);
    const double ay = fabs(zi);
    const double min = GSL_MIN(ax, ay);
    const double max = GSL_MAX(ax, ay);
    lnr->val = log(max) + 0.5 * log(1.0 + (min/max)*(min/max));
    lnr->err = 2.0 * GSL_DBL_EPSILON * fabs(lnr->val);
    theta->val = atan2(zi, zr);
    theta->err = GSL_DBL_EPSILON * fabs(lnr->val);
    return GSL_SUCCESS;
  }
  else {
    lnr->val = 0.0;
    lnr->err = 0.0;
    theta->val = 0.0;
    theta->err = 0.0;
    return GSL_EDOM;
  }
}


int
gsl_sf_log_1plusx_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= -1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(fabs(x) < GSL_ROOT6_DBL_EPSILON) {
    const double c1 = -0.5;
    const double c2 =  1.0/3.0;
    const double c3 = -1.0/4.0;
    const double c4 =  1.0/5.0;
    const double c5 = -1.0/6.0;
    const double c6 = -1.0/7.0;
    result->val = x * (1.0 + x*(c1 + x*(c2 + x*(c3 + x*(c4 + x*(c5 + x*c6))))));
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(fabs(x) < 0.5) {
    double t = 0.5*(8.0*x + 1.0)/(x+2.0);
    gsl_sf_result c;
    gsl_sf_cheb_eval_impl(&lopx_cs, t, &c);
    result->val = x * c.val;
    result->err = fabs(x * c.err);
    return GSL_SUCCESS;
  }
  else {
    result->val = log(1.0 + x);
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_log_1plusx_mx_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= -1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(fabs(x) < GSL_ROOT5_DBL_EPSILON) {
    const double c1 = -0.5;
    const double c2 =  1.0/3.0;
    const double c3 = -1.0/4.0;
    const double c4 =  1.0/5.0;
    const double c5 = -1.0/6.0;
    const double c6 = -1.0/7.0;
    result->val = x*x * (c1 + x*(c2 + x*(c3 + x*(c4 + x*(c5 + x*c6)))));
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(fabs(x) < 0.5) {
    double t = 0.5*(8.0*x + 1.0)/(x+2.0);
    gsl_sf_result c;
    gsl_sf_cheb_eval_impl(&lopxmx_cs, t, &c);
    result->val = x*x * c.val;
    result->err = x*x * c.err;
    return GSL_SUCCESS;
  }
  else {
    const double lterm = log(1.0 + x);
    result->val = lterm - x;
    result->err = GSL_DBL_EPSILON * (fabs(lterm) + fabs(x));
    return GSL_SUCCESS;
  }
}



/*-*-*-*-*-*-*-*-*-*-*-* Error Handling Versions *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_log_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_log_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_log_e", status);
  }
  return status;
}


int
gsl_sf_log_abs_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_log_abs_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_log_abs_e", status);
  }
  return status;
}


int
gsl_sf_complex_log_e(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta)
{
  int status = gsl_sf_complex_log_impl(zr, zi, lnr, theta);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_log_e", status);
  }
  return status;
}


int
gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_log_1plusx_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_log_1plusx_e", status);
  }
  return status;
}


int
gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_log_1plusx_mx_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_log_1plusx_mx_e", status);
  }
  return status;
}
