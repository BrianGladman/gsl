/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_log.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* Chebyshev expansion for log(2 + t/2), -1<t<1
 */
static double lopx_data[30] = {
  1.35428541555133985578563385941,
  0.254033307585166229641469200435,
 -0.0161332303406649185658768017408,
  0.00136612595515748913751987580588,
 -0.000130140560612475542180853964897,
  0.0000132240148253499874455398918079,
 -1.39972509401622151636330212152e-6,
  1.52390055146956508248434684778e-7,
 -1.69365655165294204665662023094e-8,
  1.91220078102081690235490603959e-9,
 -2.18593210126345323702361637506e-10,
  2.52408891745703431813433312857e-11,
 -2.93884550822808024454973150696e-12,
  3.4456829766322526784596773276e-13,
 -4.0639775588461780098847580274e-14,
  4.8177997523860592715565643205e-15,
 -5.7369450299138695018259966392e-16,
  6.8582358653331797060252480814e-17,
 -8.2271516106073998484857093291e-18,
  9.8998604320131727145526155564e-19,
 -1.19457478783351399422849426670e-19,
  1.44505611671997976193304147635e-20,
 -1.75203183738809143593277940304e-21,
  2.12861689961232292096976734884e-22,
 -2.59104387636970898051071813247e-23,
  3.15941494085994097400455648019e-24,
 -3.8586376320223228762842093164e-25,
  4.7195897206140425059171911128e-26,
 -5.7805697595818690872209693909e-27,
  7.0891039932921375654633712594e-28
};
static struct gsl_sf_cheb_series lopx_cs = {
  lopx_data,
  29,
  -1, 1,
  (double *)0,
  (double *)0
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_complex_log_impl(const double zr, const double zi, double * lnr, double * theta)
{
  if(zr != 0.0 || zi != 0.0) {
    double r2 = zr*zr + zi*zi;
    *lnr = 0.5*log(r2);
    *theta = atan2(zi, zr);
    return GSL_SUCCESS;
  }
  else {
    return GSL_EDOM;
  }
}

int gsl_sf_log_1plusx_impl(const double x, double * result)
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(fabs(x) <  0.5) {
    double t = 2.0*(x-1.0);
    *result = gsl_sf_cheb_eval(&lopx_cs, t);
    return GSL_SUCCESS;
  }
  else {
    *result = log(1.0 + x);
    return GSL_SUCCESS;
  }
}

/*-*-*-*-*-*-*-*-*-*-*-* Error Handling Versions *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_complex_log_e(const double zr, const double zi, double * lnr, double * theta)
{
  int status = gsl_sf_complex_log_impl(zr, zi, lnr, theta);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_log_e", status);
  }
  return status;
}
