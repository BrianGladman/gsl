/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_amp_phase.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_trig.h"
#include "gsl_sf_bessel.h"

#define ROOT_EIGHT (2.0*M_SQRT2)


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* based on SLATEC besj1, 1983 version, w. fullerton */

/* chebyshev expansions

 series for bj1        on the interval  0.	    to  1.60000d+01
					with weighted error   4.48e-17
					 log weighted error  16.35
			       significant figures required  15.77
				    decimal places required  16.89

*/
static double bj1_data[12] = {
  -0.11726141513332787,
  -0.25361521830790640,
   0.050127080984469569,
  -0.004631514809625081,
   0.000247996229415914,
  -0.000008678948686278,
   0.000000214293917143,
  -0.000000003936093079,
   0.000000000055911823,
  -0.000000000000632761,
   0.000000000000005840,
  -0.000000000000000044,
};
static gsl_sf_cheb_series bj1_cs = {
  bj1_data,
  11,
  -1, 1,
  (double *)0,
  (double *)0,
  8
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_J1_impl(const double x, gsl_sf_result * result)
{
  double y = fabs(x);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(y < 2.0*DBL_MIN) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(y < ROOT_EIGHT * GSL_SQRT_DBL_EPSILON) {
    result->val = 0.5*x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(y < 4.0) {
    gsl_sf_result c;
    gsl_sf_cheb_eval_impl(&bj1_cs, 0.125*y*y-1.0, &c);
    result->val = x * (0.25 + c.val);
    result->err = fabs(x * c.err);
    return GSL_SUCCESS;
  }
  else {
    /* Because the leading term in the phase is y,
     * which we assume is exactly known, the error
     * in the cos() evaluation is bounded.
     */
    const double z  = 32.0/(y*y) - 1.0;
    gsl_sf_result ca;
    gsl_sf_result ct;
    gsl_sf_result sp;
    const int stat_ca = gsl_sf_cheb_eval_impl(&_bessel_amp_phase_bm1_cs,  z, &ca);
    const int stat_ct = gsl_sf_cheb_eval_impl(&_bessel_amp_phase_bth1_cs, z, &ct);
    const int stat_sp = gsl_sf_bessel_sin_pi4_impl(y, ct.val/y, &sp);
    const double sqrty = sqrt(y);
    const double ampl  = (0.75 + ca.val) / sqrty;
    result->val  = (x < 0.0 ? -ampl : ampl) * sp.val;
    result->err  = fabs(sp.val) * ca.err/sqrty + fabs(ampl) * sp.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_J1_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_J1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_J1_e", status);
  }
  return status;
}
