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


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* based on SLATEC besj0, 1977 version, w. fullerton */

/* chebyshev expansions for Bessel functions

 series for bj0        on the interval  0.	    to  1.60000d+01
					with weighted error   7.47e-18
					 log weighted error  17.13
			       significant figures required  16.98
				    decimal places required  17.68

*/

static double bj0_data[13] = {
   0.100254161968939137, 
  -0.665223007764405132, 
   0.248983703498281314, 
  -0.0332527231700357697,
   0.0023114179304694015,
  -0.0000991127741995080,
   0.0000028916708643998,
  -0.0000000612108586630,
   0.0000000009838650793,
  -0.0000000000124235515,
   0.0000000000001265433,
  -0.0000000000000010619,
   0.0000000000000000074,
};
static gsl_sf_cheb_series bj0_cs = {
  bj0_data,
  12,
  -1, 1,
  (double *)0,
  (double *)0,
  9
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_J0_impl(const double x, gsl_sf_result * result)
{
  double y = fabs(x);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(y < 2.0*GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0;
    result->err = y*y;
    return GSL_SUCCESS;
  }
  else if(y <= 4.0) {
    return gsl_sf_cheb_eval_impl(&bj0_cs, 0.125*y*y - 1.0, result);
  }
  else {
    const double z = 32.0/(y*y) - 1.0;
    gsl_sf_result ca;
    gsl_sf_result ct;
    gsl_sf_result cp;
    const int stat_ca = gsl_sf_cheb_eval_impl(&_bessel_amp_phase_bm0_cs,  z, &ca);
    const int stat_ct = gsl_sf_cheb_eval_impl(&_bessel_amp_phase_bth0_cs, z, &ct);
    const int stat_cp = gsl_sf_bessel_cos_pi4_impl(y, ct.val/y, &cp);
    const double sqrty = sqrt(y);
    const double ampl  = (0.75 + ca.val) / sqrty;
    result->val  = ampl * cp.val;
    result->err  = fabs(cp.val) * ca.err/sqrty + fabs(ampl) * cp.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_3(stat_ca, stat_ct, stat_cp);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_J0_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_J0_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_J0_e", status);
  }
  return status;
}
