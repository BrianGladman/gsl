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

/* based on SLATEC besy0, 1980 version, w. fullerton */

/* chebyshev expansions

 series for by0        on the interval  0.	    to  1.60000d+01
					with weighted error   1.20e-17
					 log weighted error  16.92
			       significant figures required  16.15
				    decimal places required  17.48
*/

static double by0_data[13] = {
  -0.011277839392865573,
  -0.128345237560420350,
  -0.104378847997942490,
   0.023662749183969695,
  -0.002090391647700486,
   0.000103975453939057,
  -0.000003369747162423,
   0.000000077293842676,
  -0.000000001324976772,
   0.000000000017648232,
  -0.000000000000188105,
   0.000000000000001641,
  -0.000000000000000011
};
static gsl_sf_cheb_series by0_cs = {
  by0_data,
  12,
  -1, 1,
  (double *)0,
  (double *)0,
  8
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Y0_impl(const double x, gsl_sf_result * result)
{
  const double two_over_pi = 2.0/M_PI;
  const double xmax        = 1.0/GSL_DBL_EPSILON;

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if (x <= 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x < 4.0) {
    gsl_sf_result J0;
    gsl_sf_result c;
    int stat_J0 = gsl_sf_bessel_J0_impl(x, &J0);
    gsl_sf_cheb_eval_impl(&by0_cs, 0.125*x*x-1.0, &c);
    result->val = two_over_pi*(-M_LN2 + log(x))*J0.val + 0.375 + c.val;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + c.err;
    return stat_J0;
  }
  else if(x < xmax) {
    /* Leading behaviour of phase is x, which is exact,
     * so the error is bounded.
     */
    const double z  = 32.0/(x*x) - 1.0;
    gsl_sf_result c1;
    gsl_sf_result c2;
    gsl_sf_result sp;
    const int stat_c1 = gsl_sf_cheb_eval_impl(&_gsl_sf_bessel_amp_phase_bm0_cs,  z, &c1);
    const int stat_c2 = gsl_sf_cheb_eval_impl(&_gsl_sf_bessel_amp_phase_bth0_cs, z, &c2);
    const int stat_sp = gsl_sf_bessel_sin_pi4_impl(x, c2.val/x, &sp);
    const double sqrtx = sqrt(x);
    const double ampl  = (0.75 + c1.val) / sqrtx;
    result->val  = ampl * sp.val;
    result->err  = fabs(sp.val) * c1.err/sqrtx + fabs(ampl) * sp.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_3(stat_sp, stat_c1, stat_c2);
  }
  else {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EUNDRFLW;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Y0_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_Y0_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Y0_e", status);
  }
  return status;
}
