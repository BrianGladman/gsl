/* specfunc/bessel_K1.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * Copyright (C) 2016 Pavel Holoborodko, Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_bessel.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besk1(), besk1e() */

/* chebyshev expansions 

 series for bk1        on the interval  0.          to  4.00000d+00
                                        with weighted error   7.02e-18
                                         log weighted error  17.15
                               significant figures required  16.73
                                    decimal places required  17.67

 series for ak1        on the interval  1.25000d-01 to  5.00000d-01
                                        with weighted error   6.06e-17
                                         log weighted error  16.22
                               significant figures required  15.41
                                    decimal places required  16.83

 series for ak12       on the interval  0.          to  1.25000d-01
                                        with weighted error   2.58e-17
                                         log weighted error  16.59
                               significant figures required  15.22
                                    decimal places required  17.16
*/

/* from SLATEC dbesk1.f */
static double bk1_data[11] = {
  +.2530022733894777053E-1,
  -.3531559607765448757E+0,
  -.1226111808226571482E+0,
  -.6975723859639864350E-2,
  -.1730288957513052063E-3,
  -.2433406141565968235E-5,
  -.2213387630734725856E-7,
  -.1411488392633527761E-9,
  -.6666901694199329006E-12,
  -.2427449850519365934E-14,
  -.7023863479386287597E-17
};

static cheb_series bk1_cs = {
  bk1_data,
  10,
  -1, 1,
  8
};

static double ak1_data[25] = {
  +2.07996868001418246152e-01,
  +1.62581565017881475959e-01,
  -5.87070423518863640171e-03,
  +4.95021520115789501462e-04,
  -5.78958347598556986444e-05,
  +8.18614610209334725538e-06,
  -1.31604832009487277149e-06,
  +2.32546031520101213185e-07,
  -4.42206518311557986572e-08,
  +8.92163994883100360555e-09,
  -1.89046270526983427490e-09,
  +4.17568808108504701844e-10,
  -9.55912361791375793677e-11,
  +2.25769353153867757611e-11,
  -5.48128000211158482030e-12,
  +1.36386122546441925778e-12,
  -3.46936690565986409232e-13,
  +9.00354564415705941638e-14,
  -2.37950577776254431920e-14,
  +6.39447503964025335916e-15,
  -1.74498363492322044105e-15,
  +4.82994547989290473028e-16,
  -1.35460927805445605587e-16,
  +3.84604274446777234121e-17,
  -1.10456856122581316056e-17
};

static cheb_series ak1_cs = {
  ak1_data,
  24,
  -1, 1,
  9
};

/* from SLATEC dbsk1e.f */
static double ak12_data[14] = {
  +.637930834373900104E-1,
  +.283288781304972094E-1,
  -.247537067390525035E-3,
  +.577197245160724882E-5,
  -.206893921953654830E-6,
  +.973998344138180418E-8,
  -.558533614038062498E-9,
  +.373299663404618524E-10,
  -.282505196102322545E-11,
  +.237201900248414417E-12,
  -.217667738799175398E-13,
  +.215791416161603245E-14,
  -.229019693071826928E-15,
  +.258288572982327496E-16
};

static cheb_series ak12_cs = {
  ak12_data,
  13,
  -1, 1,
  7
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

/*
gsl_sf_bessel_K1_scaled_e()
  Computed scaled K1 Bessel function

Notes:
1) On [0,1] a Chebyshev expansion from SLATEC is used

2) On [1,8] a new Chebyshev series from Pavel Holoborodko is used

3) On [8,inf] a Chebyshev series from SLATEC is used
*/

int
gsl_sf_bessel_K1_scaled_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0*GSL_DBL_MIN) {
    OVERFLOW_ERROR(result);
  }
  else if(x <= 1.0) {
    const double lx = log(x);
    const double ex = exp(x);
    int stat_I1;
    gsl_sf_result I1;
    gsl_sf_result c;
    cheb_eval_e(&bk1_cs, 0.5*x*x-1.0, &c);
    stat_I1 = gsl_sf_bessel_I1_e(x, &I1);
    result->val  = ex * (log(0.5*x)*I1.val + (0.75 + c.val)/x);
    result->err  = ex * (c.err/x + fabs(lx)*I1.err);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_I1;
  }
  else if(x <= 8.0) {
    const double sx = sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak1_cs, (16.0/x-9.0)/7.0, &c);
    result->val  = (1.375 + c.val) / sx;
    result->err  = c.err / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double sx = sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak12_cs, 16.0/x-1.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = c.err / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_K1_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0*GSL_DBL_MIN) {
    OVERFLOW_ERROR(result);
  }
  else if(x <= 2.0) {
    const double lx = log(x);
    int stat_I1;
    gsl_sf_result I1;
    gsl_sf_result c;
    cheb_eval_e(&bk1_cs, 0.5*x*x-1.0, &c);
    stat_I1 = gsl_sf_bessel_I1_e(x, &I1);
    result->val  = log(0.5*x)*I1.val + (0.75 + c.val)/x;
    result->err  = c.err/x + fabs(lx)*I1.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_I1;
  }
  else {
    gsl_sf_result K1_scaled;
    int stat_K1 = gsl_sf_bessel_K1_scaled_e(x, &K1_scaled);
    int stat_e  = gsl_sf_exp_mult_err_e(-x, 0.0,
                                           K1_scaled.val, K1_scaled.err,
                                           result);
    result->err = fabs(result->val) * (GSL_DBL_EPSILON*fabs(x) + K1_scaled.err/K1_scaled.val);
    return GSL_ERROR_SELECT_2(stat_e, stat_K1);
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_bessel_K1_scaled(const double x)
{
  EVAL_RESULT(gsl_sf_bessel_K1_scaled_e(x, &result));
}

double gsl_sf_bessel_K1(const double x)
{
  EVAL_RESULT(gsl_sf_bessel_K1_e(x, &result));
}
