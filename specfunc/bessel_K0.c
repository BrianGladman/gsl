/* specfunc/bessel_K0.c
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

/* based on SLATEC bk0(), bk0e() */

/* chebyshev expansions 

 series for bk0        on the interval  0.          to  4.00000d+00
                                        with weighted error   3.57e-19
                                         log weighted error  18.45
                               significant figures required  17.99
                                    decimal places required  18.97

 series for ak0        on the interval  1.25000d-01 to  5.00000d-01
                                        with weighted error   5.34e-17
                                         log weighted error  16.27
                               significant figures required  14.92
                                    decimal places required  16.89

 series for ak02       on the interval  0.          to  1.25000d-01
                                        with weighted error   2.34e-17
                                         log weighted error  16.63
                               significant figures required  14.67
                                    decimal places required  17.20
*/

/* from SLATEC dbesk0.f */
static double bk0_data[11] = {
  -.3532739323390276872E-1,
  +.3442898999246284869E+0,
  +.3597993651536150163E-1,
  +.1264615411446925923E-2,
  +.2286212103119451786E-4,
  +.2534791079026149457E-6,
  +.1904516377220208859E-8,
  +.1034969525763362459E-10,
  +.4259816142791082577E-13,
  +.1374465435880750897E-15,
  +.3570896528508373591E-18
};

static cheb_series bk0_cs = {
  bk0_data,
  10,
  -1, 1,
  10
};

/* from SLATEC dbsk0e.f */
static double ak0_data[17] = {
  -.7643947903327941424E-1,
  -.2235652605699819052E-1,
  +.7734181154693858235E-3,
  -.4281006688886099464E-4,
  +.3081700173862974744E-5,
  -.2639367222009664974E-6,
  +.2563713036403469206E-7,
  -.2742705549900201264E-8,
  +.3169429658097499592E-9,
  -.3902353286962184142E-10,
  +.5068040698188575402E-11,
  -.6889574741007870680E-12,
  +.9744978497825917691E-13,
  -.1427332841884548505E-13,
  +.2156412571021463040E-14,
  -.3349654255149562772E-15,
  +.5335260216952911692E-16
};

static cheb_series ak0_cs = {
  ak0_data,
  16,
  -1, 1,
  10
};

/* from SLATEC dbsk0e.f */
static double ak02_data[14] = {
  -.1201869826307592240E-1,
  -.9174852691025695311E-2,
  +.1444550931775005821E-3,
  -.4013614175435709729E-5,
  +.1567831810852310673E-6,
  -.7770110438521737710E-8,
  +.4611182576179717883E-9,
  -.3158592997860565771E-10,
  +.2435018039365041128E-11,
  -.2074331387398347898E-12,
  +.1925787280589917085E-13,
  -.1927554805838956104E-14,
  +.2062198029197818278E-15,
  -.2341685117579242403E-16
};

static cheb_series ak02_cs = {
  ak02_data,
  13,
  -1, 1,
  8
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_K0_scaled_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x <= 2.0) {
    const double lx = log(x);
    const double ex = exp(x);
    int stat_I0;
    gsl_sf_result I0;
    gsl_sf_result c;
    cheb_eval_e(&bk0_cs, 0.5*x*x-1.0, &c);
    stat_I0 = gsl_sf_bessel_I0_e(x, &I0);
    result->val  = ex * (-log(0.5*x)*I0.val - 0.25 + c.val);
    result->err  = ex * ((M_LN2+fabs(lx))*I0.err + c.err);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_I0;
  }
  else if(x <= 8.0) {
    const double sx = sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak0_cs, (16.0/x-5.0)/3.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = c.err / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double sx = sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak02_cs, 16.0/x-1.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = (c.err + GSL_DBL_EPSILON) / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  } 
}


int gsl_sf_bessel_K0_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x <= 2.0) {
    const double lx = log(x);
    int stat_I0;
    gsl_sf_result I0;
    gsl_sf_result c;
    cheb_eval_e(&bk0_cs, 0.5*x*x-1.0, &c);
    stat_I0 = gsl_sf_bessel_I0_e(x, &I0);
    result->val  = (-log(0.5*x))*I0.val - 0.25 + c.val;
    result->err  = (fabs(lx) + M_LN2) * I0.err + c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_I0;
  }
  else {
    gsl_sf_result K0_scaled;
    int stat_K0 = gsl_sf_bessel_K0_scaled_e(x, &K0_scaled);
    int stat_e  = gsl_sf_exp_mult_err_e(-x, GSL_DBL_EPSILON*fabs(x),
                                           K0_scaled.val, K0_scaled.err,
                                           result);
    return GSL_ERROR_SELECT_2(stat_e, stat_K0);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_bessel_K0_scaled(const double x)
{
  EVAL_RESULT(gsl_sf_bessel_K0_scaled_e(x, &result));
}

double gsl_sf_bessel_K0(const double x)
{
  EVAL_RESULT(gsl_sf_bessel_K0_e(x, &result));
}

