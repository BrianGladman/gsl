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
static double bk0_data[16] = {
 -.353273932339027687201140060063153E-1,
 .344289899924628486886344927529213E+0,
 .359799365153615016265721303687231E-1,
 .126461541144692592338479508673447E-2,
 .228621210311945178608269830297585E-4,
 .253479107902614945730790013428354E-6,
 .190451637722020885897214059381366E-8,
 .103496952576336245851008317853089E-10,
 .425981614279108257652445327170133E-13,
 .137446543588075089694238325440000E-15,
 .357089652850837359099688597333333E-18,
 .763164366011643737667498666666666E-21,
 .136542498844078185908053333333333E-23,
 .207527526690666808319999999999999E-26,
 .271281421807298560000000000000000E-29,
 .308259388791466666666666666666666E-32
};

static cheb_series bk0_cs = {
  bk0_data,
  15,
  -1, 1,
  10
};

/* from SLATEC dbsk0e.f */
static double ak0_data[38] = {
  -.7643947903327941424082978270088E-1,
  -.2235652605699819052023095550791E-1,
  +.7734181154693858235300618174047E-3,
  -.4281006688886099464452146435416E-4,
  +.3081700173862974743650014826660E-5,
  -.2639367222009664974067448892723E-6,
  +.2563713036403469206294088265742E-7,
  -.2742705549900201263857211915244E-8,
  +.3169429658097499592080832873403E-9,
  -.3902353286962184141601065717962E-10,
  +.5068040698188575402050092127286E-11,
  -.6889574741007870679541713557984E-12,
  +.9744978497825917691388201336831E-13,
  -.1427332841884548505389855340122E-13,
  +.2156412571021463039558062976527E-14,
  -.3349654255149562772188782058530E-15,
  +.5335260216952911692145280392601E-16,
  -.8693669980890753807639622378837E-17,
  +.1446404347862212227887763442346E-17,
  -.2452889825500129682404678751573E-18,
  +.4233754526232171572821706342400E-19,
  -.7427946526454464195695341294933E-20,
  +.1323150529392666866277967462400E-20,
  -.2390587164739649451335981465599E-21,
  +.4376827585923226140165712554666E-22,
  -.8113700607345118059339011413333E-23,
  +.1521819913832172958310378154666E-23,
  -.2886041941483397770235958613333E-24,
  +.5530620667054717979992610133333E-25,
  -.1070377329249898728591633066666E-25,
  +.2091086893142384300296328533333E-26,
  -.4121713723646203827410261333333E-27,
  +.8193483971121307640135680000000E-28,
  -.1642000275459297726780757333333E-28,
  +.3316143281480227195890346666666E-29,
  -.6746863644145295941085866666666E-30,
  +.1382429146318424677635413333333E-30,
  -.2851874167359832570811733333333E-31
};

static cheb_series ak0_cs = {
  ak0_data,
  37,
  -1, 1,
  16
};

/* from SLATEC dbsk0e.f */
static double ak02_data[33] = {
  -.1201869826307592239839346212452E-1,
  -.9174852691025695310652561075713E-2,
  +.1444550931775005821048843878057E-3,
  -.4013614175435709728671021077879E-5,
  +.1567831810852310672590348990333E-6,
  -.7770110438521737710315799754460E-8,
  +.4611182576179717882533130529586E-9,
  -.3158592997860565770526665803309E-10,
  +.2435018039365041127835887814329E-11,
  -.2074331387398347897709853373506E-12,
  +.1925787280589917084742736504693E-13,
  -.1927554805838956103600347182218E-14,
  +.2062198029197818278285237869644E-15,
  -.2341685117579242402603640195071E-16,
  +.2805902810643042246815178828458E-17,
  -.3530507631161807945815482463573E-18,
  +.4645295422935108267424216337066E-19,
  -.6368625941344266473922053461333E-20,
  +.9069521310986515567622348800000E-21,
  -.1337974785423690739845005311999E-21,
  +.2039836021859952315522088960000E-22,
  -.3207027481367840500060869973333E-23,
  +.5189744413662309963626359466666E-24,
  -.8629501497540572192964607999999E-25,
  +.1472161183102559855208038400000E-25,
  -.2573069023867011283812351999999E-26,
  +.4601774086643516587376640000000E-27,
  -.8411555324201093737130666666666E-28,
  +.1569806306635368939301546666666E-28,
  -.2988226453005757788979199999999E-29,
  +.5796831375216836520618666666666E-30,
  -.1145035994347681332155733333333E-30,
  +.2301266594249682802005333333333E-31
};

static cheb_series ak02_cs = {
  ak02_data,
  32,
  -1, 1,
  13
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

