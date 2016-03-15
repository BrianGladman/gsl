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
static double bk1_data[16] = {
  +.25300227338947770532531120868533E-1,
  -.35315596077654487566723831691801E+0,
  -.12261118082265714823479067930042E+0,
  -.69757238596398643501812920296083E-2,
  -.17302889575130520630176507368979E-3,
  -.24334061415659682349600735030164E-5,
  -.22133876307347258558315252545126E-7,
  -.14114883926335277610958330212608E-9,
  -.66669016941993290060853751264373E-12,
  -.24274498505193659339263196864853E-14,
  -.70238634793862875971783797120000E-17,
  -.16543275155100994675491029333333E-19,
  -.32338347459944491991893333333333E-22,
  -.53312750529265274999466666666666E-25,
  -.75130407162157226666666666666666E-28,
  -.91550857176541866666666666666666E-31
};

static cheb_series bk1_cs = {
  bk1_data,
  15,
  -1, 1,
  8
};

/* from SLATEC dbsk1e.f */
static double ak1_data[38] = {
  +.27443134069738829695257666227266E+0,
  +.75719899531993678170892378149290E-1,
  -.14410515564754061229853116175625E-2,
  +.66501169551257479394251385477036E-4,
  -.43699847095201407660580845089167E-5,
  +.35402774997630526799417139008534E-6,
  -.33111637792932920208982688245704E-7,
  +.34459775819010534532311499770992E-8,
  -.38989323474754271048981937492758E-9,
  +.47208197504658356400947449339005E-10,
  -.60478356628753562345373591562890E-11,
  +.81284948748658747888193837985663E-12,
  -.11386945747147891428923915951042E-12,
  +.16540358408462282325972948205090E-13,
  -.24809025677068848221516010440533E-14,
  +.38292378907024096948429227299157E-15,
  -.60647341040012418187768210377386E-16,
  +.98324256232648616038194004650666E-17,
  -.16284168738284380035666620115626E-17,
  +.27501536496752623718284120337066E-18,
  -.47289666463953250924281069568000E-19,
  +.82681500028109932722392050346666E-20,
  -.14681405136624956337193964885333E-20,
  +.26447639269208245978085894826666E-21,
  -.48290157564856387897969868800000E-22,
  +.89293020743610130180656332799999E-23,
  -.16708397168972517176997751466666E-23,
  +.31616456034040694931368618666666E-24,
  -.60462055312274989106506410666666E-25,
  +.11678798942042732700718421333333E-25,
  -.22773741582653996232867840000000E-26,
  +.44811097300773675795305813333333E-27,
  -.88932884769020194062336000000000E-28,
  +.17794680018850275131392000000000E-28,
  -.35884555967329095821994666666666E-29,
  +.72906290492694257991679999999999E-30,
  -.14918449845546227073024000000000E-30,
  +.30736573872934276300799999999999E-31
};

static cheb_series ak1_cs = {
  ak1_data,
  37,
  -1, 1,
  16
};

/* from SLATEC dbsk1e.f */
static double ak12_data[33] = {
  +.6379308343739001036600488534102E-1,
  +.2832887813049720935835030284708E-1,
  -.2475370673905250345414545566732E-3,
  +.5771972451607248820470976625763E-5,
  -.2068939219536548302745533196552E-6,
  +.9739983441381804180309213097887E-8,
  -.5585336140380624984688895511129E-9,
  +.3732996634046185240221212854731E-10,
  -.2825051961023225445135065754928E-11,
  +.2372019002484144173643496955486E-12,
  -.2176677387991753979268301667938E-13,
  +.2157914161616032453939562689706E-14,
  -.2290196930718269275991551338154E-15,
  +.2582885729823274961919939565226E-16,
  -.3076752641268463187621098173440E-17,
  +.3851487721280491597094896844799E-18,
  -.5044794897641528977117282508800E-19,
  +.6888673850418544237018292223999E-20,
  -.9775041541950118303002132480000E-21,
  +.1437416218523836461001659733333E-21,
  -.2185059497344347373499733333333E-22,
  +.3426245621809220631645388800000E-23,
  -.5531064394246408232501248000000E-24,
  +.9176601505685995403782826666666E-25,
  -.1562287203618024911448746666666E-25,
  +.2725419375484333132349439999999E-26,
  -.4865674910074827992378026666666E-27,
  +.8879388552723502587357866666666E-28,
  -.1654585918039257548936533333333E-28,
  +.3145111321357848674303999999999E-29,
  -.6092998312193127612416000000000E-30,
  +.1202021939369815834623999999999E-30,
  -.2412930801459408841386666666666E-31
};

static cheb_series ak12_cs = {
  ak12_data,
  32,
  -1, 1,
  13
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_K1_scaled_e(const double x, gsl_sf_result * result)
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
    cheb_eval_e(&ak1_cs, (16.0/x-5.0)/3.0, &c);
    result->val  = (1.25 + c.val) / sx;
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
