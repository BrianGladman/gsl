/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_test.h>
#include <gsl_sf.h>

#include "test_sf.h"


double test_sf_frac_diff(double x1, double x2)
{
  if(x1 == 0.0 && x2 == 0.0) {
    return 0.0;
  }
  else if(x1 <= DBL_MAX && x2 <= DBL_MAX && (x1 + x2 != 0.0))
    return fabs((x1-x2)/(x1+x2));
  else
    return 1.0;
}


int
test_sf_analyze(const char * message, gsl_sf_result r,
                double val, double tol, int return_code,
		int expect_return_code)
{
  unsigned int s = 0;
  double f = test_sf_frac_diff(val, r.val);

  s |=     ( fabs(val - r.val) > r.err );         /* consistent within error */
  s |= 2 * ( r.err <= 0.0 );                      /* error positive	     */
  s |= 4 * ( f > tol );                           /* within tolerance	     */
  s |= 8 * ( return_code != expect_return_code ); /* check return value      */

  gsl_test(s, message);

  if(s != 0) {
    printf("  expected: %20.16g\n", ref.val);
    printf("  obtained: %20.16g   %20.16g  %g\n", r.val, r.err, r.err/(fabs(r.val) + r.err));
    printf("  fracdiff: %20.16g\n", f);
  }
  if(s & 1) {
    printf("  value/expected not consistent within reported error\n");
  }
  if(s & 2) {
    printf("  reported error non-positive\n");
  }
  if(s & 4) {
    printf("  value not within tolerance of expected value\n");
  }
  if(s & 8) {
    printf("  nonzero return code: %d\n", return_code);
  }

  return (int) s;  
}


int check_bessel(void)
{
  gsl_sf_result r;
  double J[100];
  double Y[100];
  double I[100];
  double K[100];
  int s = 0;

  TEST_SF(s, gsl_sf_bessel_J0_e, (0.1, &r),     0.99750156206604003230,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_J0_e, (2.0, &r),     0.22389077914123566805,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_J0_e, (100.0, &r),   0.019985850304223122424,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_J0_e, (1.0e+10, &r), 2.1755917502468917269e-06, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_J1_e, (0.1, &r),      0.04993752603624199756,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_J1_e, (2.0, &r),      0.57672480775687338720,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_J1_e, (100.0, &r),   -0.07714535201411215803,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_J1_e, (1.0e+10, &r), -7.676508175684157103e-06, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Jn_e, (4, &r),     2.6028648545684032338e-07,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Jn_e, (5, &r),     0.007039629755871685484,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Jn_e, (100, &r),   0.09636667329586155967,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Jn_e, (1000, &r), -0.000011612675135233541081, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Y0_e, (0.1, &r),	    -1.5342386513503668441,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Y0_e, (2, &r),            0.5103756726497451196,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Y0_e, (256.0, &r),	    -0.03381290171792454909 ,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Y0_e, (4294967296.0, &r), 3.657903190017678681e-06, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Y1_e, (0.1, &r),         -6.45895109470202698800,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Y1_e, (2, &r),           -0.10703243154093754689,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Y1_e, (100.0, &r),       -0.020372312002759793305,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Y1_e, (4294967296.0, &r), 0.000011612249378370766284, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Yn_e, (4, 0.1, &r),	          -305832.29793353160319,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Yn_e, (5, 2, &r),              -9.935989128481974981,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Yn_e, (100, 100.0, &r),        -0.16692141141757650654,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Yn_e, (100, 4294967296.0, &r),  3.657889671577715808e-06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Yn_e, (1000, 4294967296.0, &r), 3.656551321485397501e-06, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_I0_scaled_e, (0.1, &r),     0.90710092578230109640,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_I0_scaled_e, (2, &r),       0.30850832255367103953,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_I0_scaled_e, (100.0, &r),   0.03994437929909668265,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_I0_scaled_e, (65536.0, &r), 0.0015583712551952223537, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_I1_scaled_e, (0.1, &r),     0.04529844680880932501,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_I1_scaled_e, (2, &r),       0.21526928924893765916,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_I1_scaled_e, (100.0, &r),   0.03974415302513025267,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_I1_scaled_e, (65536.0, &r), 0.0015583593657207350452, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_In_scaled_e, (  -4,    0.1), 2.3575258620054605307e-07, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_In_scaled_e, (   4,    0.1, &r), 2.3575258620054605307e-07, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_In_scaled_e, (   5,    2.0, &r), 0.0013297610941881578142, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_In_scaled_e, ( 100,  100.0, &r), 1.7266862628167695785e-22, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_I0_e, (0.1, &r), 1.0025015629340956014, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_I0_e, (2.0, &r), 2.2795853023360672674, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_I0_e, (100.0, &r), 1.0737517071310738235e, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_I1_e, (0.1, &r), 0.05006252604709269211, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_I1_e, (2.0, &r), 1.59063685463732906340, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_I1_e, (100.0, &r), 1.0683693903381624812e, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_In_e, (   4,    0.1, &r), 2.6054690212996573677e-07, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_In_e, (   5,    2.0, &r), 0.009825679323131702321, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_In_e, ( 100,  100.0, &r), 4.641534941616199114e, TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_bessel_K0_scaled_e, (0.1, &r), 2.6823261022628943831, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_K0_scaled_e, (2.0, &r), 0.8415682150707714179, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_K0_scaled_e, (100.0, &r), 0.1251756216591265789, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_K1_scaled_e, (0.1, &r), 10.890182683049696574, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_K1_scaled_e, (2.0, &r), 1.0334768470686885732, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_K1_scaled_e, (100.0, &r), 0.1257999504795785293, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Kn_scaled_e, (   4,    0.1, &r), 530040.2483725626207, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Kn_scaled_e, (   5,    2.0, &r), 69.68655087607675118, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Kn_scaled_e, ( 100,  100.0, &r), 2.0475736731166756813e+19, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_K0_e, (0.1, &r), 2.4270690247020166125, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_K0_e, (2.0, &r), 0.11389387274953343565, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_K0_e, (100.0, &r), 4.656628229175902019e-45, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_K1_e, (0.1, &r), 9.853844780870606135, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_K1_e, (2.0, &r), 0.13986588181652242728, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_K1_e, (100.0, &r), 4.679853735636909287e-45, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Kn_e, (   4,    0.1, &r), 479600.2497925682849, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Kn_e, (   5,    2.0, &r), 9.431049100596467443, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Kn_e, ( 100,  100.0, &r), 7.617129630494085416e-25, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_bessel_j0_e, (-10.0, &r), -0.05440211108893698134, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j0_e, (0.001, &r), 0.9999998333333416667, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j0_e, (  1.0, &r), 0.84147098480789650670, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j0_e, ( 10.0, &r), -0.05440211108893698134, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j0_e, (100.0, &r), -0.005063656411097587937, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j0_e, (1048576.0, &r), 3.1518281938718287624e-07, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_bessel_j1_e, (-10.0, &r), -0.07846694179875154709, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j1_e, (0.01, &r), 0.003333300000119047399, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j1_e, (  1.0, &r), 0.30116867893975678925, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j1_e, ( 10.0, &r), 0.07846694179875154709, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j1_e, (100.0, &r), -0.008673825286987815220, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j1_e, (1048576.0, &r), -9.000855242905546158e-07, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_bessel_j2_e, (-10.0, &r), 0.07794219362856244547, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j2_e, (0.01, &r), 6.666619047751322551e-06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j2_e, (  1.0, &r), 0.06203505201137386110, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j2_e, ( 10.0, &r), 0.07794219362856244547, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j2_e, (100.0, &r), 0.004803441652487953480, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_j2_e, (1048576.0, &r), -3.1518539455252413111e-07, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_bessel_jl_e, (1,       10.0, &r), 0.07846694179875154709000, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_jl_e, (5,        1.0, &r), 0.00009256115861125816357, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_jl_e, (10,      10.0, &r), 0.06460515449256426427, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_jl_e, (100,    100.0, &r), 0.010880477011438336539, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_jl_e, (2000, 1048576.0, &r), 7.449384239168568534e-07, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_bessel_y0_e, (0.001, &r), -999.99950000004166670, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y0_e, (  1.0, &r), -0.5403023058681397174, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y0_e, ( 10.0, &r), 0.08390715290764524523, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y0_e, (100.0, &r), -0.008623188722876839341, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y0_e, (65536.0, &r), 0.000011014324202158573930, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y0_e, (4294967296.0, &r), 2.0649445131370357007e-10, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_bessel_y1_e, ( 0.01, &r), -10000.499987500069444, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y1_e, (  1.0, &r), -1.3817732906760362241, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y1_e, ( 10.0, &r), 0.06279282637970150586, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y1_e, (100.0, &r), 0.004977424523868819543, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y1_e, (4294967296.0, &r), 1.0756463271573404688e-10, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_bessel_y2_e, ( 0.01, &r), -3.0000500012499791668e+06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y2_e, (  1.0, &r), -3.605017566159968955, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y2_e, ( 10.0, &r), -0.06506930499373479347, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y2_e, (100.0, &r), 0.008772511458592903927, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_bessel_y2_e, (4294967296.0, &r), -2.0649445123857054207e-10, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_i0_scaled_e, (0.1, &r), 0.9063462346100907067, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_i0_scaled_e, (2.0, &r), 0.24542109027781645493, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_i0_scaled_e, (100.0, &r), 0.005000000000000000000, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_i1_scaled_e, (0.1, &r), 0.030191419289002226846, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_i1_scaled_e, (2.0, &r), 0.131868364583275317610, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_i1_scaled_e, (100.0, &r), 0.004950000000000000000, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_i2_scaled_e, (0.1, &r), 0.0006036559400239012567, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_i2_scaled_e, (2.0, &r), 0.0476185434029034785100, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_i2_scaled_e, (100.0, &r), 0.0048515000000000000000, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_il_scaled_e, (   4, 0.001, &r), 1.0571434341190365013e-15, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_il_scaled_e, (   4,   0.1, &r), 9.579352242057134927e-08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_il_scaled_e, (   5,   2.0, &r), 0.0004851564602127540059, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_il_scaled_e, (   5, 100.0, &r), 0.004300446777500000000, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_il_scaled_e, ( 100, 100.0, &r), 1.3898161964299132789e-23, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_k0_scaled_e, (0.1, &r), 15.707963267948966192, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_k0_scaled_e, (2.0, &r), 0.7853981633974483096, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_k0_scaled_e, (100.0, &r), 0.015707963267948966192, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_k1_scaled_e, (0.1, &r), 172.78759594743862812, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_k1_scaled_e, (2.0, &r), 1.1780972450961724644, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_k1_scaled_e, (100.0, &r), 0.015865042900628455854, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_k2_scaled_e, (0.1, &r), 5199.335841691107810, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_k2_scaled_e, (2.0, &r), 2.5525440310417070063, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_k2_scaled_e, (100.0, &r), 0.016183914554967819868, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_kl_scaled_e, (   4, 1.0/256.0, &r), 1.8205599816961954439e+14, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_kl_scaled_e, (   4, 1.0/8.0, &r), 6.1173217814406597530e+06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_kl_scaled_e, (   5,   2.0, &r), 138.10735829492005119, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_kl_scaled_e, ( 100, 100.0, &r), 3.985930768060258219e+18, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Jnu_e, (0.0001,10.0, &r), -0.2459270166445205, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Jnu_e, ( 1.0, 0.001, &r), 0.0004999999375000026, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Jnu_e, ( 1.0,   1.0, &r), 0.4400505857449335160, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Jnu_e, (30.0,   1.0, &r), 3.482869794251482902e-42, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Jnu_e, (30.0, 100.0, &r), 0.08146012958117222297, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Jnu_e, (10.0,   1.0, &r), 2.6306151236874532070e-10, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Jnu_e, (10.0, 100.0, &r), -0.05473217693547201474, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Jnu_e, (10.2, 100.0, &r), -0.03548919161046526864, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Ynu_e, (0.0001,10.0, &r), 0.05570979797521875261, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Ynu_e, ( 1.0, 0.001, &r), -636.6221672311394281, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Ynu_e, ( 1.0,   1.0, &r), -0.7812128213002887165, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Ynu_e, (30.0,   1.0, &r), -3.0481287832256432162e+39, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Ynu_e, (30.0, 100.0, &r), 0.006138839212010033452, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Ynu_e, (10.0,   1.0, &r), -1.2161801427868918929e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Ynu_e, (10.0, 100.0, &r), 0.05833157423641492875, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Ynu_e, (10.2, 100.0, &r), 0.07169383985546287091, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Inu_scaled_e, (0.0001,10.0, &r), 0.12783333709581669672, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_scaled_e, ( 1.0, 0.001, &r), 0.0004995003123542213370, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_scaled_e, ( 1.0,   1.0, &r), 0.20791041534970844887, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_scaled_e, (30.0,   1.0, &r), 1.3021094983785914437e-42, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_scaled_e, (30.0, 100.0, &r), 0.0004486987756920986146, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_scaled_e, (10.0,   1.0, &r), 1.0127529864692066036e-10, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_scaled_e, (10.0, 100.0, &r), 0.024176682718258828365, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_scaled_e, (10.2, 100.0, &r), 0.023691628843913810043, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Inu_e, (0.0001,10.0, &r), 2815.7166269770030352, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_e, ( 1.0, 0.001, &r), 0.0005000000625000026042, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_e, ( 1.0,   1.0, &r), 0.5651591039924850272, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_e, (30.0,   1.0, &r), 3.539500588106447747e-42, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_e, (30.0, 100.0, &r), 1.2061548704498434006e+40, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_e, (10.0,   1.0, &r), 2.7529480398368736252e-10, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_e, (10.0, 100.0, &r), 6.498975524720147799e+41, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Inu_e, (10.2, 100.0, &r), 6.368587361287030443e+41, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Knu_scaled_e, (0.0001,10.0, &r), 0.3916319346235421817, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_scaled_e, ( 1.0, 0.001, &r), 1000.9967345590684524, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_scaled_e, ( 1.0,   1.0, &r), 1.6361534862632582465, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_scaled_e, (30.0,   1.0, &r), 1.2792629867539753925e+40, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_scaled_e, (30.0, 100.0, &r), 10.673443449954850040, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_scaled_e, (10.0,   1.0, &r), 4.912296520990198599e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_scaled_e, (10.0, 100.0, &r), 0.20578687173955779807, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_scaled_e, (10.0, 1000.0, &r), 0.04165905142800565788, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_scaled_e, (10.0, 1.0e+8, &r), 0.00012533147624060789938, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_scaled_e, (10.2, 100.0, &r), 0.20995808355244385075, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_Knu_e, (0.0001,0.001, &r), 7.023689431812884141, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_e, (0.0001,10.0, &r), 0.000017780062324654874306, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_e, ( 1.0, 0.001, &r), 999.9962381560855743, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_e, ( 1.0,   1.0, &r), 0.6019072301972345747, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_e, (10.0, 0.001, &r), 1.8579455483904008064e+38, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_e, (10.0,   1.0, &r), 1.8071328990102945469e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_e, (10.0, 100.0, &r), 7.655427977388100611e-45, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_e, (10.2, 100.0, &r), 7.810600225948217841e-45, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_e, (30.0,   1.0, &r), 4.706145526783626883e+39, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_Knu_e, (30.0, 100.0, &r), 3.970602055959398739e-43, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (0.0001,1.0e-100, &r), 5.439794449319847, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (0.0001,0.0001, &r), 2.232835507214331, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (0.0001,10.0, &r), -10.93743282256098, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, ( 1.0, 1.0e-100, &r), 230.2585092994045, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, ( 1.0, 1.0e-10, &r), 23.025850929940456840, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, ( 1.0, 0.001, &r), 6.907751517131146, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, ( 1.0,   1.0, &r), -0.5076519482107523309, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (30.0, 1.0e-100, &r), 6999.113586185543475, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (30.0,   1.0, &r), 91.34968784026325464, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (30.0, 100.0, &r), -97.63224126416760932, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (100.0, 1.0e-100, &r), 23453.606706185466825, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (100.0, 1.0, &r), 427.7532510250188083, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (100.0, 100.0, &r), -55.53422771502921431, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (1000.0, 1.0e-100, &r), 236856.183755993135, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_bessel_lnKnu_e, (10000.0, 1.0e-100, &r), 2.39161558914890695e+06, TEST_TOL0, GSL_SUCCESS);

  gsl_sf_bessel_Jn_array_impl(3, 38, 1.0, J);
  s += ( frac_diff(J[0],  0.0195633539826684059190  ) >  TEST_TOL0);
  s += ( frac_diff(J[1],  0.0024766389641099550438  ) >  TEST_TOL0);
  s += ( frac_diff(J[10], 1.9256167644801728904e-14 ) >  TEST_TOL0);
  s += ( frac_diff(J[35], 6.911232970971626272e-57  ) >  TEST_TOL0);
  gsl_test(s, "  gsl_sf_bessel_Jn_array_impl");
  status += s;

  gsl_sf_bessel_Yn_array_impl(3, 38, 1.0, Y);
  s += ( frac_diff(Y[0],  -5.821517605964728848      ) > TOL );
  s += ( frac_diff(Y[1],  -33.27842302897211870      ) > TOL );
  s += ( frac_diff(Y[10], -1.2753618701519837595e+12 ) > TOL );
  s += ( frac_diff(Y[35], -1.2124435663593357154e+54 ) > TOL );
  gsl_test(s, "  gsl_sf_bessel_Yn_array_impl");
  status += s;

  gsl_sf_bessel_In_scaled_array_impl(3, 38, 1.0, I);
  s += ( frac_diff(I[0],  0.0081553077728142938170  ) > TOL );
  s += ( frac_diff(I[1],  0.0010069302573377758637  ) > TOL );
  s += ( frac_diff(I[10], 7.341518665628926244e-15  ) > TOL );
  s += ( frac_diff(I[35], 2.5753065298357542893e-57 ) > TOL );
  gsl_test(s, "  gsl_sf_bessel_In_scaled_array_impl");
  status += s;

  gsl_sf_bessel_In_array_impl(3, 38, 1.0, Y);
  s += ( frac_diff(Y[0],   0.0221684249243319024760  ) > TOL );
  s += ( frac_diff(Y[1],   0.0027371202210468663251  ) > TOL );
  s += ( frac_diff(Y[10],  1.9956316782072007564e-14 ) > TOL );
  s += ( frac_diff(Y[35],  7.000408942764452901e-57  ) > TOL );
  gsl_test(s, "  gsl_sf_bessel_In_array_impl");
  status += s;

  gsl_sf_bessel_Kn_array_impl(3, 38, 1.0, K);
  s += ( frac_diff(K[0],  7.101262824737944506 ) > TOL );
  s += ( frac_diff(K[1],  44.23241584706284452 ) > TOL );
  s += ( frac_diff(K[10], 1.9215763927929940846e+12 ) > TOL );
  s += ( frac_diff(K[35], 1.8789385023806051223e+54 ) > TOL );
  gsl_test(s, "  gsl_sf_bessel_Kn_array_impl");
  status += s;

  gsl_sf_bessel_Kn_scaled_array_impl(3, 38, 1.0, K);
  s += ( frac_diff(K[0],  19.303233695596904277 ) > TOL );
  s += ( frac_diff(K[1],  120.23617222591483717 ) > TOL );
  s += ( frac_diff(K[10], 5.223386190525076473e+12 ) > TOL );
  s += ( frac_diff(K[35], 5.107484387813251411e+54 ) > TOL );
  gsl_test(s, "  gsl_sf_bessel_Kn_scaled_array_impl");
  status += s;

  gsl_sf_bessel_jl_array_impl(50, 1.0, J);
  s += ( frac_diff(J[0],  0.84147098480789650670   ) > TOL );
  s += ( frac_diff(J[1],  0.30116867893975678925   ) > TOL );
  s += ( frac_diff(J[10], 7.116552640047313024e-11 ) > TOL );
  s += ( frac_diff(J[50], 3.615274717489787311e-81 ) > TOL );
  gsl_test(s, "  gsl_sf_bessel_jl_array_impl");
  status += s;

  gsl_sf_bessel_jl_steed_array_impl(50, 1.0, J);
  s += ( frac_diff(J[0],  0.84147098480789650670   ) > TOL );
  s += ( frac_diff(J[1],  0.30116867893975678925   ) > TOL );
  s += ( frac_diff(J[10], 7.116552640047313024e-11 ) > TOL );
  s += ( frac_diff(J[50], 3.615274717489787311e-81 ) > TOL );
  gsl_test(s, "  gsl_sf_bessel_jl_steed_array_impl");
  status += s;

  gsl_sf_bessel_yl_array_impl(50, 1.0, Y);
  s += ( frac_diff(Y[0],  -0.5403023058681397174 ) > TOL );
  s += ( frac_diff(Y[1],  -1.3817732906760362241 ) > TOL );
  s += ( frac_diff(Y[10], -6.722150082562084436e+08  ) > TOL );
  s += ( frac_diff(Y[50], -2.7391922846297571576e+78 ) > TOL );
  gsl_test(s, "  gsl_sf_bessel_yl_array_impl");
  status += s;

  gsl_sf_bessel_il_scaled_array_impl(50, 1.0, I);
  s += ( frac_diff(I[0],  0.43233235838169365410 ) > TOL );
  s += ( frac_diff(I[1],  0.13533528323661269189 ) > TOL );
  s += ( frac_diff(I[10], 2.7343719371837065460e-11 ) > TOL );
  s += ( frac_diff(I[50], 1.3429606061892023653e-81 ) > TOL );
  gsl_test(s, "  gsl_sf_bessel_il_scaled_array_impl");
  status += s;
  
  gsl_sf_bessel_kl_scaled_array_impl(50, 1.0, K);
  s += ( frac_diff(K[0],  1.5707963267948966192     ) > TOL );
  s += ( frac_diff(K[1],  3.1415926535897932385     ) > TOL );
  s += ( frac_diff(K[10], 2.7231075458948147010e+09 ) > TOL );
  s += ( frac_diff(K[50], 1.1578440432804522544e+79 ) > TOL );
  gsl_test(s, "  gsl_sf_bessel_kl_scaled_array_impl");
  status += s;

  return status;
}

int check_cheb(void)
{
  gsl_sf_result r;
  double x;
  double f;
  int status = 0;
  int s;

  gsl_sf_cheb_series * cs = gsl_sf_cheb_new(sin, -M_PI, M_PI, 40);

  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    gsl_sf_cheb_eval_impl(cs, x, &r);
    f += fabs(r.val - sin(x));
  }
  s = 0;
  s += ( f > 100.0 * TOL );
  gsl_test(s, "  gsl_sf_cheb_eval_impl()");
  status += s;
  
  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    gsl_sf_cheb_eval_n_impl(cs, 25, x, &r);
    f += fabs(r.val - sin(x));
  }
  s = 0;
  s += ( f > 100.0 * TOL );
  gsl_test(s, "  gsl_sf_cheb_eval_n_impl()");
  status += s;

  gsl_sf_cheb_calc_impl(cs, sin);
  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    gsl_sf_cheb_eval(cs, x, &r);
    f += fabs(r.val - sin(x));
  }
  s = 0;
  s += ( f > 100.0 * TOL );
  gsl_test(s, "  gsl_sf_cheb_calc_impl()");
  status += s;

  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    gsl_sf_cheb_eval_deriv_impl(cs, x, &r);
    f += fabs(r.val - cos(x));
  }
  s = 0;
  s += ( f > 100.0 * 10.0 * TOL );
  gsl_test(s, "  gsl_sf_cheb_eval_deriv_impl()");
  status += s;

  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    gsl_sf_cheb_eval_integ(cs, x, &r);
    f += fabs(r.val + (1.0 + cos(x)));
  }
  s = 0;
  s += ( f > 100.0 * TOL );
  gsl_test(s, "  gsl_sf_cheb_eval_integ()");
  status += s;

  gsl_sf_cheb_free(cs);

  return status;
}


int check_clausen(void)
{
  gsl_sf_result r;
  double y;
  int status = 0;
  int s;

  TEST_SF(s,  gsl_sf_clausen_e, (M_PI/20.0, &r), 0.4478882448133546, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_clausen_e, (M_PI/6.0, &r), 0.8643791310538927, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_clausen_e, (M_PI/3.0, &r), 1.0149416064096535, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_clausen_e, (  2.0*M_PI + M_PI/3.0, &r), 1.0149416064096535, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_clausen_e, (100.0*M_PI + M_PI/3.0, &r), 1.0149416064096535, TEST_TOL0, GSL_SUCCESS);

  return status;
}


#define PRINT(n) printf("%22.18g  %22.18g  %22.18g  %22.18g\n", F[n], Fp[n], G[n], Gp[n])

int check_coulomb(void)
{
  int status = 0;
  int s;
  
  const int kmax = 20;
  double F[kmax+1], Fp[kmax+1], G[kmax+1], Gp[kmax+1];
  double Fe, Ge;
  double lam_min;
  double lam_F;
  double eta, x;
  int k_G;

  TEST_SF(s, gsl_sf_hydrogenicR_1_e, (3.0, 2.0, &r), 0.025759948256148471036, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_1_e, (3.0, 10.0, &r), 9.724727052062819704e-13, TEST_TOL0, GSL_SUCCESS);  status += s;

  TEST_SF(s, gsl_sf_hydrogenicR_e, (4, 0, 3.0, 2.0, &r), -0.03623182256981820062, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_e, (4, 1, 3.0, 2.0, &r), -0.028065049083129581005, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_e, (4, 2, 3.0, 2.0, &r), 0.14583027278668431009, TEST_TOL0, GSL_SUCCESS);  status += s;

  TEST_SF(s, gsl_sf_hydrogenicR_e, (100,  0, 3.0, 2.0, &r), -0.00007938950980052281367, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_e, (100, 10, 3.0, 2.0, &r), 7.112823375353605977e-12, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_e, (100, 90, 3.0, 2.0, &r), 5.845231751418131548e-245, TEST_TOL0, GSL_SUCCESS);

#if 0
  lam_F = 0.0;
  k_G   = 0;
  eta = 1.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  0.6849374120059439677 ) > 1.e-10 );
  s += ( frac_diff( Fp[0], -0.7236423862556063963 ) > 1.e-10 );
  s += ( frac_diff(  G[0], -0.8984143590920205487 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -0.5108047585190350106 ) > 1.e-10 );
  gsl_test(s,"  gsl_sf_coulomb_wave_FG_impl(1.0, 5.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 10.0;
  k_G   = 2;
  eta = 1.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  0.0006423773354915823698 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  0.0013299570958719702545 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  33.27615734455096130     ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -45.49180102261540580     ) > 1.e-10 );
  gsl_test(s,"  gsl_sf_coulomb_wave_FG_impl(1.0, 5.0, lam_F=10, lam_G=8)");
  status += s;

  lam_F = 4.0;
  k_G   = 2;
  eta = 50.0;
  x = 120.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0], 0.0735194711823798495 ) > 1.e-10 );
  s += ( frac_diff( Fp[0], 0.6368149124126783325 ) > 1.e-10 );
  /*
  s += ( frac_diff(  G[0],  ) > 1.e-10 );
  s += ( frac_diff( Gp[0],  ) > 1.e-10 );
  */
  gsl_test(s,"  gsl_sf_coulomb_wave_FG_impl(50.0, 120.0, lam_F=4, lam_G=2)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = -1000.0;
  x = 1.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  9.68222518991e-02 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  5.12063396274e+00 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.13936784380e-01 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -4.30243486522e+00 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(-1000.0, 1.0, lam_F=0, lam_G=0)");
  status += s;

  lam_min = 0.0;
  eta = -50.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  1.52236975714e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  2.03091041166e+00 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  4.41680690236e-01 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -6.76485374767e-01 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(-50.0, 5.0, lam_F=0, lam_G=0)");
  status += s;

  lam_min = 0.0;
  eta = -50.0;
  x = 1000.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0], -0.2267212182760888523 ) > 1.e-10 );
  s += ( frac_diff( Fp[0], -0.9961306810018401525 ) > 1.e-10 );
  s += ( frac_diff(  G[0], -0.9497684438900352186 ) > 1.e-10 );
  s += ( frac_diff( Gp[0],  0.2377656295411961399 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(-50.0, 1000.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 10.0;
  k_G = 0;
  eta = -50.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0], -3.68114360218e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  1.33846751032e+00 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  3.31588324611e-01 ) > 1.e-10 );
  s += ( frac_diff( Gp[0],  1.51088862814e+00 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(-50.0, 5.0, lam_F=10, lam_G=10)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = -4.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  4.07862723006e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  1.09821233636e+00 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  6.74327035383e-01 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -6.36110427280e-01 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(-4.0, 5.0, lam_F=0, lam_G=0");
  status += s;

  lam_F = 3.0;
  k_G = 0;
  eta = -4.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0], -2.56863093558e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  1.14322942201e+00 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  7.87989922393e-01 ) > 1.e-10 );
  s += ( frac_diff( Gp[0],  3.85990587811e-01 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(-4.0, 5.0, lam_F=3, lam_G=3");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 1.0;
  x = 2.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  6.61781613833e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  4.81557455710e-01 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.27577878477e+00 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -5.82728813097e-01 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(1.0, 2.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 1.0;
  x = 0.5;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  0.08315404535022023302  ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  0.22693874616222787568  ) > 1.e-10 );
  s += ( frac_diff(  G[0],  3.1060069279548875140   ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -3.549156038719924236    ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(1.0, 0.5, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.5;
  k_G = 0;
  eta = 1.0;
  x = 0.5;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  0.04049078073829290935 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  0.14194939168094778795 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  4.720553853049677897   ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -8.148033852319180005   ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(1.0, 0.5, lam_F=0.5, lam_G=0.5)");
  status += s;

  lam_F = 0.1;
  k_G = 0;
  eta = 1.0;
  x = 0.5;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  0.07365466672379703418 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  0.21147121807178518647 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  3.306705446241024890   ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -4.082931670935696644   ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(1.0, 0.5, lam_F=0.1, lam_G=0.1)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 8.0;
  x = 1.05;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  9.882706082810274357e-09 ) > 1.e-10);
  s += ( frac_diff( Fp[0],  4.005167028235547770e-08 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.333127992006686320e+07 ) > 1.e-6 );
  s += ( frac_diff( Gp[0], -4.715914530842402330e+07 ) > 1.e-6 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(8.0, 1.05, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.1;
  k_G = 0;
  eta = 8.0;
  x = 1.05;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  9.611416736061987761e-09  ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  3.909628126126824140e-08  ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.3659284642192625811e+07 ) > 1.e-6 );
  s += ( frac_diff( Gp[0], -4.848117385783386850e+07  ) > 1.e-6 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(8.0, 1.05, lam_F=0.1, lam_G=0.1)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 50.0;
  x = 0.1;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  2.807788027954216071e-67 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  9.677600748751576606e-66 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  5.579810686998358766e+64 ) > 1.e-8 );
  s += ( frac_diff( Gp[0], -1.638329512756321424e+66 ) > 1.e-8 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(50.0, 0.1, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 10.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  1.7207454091787930614e-06 ) > 1.e-03 );
  s += ( frac_diff( Fp[0],  3.0975994706405458046e-06 ) > 1.e-03 );
  s += ( frac_diff(  G[0],  167637.56609459967623     ) > 1.e-03 );
  s += ( frac_diff( Gp[0], -279370.76655361803075     ) > 1.e-03 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(10.0, 5.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 25.0;
  x = 10.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  1.5451274501076114315e-16 ) > 1.e-03 );
  s += ( frac_diff( Fp[0],  3.1390869393378630928e-16 ) > 1.e-03 );
  s += ( frac_diff(  G[0],  1.6177129008336318136e+15 ) > 1.e-03 );
  s += ( frac_diff( Gp[0], -3.1854062013149740860e+15 ) > 1.e-03 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(25.0, 10.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 1.0;
  x = 9.2;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0], -0.25632012319757955655 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  0.91518792286724220370 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.03120585918973466110 ) > 1.e-10 );
  s += ( frac_diff( Gp[0],  0.21946326717491250193 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(1.0, 9.2, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  eta = 10.0;
  x = 10.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  0.0016262711250135878249 ) > 1.e-03 );
  s += ( frac_diff( Fp[0],  0.0017060476320792806014 ) > 1.e-03 );
  s += ( frac_diff(  G[0],  307.87321661090837987    ) > 1.e-03 );
  s += ( frac_diff( Gp[0], -291.92772380826822871    ) > 1.e-03 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(10.0, 10.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  eta = 100.0;
  x = 1.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  8.999367996930662705e-126 ) > 1.e-03 );
  s += ( frac_diff( Fp[0],  1.292746745757069321e-124 ) > 1.e-03 );
  s += ( frac_diff(  G[0],  3.936654148133683610e+123 ) > 1.e-03 );
  s += ( frac_diff( Gp[0], -5.456942268061526371e+124 ) > 1.e-03 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(100.0, 1.0, lam_F=0, lam_G=0)");
  status += s;
#endif

  return status;
}

int check_coupling(void)
{
  gsl_sf_result y;
  int s = 0;

  TEST_SF(s, gsl_sf_coupling_3j_impl, (0, 1, 1, 0, 1, -1, &y), sqrt(1.0/2.0), TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_coupling_3j_impl, (1, 1, 2, 1, -1, 0, &y), sqrt(1.0/6.0), TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_coupling_3j_impl, (2, 4, 6, 0, 2, -2, &y), sqrt(8.0/105.0), TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_coupling_6j_impl, (2, 2, 4, 2, 2, 2, &y), 1.0/6.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_coupling_9j_impl, (4, 2, 4, 3, 3, 2, 1, 1, 2, &y), -0.040824829046386, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_coupling_9j_impl, (8, 4, 10, 7, 3, 8, 1, 1, 2, &y), 0.025458753860866, TEST_TOL0, GSL_SUCCESS);

  return status;
}

int check_dawson(void)
{
  gsl_sf_result r;
  int status = 0;
  int s;

  TEST_SF(s,  gsl_sf_dawson_e, (1.0e-15, &r), 1.0e-15, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dawson_e, (0.5, &r), 0.4244363835020222959, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dawson_e, (2.0, &r), 0.30134038892379196603, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dawson_e, (1000.0, &r), 0.0005000002500003750009, TEST_TOL0, GSL_SUCCESS);
  
  return status;
}

int check_debye(void)
{
  gsl_sf_result y;
  int s = 0;

  /* FIXME: I do not have more accurate test values than these.
   */

  TEST_SF(s, gsl_sf_debye_1_impl, (0.1, &y), 0.975278, 1.0e-5, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_1_impl, (1.0, &y), 0.777505, 1.0e-5, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_1_impl, (10.0, &y), 0.164443, 1.0e-5, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_debye_2_impl, (0.1, &y), 0.967083, 1.0e-5, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_2_impl, (1.0, &y), 0.707878, 1.0e-5, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_2_impl, (10.0, &y), 0.047971, 1.0e-5, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_debye_3_impl, (0.1, &y), 0.963000, 1.0e-5, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_3_impl, (1.0, &y), 0.674416, 1.0e-5, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_3_impl, (10.0, &y), 0.019296, 1.0e-5, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_debye_4_impl, (0.1, &y), 0.960555, 1.0e-5, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_4_impl, (1.0, &y), 0.654874, 1.0e-5, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_4_impl, (10.0, &y), 0.009674, 1.0e-5, GSL_SUCCESS);

  return status;
}


int check_dilog(void)
{
  double x, y;
  int status = 0;
  int s;

  /* real dilog */

  TEST_SF(s,  gsl_sf_dilog_e, (-3.0, &r), -1.9393754207667089531, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (-0.5, &r), -0.4484142069236462024, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (-0.001, &r), -0.0009997501110486510834, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (0.1, &r), 0.1026177910993911, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (0.7, &r), 0.8893776242860387386, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (1.0, &r), 1.6449340668482260, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (1.5, &r), 2.3743952702724802007, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (2.0, &r), 2.4674011002723397, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, ( 5.0, &r), 1.7837191612666306277, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, ( 11.0, &r), 0.3218540439999117111, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (12.59, &r), 0.0010060918167266208634, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (12.595, &r), 0.00003314826006436236810, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (13.0, &r), -0.07806971248458575855, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (20.0, &r), -1.2479770861745251168, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (150.0, &r), -9.270042702348657270, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_e, (1100.0, &r), -21.232504073931749553, TEST_TOL0, GSL_SUCCESS);


  /* complex dilog */
  /* FIXME: probably need more tests here... 
   * also need to work on accuracy for r->1; need to
   * adjust the switch-over point I suppose.
   */

  s = 0;
  gsl_sf_complex_dilog_impl(1.00001, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.20562022409960237363 ) > TOL );
  s += ( frac_diff( y,  0.91597344814458309320 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(1.00001, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.99999, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.20561329262779687646 ) > TOL );
  s += ( frac_diff( y,  0.91595774018131512060 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.99999, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.991, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.20250384721077806127 ) > 1.0e-06 );
  s += ( frac_diff( y,  0.90888544355846447810 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.991, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.98, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.19871638377785918403 ) > 1.0e-05 );
  s += ( frac_diff( y,  0.90020045882981847610 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.98, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.95, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.18848636456893572091 ) > TOL );
  s += ( frac_diff( y,  0.87633754133420277830 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.95, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.8, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.13980800855429037810 ) > TOL );
  s += ( frac_diff( y,  0.75310609092419884460 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.8, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.5, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.05897507442156586346 ) > TOL );
  s += ( frac_diff( y,  0.48722235829452235710 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.5, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.01, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.000024999375027776215378 ) > 1.0e-12 );
  s += ( frac_diff( y,  0.009999888892888684820    ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.01, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(10.0, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -3.0596887943287347304 ) > TOL );
  s += ( frac_diff( y,  3.7167814930680685900 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(10.0, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(100.0, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -11.015004738293824854 ) > TOL );
  s += ( frac_diff( y,  7.2437843013083534970 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(100.0, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.9, M_PI/8.0, &x, &y);
  s += ( frac_diff( x, 0.9659561692810778695 ) > TOL );
  s += ( frac_diff( y, 0.6295758865676813951 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.9, Pi/8)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.99, M_PI/8.0, &x, &y);
  s += ( frac_diff( x, 1.0571539648820244720 ) > TOL );
  s += ( frac_diff( y, 0.7469145254610851318 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.99, Pi/8)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.99, M_PI/64.0, &x, &y);
  s += ( frac_diff( x, 1.5381800285902999666 ) > TOL );
  s += ( frac_diff( y, 0.1825271634987756651 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.99, Pi/64)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.9, 3.0*M_PI/4.0, &x, &y);
  s += ( frac_diff( x, -0.6062840301356530985 ) > TOL );
  s += ( frac_diff( y,  0.4836632833122775721 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_dilog(0.9, 3Pi/4)");
  status += s;

  return status;
}


int check_elementary(void)
{
  double y;
  double x = 0.2*DBL_MAX;
  int s = 0;

  TEST_SF(s,  gsl_sf_multiply_e, (-3.0,2.0, &r), -6.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_multiply_e, (x, 1.0/x, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_multiply_e, (x, 0.2, &r), 0.2, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_multiply_e, (x, 2.0, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  s += ( gsl_sf_multiply_impl(DBL_MAX, 1.1, &r) != GSL_EOVRFLW);
  s += ( gsl_sf_multiply_impl(DBL_MIN, 0.9, &r) != GSL_EUNDRFLW);

  return s;
}


int check_ellint(void)
{
  gsl_prec_t goal = GSL_DOUBLE_PREC;
  int status = 0;
  int s = 0;
  
  TEST_SF(s,  gsl_sf_ellint_Kcomp_e, ( 0.99, goal, &r), 3.3566005233611923760, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_Kcomp_e, ( 0.50, goal, &r), 1.6857503548125960429, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_Kcomp_e, (0.010, goal, &r), 1.5708355989121522360, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_Ecomp_e, (0.99, goal, &r), 1.0284758090288040010, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_Ecomp_e, (0.50, goal, &r), 1.4674622093394271555, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_Ecomp_e, (0.01, goal, &r), 1.5707570561503852873, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_F_e, (M_PI/3.0, 0.99, goal, &r), 1.3065333392738766762, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_F_e, (M_PI/3.0, 0.50, goal, &r), 1.0895506700518854093, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_F_e, (M_PI/3.0, 0.01, goal, &r), 1.0472129063770918952, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_E_e, (M_PI/3.0, 0.99, goal, &r), 0.8704819220377943536, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_E_e, (M_PI/3.0, 0.50, goal, &r), 1.0075555551444720293, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_E_e, (M_PI/3.0, 0.01, goal, &r), 1.0471821963889481104, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_P_e, (M_PI/3.0, 0.99, 0.5, goal, &r), 1.1288726598764099882, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_P_e, (M_PI/3.0, 0.50, 0.5, goal, &r), 0.9570574331323584890, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_P_e, (M_PI/3.0, 0.01, 0.5, goal, &r), 0.9228868127118118465, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RF_e, (5.0e-11, 1.0e-10, 1.0, goal, &r), 12.36441982979439, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RF_e, (1.0, 2.0, 3.0, goal, &r), 0.7269459354689082, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RD_e, (5.0e-11, 1.0e-10, 1.0, goal, &r), 34.0932594919337362, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RD_e, (1.0, 2.0, 3.0, goal, &r), 0.2904602810289906, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RC_e, (1.0, 2.0, goal, &r), 0.7853981633974482, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RJ_e, (2.0, 3.0, 4.0, 5.0, goal, &r), 0.1429757966715675, TEST_TOL0, GSL_SUCCESS);


  return s;
}


int check_erf(void)
{
  gsl_sf_result r;
  int s = 0;

  TEST_SF(s, gsl_sf_erfc_e, (-10.0, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_e, (-1.0, &r), 1.8427007929497148693, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_e, (-0.5, &r), 1.5204998778130465377, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_e, (1.0, &r), 0.15729920705028513066, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_e, (3.0, &r), 0.000022090496998585441373, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_e, (7.0, &r), 4.183825607779414399e-23, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_e, (10.0, &r), 2.0884875837625447570e-45, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log_erfc_e, (-1.0, &r), log(1.842700792949714869) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_erfc_e, (1.0, &r), log(0.15729920705028513066) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_erfc_e, (10.0, &r), log(2.0884875837625447570e-45) , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_erf_e, (-10.0, &r), -1.0000000000000000000, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erf_e, (0.5, &r), 0.5204998778130465377, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erf_e, (1.0, &r), 0.8427007929497148693, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erf_e, (10.0, &r), 1.0000000000000000000, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_erf_Z_e, (1.0, &r), 0.24197072451914334980, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erf_Q_e, (10.0, &r), 0.5 * gsl_sf_erfc(10.0/M_SQRT2) , TEST_TOL0, GSL_SUCCESS);

  return s;
}


int check_exp(void)
{
  gsl_sf_result r;
  double x;
  int status = 0;
  int s = 0;

  TEST_SF(s,  gsl_sf_exp_e, (-10.0, &r), exp(-10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exp_e, ( 10.0, &r), exp( 10.0) , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exp_sgn_e, (-10.0, 2.0, &r), exp(-10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exp_sgn_e, ( 10.0, 3.7, &r), exp( 10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exp_sgn_e, (-10.0, -2.0, &r), -exp(-10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exp_sgn_e, ( 10.0, -3.7, &r), -exp( 10.0) , TEST_TOL0, GSL_SUCCESS);

  x = 0.8*GSL_LOG_DBL_MAX;
  TEST_SF(s, gsl_sf_exp_mult_e, (-10.0,  1.0e-06, &r), 1.0e-06*exp(-10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (-10.0,  2.0, &r), 2.0*exp(-10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (-10.0, -2.0, &r), -2.0*exp(-10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, ( 10.0,  1.0e-06, &r), 1.0e-06*exp( 10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, ( 10.0, -2.0, &r), -2.0*exp( 10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, 1.00001, &r), 1.00001*exp(x) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, 1.000001, &r), 1.000001*exp(x) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, 1.000000001, &r), 1.000000001*exp(x) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, 100.0, &r), 100.0*exp(x) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, 1.0e+20, &r), 1.0e+20*exp(x) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, exp(-x)*exp(M_LN2)),       2.0, 1.0e-13, GSL_SUCCESS );

  TEST_SF(s,  gsl_sf_expm1_e, (-10.0, &r), exp(-10.0)-1.0          , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expm1_e, (-0.001, &r), -0.00099950016662500845, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expm1_e, (-1.0e-8, &r), -1.0e-08 + 0.5e-16      , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expm1_e, ( 1.0e-8, &r), 1.0e-08 + 0.5e-16      , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expm1_e, ( 0.001, &r), 0.0010005001667083417  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expm1_e, ( 10.0, &r), exp(10.0)-1.0           , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exprel_e, (-10.0, &r), 0.0999954600070237515 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_e, (-0.001, &r), 0.9995001666250084    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_e, (-1.0e-8, &r), 1.0 - 0.5e-08       , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_e, ( 1.0e-8, &r), 1.0 + 0.5e-08       , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_e, ( 0.001, &r), 1.0005001667083417   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_e, ( 10.0, &r), 2202.5465794806716517 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exprel_2_e, (-10.0, &r), 0.18000090799859524970 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_2_e, (-0.001, &r), 0.9996667499833361107  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_2_e, (-1.0e-8, &r), 0.9999999966666666750  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_2_e, ( 1.0e-8, &r), 1.0000000033333333417  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_2_e, ( 0.001, &r), 1.0003334166833361115  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_2_e, ( 10.0, &r), 440.3093158961343303   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exprel_n_e, (3, -1000.0, &r), 0.00299400600000000000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3, -100.0, &r), 0.02940600000000000000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3, -10.0, &r), 0.24599972760042142509 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3, -3.0, &r), 0.5444917625849191238  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3, -0.001, &r), 0.9997500499916678570  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3, -1.0e-8, &r), 0.9999999975000000050  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3,  1.0e-8, &r), 1.0000000025000000050  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3,  0.001, &r), 1.0002500500083345240  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3,  3.0, &r), 2.5745637607083706091  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3,  3.1, &r), 2.6772417068460206247  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3,  10.0, &r), 131.79279476884029910  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (3,  100.0, &r), 1.6128702850896812690e+38 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exprel_n_e, (50, -1000.0, &r), 0.04766231609253975959 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50, -100.0, &r), 0.3348247572345889317  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50, -10.0, &r), 0.8356287051853286482  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50, -3.0, &r), 0.9443881609152163615  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50, -1.0, &r), 0.980762245565660617   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50, -1.0e-8, &r), 1.0 -1.0e-8/51.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50,  1.0e-8, &r), 1.0 +1.0e-8/51.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50,  1.0, &r), 1.01999216583666790   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50,  3.0, &r), 1.0624205757460368307 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50,  48.0, &r), 7.499573876877194416  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50,  50.1, &r), 9.311803306230992272  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50,  100.0, &r), 8.175664432485807634e+07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (50,  500.0, &r), 4.806352370663185330e+146 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exprel_n_e, (500, -1000.0, &r), 0.3334815803127619256 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500, -100.0, &r), 0.8335646217536183909 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500, -10.0, &r), 0.9804297803131823066 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500, -3.0, &r), 0.9940475488850672997 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500, -1.0, &r), 0.9980079602383488808 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500, -1.0e-8, &r), 1.0 -1.0e-8/501.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500,  1.0e-8, &r), 1.0 +1.0e-8/501.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500,  1.0, &r), 1.0019999920160634252 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500,  3.0, &r), 1.0060240236632444934 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500,  48.0, &r), 1.1059355517981272174 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500,  100.0, &r), 1.2492221464878287204 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500,  500.0, &r), 28.363019877927630858 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500,  1000.0, &r), 2.4037563160335300322e+68 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_e, (500,  1600.0, &r), 7.899293535320607403e+226 , TEST_TOL0, GSL_SUCCESS);

  return s;
}


int check_expint(void)
{
  gsl_sf_result r;
  double y;
  int status = 0;
  int s;

  TEST_SF(s,  gsl_sf_expint_E1_e, (-1.0, &r), -1.8951178163559367555  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_e, (1.0e-10, &r), 22.448635265138923980  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_e, (1.0e-05, &r), 10.935719800043695615  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_e, (0.1, &r), 1.82292395841939066610 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_e, (1.0, &r), 0.21938393439552027368 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_e, (10.0, &r), 4.156968929685324277e-06  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_e, (50.0, &r), 3.783264029550459019e-24  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_e, (300.0, &r), 1.710384276804510115e-133 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_expint_E2_e, (-1.0, &r), 0.8231640121031084799  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_e, (1.0/4294967296.0, &r), 0.9999999947372139168  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_e, (1.0/65536.0, &r), 0.9998243233207178845  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_e, (0.1, &r), 0.7225450221940205066  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_e, (1.0, &r), 0.14849550677592204792 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_e, (10.0, &r), 3.830240465631608762e-06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_e, (50.0, &r), 3.711783318868827367e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_e, (300.0, &r), 1.7047391998483433998e-133 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_expint_Ei_e, (-1.0, &r), -0.21938393439552027368 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_Ei_e, (1.0/4294967296.0, &r), -21.603494112783886397  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_Ei_e, (1.0, &r), 1.8951178163559367555  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_Shi_e, (-1.0, &r), -1.0572508753757285146     , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_e, (1.0/4294967296.0, &r), 2.3283064365386962891e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_e, (1.0/65536.0, &r), 0.00001525878906269737298 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_e, (0.1, &r), 0.1000555722250569955 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_e, (1.0, &r), 1.0572508753757285146 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_e, (10.0, &r), 1246.1144901994233444 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_e, (50.0, &r), 5.292818448565845482e+19  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_e, (300.0, &r), 3.248241254044332895e+127 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_Chi_e, (-1.0, &r), 0.8378669409802082409 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_e, (1.0/4294967296.0, &r), -21.603494113016717041 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_e, (1.0/65536.0, &r), -10.513139223999384429 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_e, (1.0/8.0, &r), -1.4983170827635760646 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_e, (1.0, &r), 0.8378669409802082409 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_e, (10.0, &r), 1246.1144860424544147 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_e, (50.0, &r), 5.292818448565845482e+19  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_e, (300.0, &r), 3.248241254044332895e+127 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_expint_3_e, (1.0e-10, &r), 1.0e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_e, (1.0e-05, &r), 9.9999999999999975e-06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_e, (0.1, &r), 0.09997500714119079665122 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_e, (0.5, &r), 0.48491714311363971332427 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_e, (1.0, &r), 0.80751118213967145285833 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_e, (2.0, &r), 0.89295351429387631138208 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_e, (5.0, &r), 0.89297951156924921121856 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_e, (10.0, &r), 0.89297951156924921121856 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_e, (100.0, &r), 0.89297951156924921121856 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_Si_e, (-1.0, &r), -0.9460830703671830149 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_e, (1.0e-10, &r), 1.0e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_e, (1.0e-05, &r), 9.999999999944444444e-06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_e, (0.1, &r), 0.09994446110827695016   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_e, (1.0, &r), 0.9460830703671830149    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_e, (10.0, &r), 1.6583475942188740493    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_e, (50.0, &r), 1.5516170724859358947 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_e, (300.0, &r), 1.5708810882137495193 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_e, (1.0e+20, &r), 1.5707963267948966192 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_Ci_e, (1.0/4294967296.0, &r), -21.603494113016717041   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_e, (1.0/65536.0, &r), -10.513139224115799751   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_e, (1.0/8.0, &r), -1.5061295845296396649   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_e, (1.0, &r), 0.3374039229009681347   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_e, (10.0, &r), -0.04545643300445537263  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_e, (50.0, &r), -0.005628386324116305440 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_e, (300.0, &r), -0.003332199918592111780 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_e, (65536.0, &r), 0.000010560248837656279453 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_e, (4294967296.0, &r), -1.0756463261957757485e-10  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_e, (1099511627776.0, &r), -3.689865584710764214e-13   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_atanint_e, (1.0e-10, &r), 1.0e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_e, (1.0e-05, &r), 9.99999999988888888889e-06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_e, (0.1, &r), 0.09988928686033618404 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_e, (1.0, &r), 0.91596559417721901505 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_e, (2.0, &r), 1.57601540344632342236 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_e, (10.0, &r), 3.71678149306806859029 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_e, (50.0, &r), 6.16499047850274874222 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_e, (300.0, &r), 8.96281388924518959990 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_e, (1.0e+5, &r), 18.084471031038661920  , TEST_TOL0, GSL_SUCCESS);


  return s;
}


int check_fermidirac(void)
{
  int status = 0;
  int s = 0;

  TEST_SF(s, gsl_sf_fermi_dirac_m1_e, (-10.0, &r), 0.00004539786870243439450 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_m1_e, ( -1.0, &r), 0.26894142136999512075 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_m1_e, (  1.0, &r), 0.7310585786300048793  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_m1_e, ( 10.0, &r), 0.9999546021312975656  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_0_e, (-10.0, &r), 0.00004539889921686464677 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_0_e, ( -1.0, &r), 0.31326168751822283405 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_0_e, (  1.0, &r), 1.3132616875182228340  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_0_e, ( 10.0, &r), 10.000045398899216865  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_1_e, (-10.0, &r), 0.00004539941448447633524 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, ( -2.0, &r), 0.13101248471442377127 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, ( -1.0, &r), 0.3386479964034521798  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, ( -0.4, &r), 0.5825520806897909028  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, (  0.4, &r), 1.1423819861584355337  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, (  1.0, &r), 1.8062860704447742567  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, (  1.5, &r), 2.5581520872227806402  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, (  2.5, &r), 4.689474797599761667   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, ( 10.0, &r), 51.64488866743374196   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, ( 12.0, &r), 73.64492792264531092   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, ( 20.0, &r), 201.64493406478707282  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_e, ( 50.0, &r), 1251.6449340668482264  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_2_e, (-10.0, &r), 0.00004539967212174776662 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, ( -2.0, &r), 0.13313272938565030508 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, ( -1.0, &r), 0.3525648792978077590  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, ( -0.4, &r), 0.6229402647001272120  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, (  0.4, &r), 1.2915805581060844533  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, (  1.0, &r), 2.1641656128127008622  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, (  1.5, &r), 3.247184513920792475   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, (  2.5, &r), 6.797764392735056317   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, ( 10.0, &r), 183.11605273482105278  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, ( 12.0, &r), 307.73921494638635166  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, ( 20.0, &r), 1366.2320146723590157  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, ( 50.0, &r), 20915.580036675744655  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_e, (200.0, &r), 1.3336623201467029786e+06 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, (-10.0, &r), 0.00004539847236080549532 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, ( -2.0, &r), 0.12366562180120994266 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, ( -1.0, &r), 0.29402761761145122022 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, ( -0.4, &r), 0.4631755336886027800 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, (  0.4, &r), 0.7654084737661656915 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, (  1.0, &r), 1.0270571254743506890 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, (  1.5, &r), 1.2493233478527122008 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, (  2.5, &r), 1.6663128834358313625 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, ( 10.0, &r), 3.552779239536617160 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, ( 12.0, &r), 3.897268231925439359 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, ( 20.0, &r), 5.041018507535328603 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_e, ( 50.0, &r), 7.977530858581869960 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_half_e, (-10.0, &r), 0.00004539920105264132755 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, ( -2.0, &r), 0.12929851332007559106 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, ( -1.0, &r), 0.3277951592607115477 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, ( -0.4, &r), 0.5522452153690688947 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, (  0.4, &r), 1.0386797503389389277 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, (  1.0, &r), 1.5756407761513002308 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, (  1.5, &r), 2.1448608775831140360 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, (  2.5, &r), 3.606975377950373251  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, ( 10.0, &r), 24.084656964637653615 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, ( 12.0, &r), 31.540203287044242593 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, ( 20.0, &r), 67.49151222165892049  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_e, ( 50.0, &r), 266.09281252136259343 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, (-10.0, &r), 0.00004539956540456176333 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, ( -2.0, &r), 0.13224678225177236685 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, ( -1.0, &r), 0.3466747947990574170  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, ( -0.4, &r), 0.6056120213305040910  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, (  0.4, &r), 1.2258236403963668282  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, (  1.0, &r), 2.0022581487784644573  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, (  1.5, &r), 2.9277494127932173068  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, (  2.5, &r), 5.768879312210516582   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, ( 10.0, &r), 101.00510084332600020  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, ( 12.0, &r), 156.51518642795728036  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, ( 20.0, &r), 546.5630100657601959   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_e, ( 50.0, &r), 5332.353566687145552   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (3,  -2.0, &r), 0.1342199155038680215 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (3,   0.0, &r), 0.9470328294972459176 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (3,   0.1, &r), 1.0414170610956165759 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (3,   1.0, &r), 2.3982260822489407070 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (3,   3.0, &r), 12.621635313399690724 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (3, 100.0, &r), 4.174893231066566793e+06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (3, 500.0, &r), 2.604372285319088354e+09 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (5,  -2.0, &r), 0.13505242246823676478 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (5,   0.0, &r), 0.9855510912974351041  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (5,   0.1, &r), 1.0876519750101492782  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (5,   1.0, &r), 2.6222337848692390539  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (5,   3.0, &r), 17.008801618012113022  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (5, 100.0, &r), 1.3957522531334869874e+09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (5, 500.0, &r), 2.1705672808114817955e+13 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (7,  -2.0, &r), 0.1352641105671255851 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (7,   0.0, &r), 0.9962330018526478992 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (7,   0.1, &r), 1.1005861815180315485 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (7,   1.0, &r), 2.6918878172003129203 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (7,   3.0, &r), 19.033338976999367642 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (7,  10.0, &r), 5654.530932873610014  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (7,  50.0, &r), 1.005005069985066278e+09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (7, 500.0, &r), 9.691690268341569514e+16 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (9,  -2.0, &r), 0.1353174385330242691 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (9,   0.0, &r), 0.9990395075982715656 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (9,   0.1, &r), 1.1039997234712941212 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (9,   1.0, &r), 2.7113648898129249947 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (9,   3.0, &r), 19.768544008138602223 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (9,  10.0, &r), 10388.990167312912478 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (9,  50.0, &r), 2.85466960802601649e+10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (9, 500.0, &r), 2.69273849842695876e+20 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (10,  -2.0, &r), 0.13532635396712288092 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (10,   0.0, &r), 0.9995171434980607541 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (10,   0.1, &r), 1.1045818238852612296 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (10,   1.0, &r), 2.7147765350346120647 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (10,   3.0, &r), 19.917151938411675171 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (10,  10.0, &r), 12790.918595516495955 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (10,  50.0, &r), 1.3147703201869657654e+11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (10, 500.0, &r), 1.2241331244469204398e+22 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (11,  -2.0, &r), 0.1353308162894847149 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (11,   0.0, &r), 0.9997576851438581909 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (11,   0.1, &r), 1.1048751811565850418 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (11,   1.0, &r), 2.7165128749007313436 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (11,   3.0, &r), 19.997483022044603065 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (11,  10.0, &r), 14987.996005901818036 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (11,  50.0, &r), 5.558322924078990628e+11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (11, 500.0, &r), 5.101293089606198280e+23 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (20,  -2.0, &r), 0.13533527450327238373 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (20,   0.0, &r), 0.9999995232582155428  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (20,   0.1, &r), 1.1051703357941368203  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (20,   1.0, &r), 2.7182783069905721654  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (20,   3.0, &r), 20.085345296028242734  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (20,  10.0, &r), 21898.072920149606475  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (20,  50.0, &r), 1.236873256595717618e+16 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_e, (20, 500.0, &r), 9.358938204369557277e+36 , TEST_TOL0, GSL_SUCCESS);


  return s;
}


int check_gamma(void)
{
  gsl_sf_result r;
  double y;
  double zr, zi, lg_r, lg_i, lg, sgn;
  int status = 0;
  int s;

  TEST_SF(s,  gsl_sf_lngamma_e, (-0.1, &r), 2.368961332728788655 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lngamma_e, (-1.0/256.0, &r), 5.547444766967471595  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lngamma_e, (1.0e-08, &r), 18.420680738180208905 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lngamma_e, (0.1, &r), 2.252712651734205 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lngamma_e, (1.0 + 1.0/256.0, &r), -0.0022422226599611501448 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lngamma_e, (2.0 + 1.0/256.0, &r), 0.0016564177556961728692 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lngamma_e, (100.0, &r), 359.1342053695753 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lngamma_e, (-1.0-1.0/65536.0, &r), 11.090348438090047844 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lngamma_e, (-1.0-1.0/268435456.0, &r), 19.408121054103474300 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lngamma_e, (-100.5, &r), -364.9009683094273518 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lngamma_e, (-100-1.0/65536.0, &r), -352.6490910117097874 , TEST_TOL0, GSL_SUCCESS);

  s = 0;
  gsl_sf_lngamma_sgn_impl(0.7, &lg, &sgn);
  s += ( frac_diff( lg,  0.26086724653166651439 ) > TOL );
  s += ( sgn != 1.0 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(0.7)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(0.1, &lg, &sgn);
  s += ( frac_diff( lg,  2.2527126517342059599 ) > TOL );
  s += ( sgn != 1.0 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(0.1)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-0.1, &lg, &sgn);
  s += ( frac_diff( lg,   2.368961332728788655 ) > TOL );
  s += ( sgn != -1.0 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-0.1)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-1.0-1.0/65536.0, &lg, &sgn);
  s += ( frac_diff( lg,  11.090348438090047844 ) > TOL );
  s += ( frac_diff( sgn, 1.0 ) > TOL );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-1.0-1.0/65536.0)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-2.0-1.0/256.0, &lg, &sgn);
  s += ( frac_diff( lg, 4.848447725860607213 ) > TOL );
  s += ( frac_diff( sgn, -1.0 ) > TOL );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-2.0-1.0/256.0)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-2.0-1.0/65536.0, &lg, &sgn);
  s += ( frac_diff( lg, 10.397193628164674967 ) > TOL );
  s += ( frac_diff( sgn, -1.0 ) > TOL );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-2.0-1.0/65536.0)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-3.0-1.0/8.0, &lg, &sgn);
  s += ( frac_diff( lg, 0.15431112768404182427 ) > 2.0*TOL );
  s += ( sgn != 1.0 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-3.0-1.0/8.0)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-100.5, &lg, &sgn);
  s += ( frac_diff( lg, -364.9009683094273518 ) > TOL );
  s += ( sgn != -1.0 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-100.5)");
  status += s;


  TEST_SF(s,  gsl_sf_gamma_e, (1.0 + 1.0/4096.0, &r), 0.9998591371459403421 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gamma_e, (1.0 + 1.0/32.0, &r), 0.9829010992836269148 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gamma_e, (2.0 + 1.0/256.0, &r), 1.0016577903733583299 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gamma_e, (9.0, &r), 40320.0                   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gamma_e, (10.0, &r), 362880.0                  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gamma_e, (100.0, &r), 9.332621544394415268e+155 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gamma_e, (170.0, &r), 4.269068009004705275e+304 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gamma_e, (-10.5, &r), -2.640121820547716316e-07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gamma_e, (-11.25, &r), 6.027393816261931672e-08  , TEST_TOL0, GSL_SUCCESS); /* exp()... not my fault */
  TEST_SF(s,  gsl_sf_gamma_e, (-1.0+1.0/65536.0, &r), -65536.42280587818970 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_gammastar_e, (1.0e-08, &r), 3989.423555759890865  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammastar_e, (1.0e-05, &r), 126.17168469882690233 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammastar_e, (0.001, &r), 12.708492464364073506 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammastar_e, (1.5, &r), 1.0563442442685598666 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammastar_e, (3.0, &r), 1.0280645179187893045 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammastar_e, (9.0, &r), 1.0092984264218189715 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammastar_e, (11.0, &r), 1.0076024283104962850 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammastar_e, (100.0, &r), 1.0008336778720121418 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammastar_e, (1.0e+05, &r), 1.0000008333336805529 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammastar_e, (1.0e+20, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_gammainv_e, (10.0, &r), 1.0/362880.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammainv_e, (100.0, &r), 1.0715102881254669232e-156    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammainv_e, (-10.5, &r), -1.0/2.640121820547716316e-07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammainv_e, (-11.25, &r), 1.0/6.027393816261931672e-08 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gammainv_e, (-1.0+1.0/65536.0, &r), -1.0/65536.42280587818970 , TEST_TOL0, GSL_SUCCESS);


  s = 0;
  zr = 5.0;
  zi = 2.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, 2.7487017561338026749 ) > TOL );
  s += ( frac_diff( lg_i, 3.0738434100497007915 ) > TOL );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(5 + 2 I)");
  status += s;

  s = 0;
  zr = 100.0;
  zi = 100.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, 315.07804459949331323 ) > TOL );
  s += ( frac_diff( lg_i, 2.0821801804113110099 ) > TOL );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(100 + 100 I)");
  printf("%22.18g\n %22.18g\n", lg_r, lg_i);
  status += s;

  s = 0;
  zr =   100.0;
  zi = -1000.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, -882.3920483010362817000 ) > TOL );
  s += ( frac_diff( lg_i,   -2.1169293725678813270 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(100 - 1000 I)");
  printf("%22.18g\n %22.18g\n", lg_r, lg_i);
  status += s;

  s = 0;
  zr = -100.0;
  zi =   -1.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, -365.0362469529239516000 ) > TOL );
  s += ( frac_diff( lg_i,   -3.0393820262864361140 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(-1000 - I)");
  printf("%22.18g\n %22.18g\n", lg_r, lg_i);
  status += s;

  TEST_SF(s,  gsl_sf_taylorcoeff_e, (10,   1.0/1048576.0, &r), 1.7148961854776073928e-67  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_taylorcoeff_e, (10,   1.0/1024.0, &r), 2.1738891788497900281e-37  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_taylorcoeff_e, (10,   1.0, &r), 2.7557319223985890653e-07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_taylorcoeff_e, (10,   5.0, &r), 2.6911444554673721340      , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_taylorcoeff_e, (10,   500.0, &r), 2.6911444554673721340e+20  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_taylorcoeff_e, (100,  100.0, &r), 1.0715102881254669232e+42  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_taylorcoeff_e, (1000, 200.0, &r), 2.6628790558154746898e-267 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_taylorcoeff_e, (1000, 500.0, &r), 2.3193170139740855074e+131 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,   gsl_sf_fact_e, (0, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,   gsl_sf_fact_e, (1, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,   gsl_sf_fact_e, (7, &r), 5040.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_fact_e, (33, &r), 8.683317618811886496e+36 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,   gsl_sf_doublefact_e, (0, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,   gsl_sf_doublefact_e, (1, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,   gsl_sf_doublefact_e, (7, &r), 105.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_doublefact_e, (33, &r), 6.332659870762850625e+18 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_lnfact_e, (0, &r), 0.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnfact_e, (1, &r), 0.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnfact_e, (7, &r), 8.525161361065414300 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnfact_e, (33, &r), 85.05446701758151741 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_lndoublefact_e, (0, &r), 0.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lndoublefact_e, (7, &r), 4.653960350157523371  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lndoublefact_e, (33, &r), 43.292252022541719660 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lndoublefact_e, (34, &r), 45.288575519655959140 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lndoublefact_e, (1034, &r), 3075.6383796271197707 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lndoublefact_e, (1035, &r), 3078.8839081731809169 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_lnchoose_e, (7,3, &r), 3.555348061489413680 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnchoose_e, (5,2, &r), 2.302585092994045684 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_choose_e, (7,3, &r), 35.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_choose_e, (5,2, &r), 10.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_choose_e, (500,200, &r), 5.054949849935532221e+144 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_lnpoch_e, (5, 1.0/65536.0, &r), 0.000022981557571259389129 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnpoch_e, (5, 1.0/256.0, &r), 0.005884960217985189004    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnpoch_e, (7,3, &r), 6.222576268071368616 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnpoch_e, (5,2, &r), 3.401197381662155375 , TEST_TOL0, GSL_SUCCESS);

  s = 0;
  gsl_sf_lnpoch_sgn_impl(-4.5, 0.25, &y, &sgn);
  s += ( frac_diff(y, 0.7430116475119920117 ) > TOL );
  s += ( sgn != 1.0 );
  gsl_test(s, "  gsl_sf_lnpoch_sgn_impl(-4.5, 0.25)");
  status += s;

  s = 0;
  gsl_sf_lnpoch_sgn_impl(-4.5, 1.25, &y, &sgn);
  s += ( frac_diff(y, 2.1899306304483174731 ) > TOL );
  s += ( sgn != -1.0 );
  gsl_test(s, "  gsl_sf_lnpoch_sgn_impl(-4.5, 1.25)");
  status += s;

  TEST_SF(s,  gsl_sf_poch_e, (7,3, &r), 504.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (5,2, &r), 30.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (5,1.0/256.0, &r), 1.0059023106151364982 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_pochrel_e, (7,3, &r), 503.0/3.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (5,2, &r), 29.0/2.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (5,0.01, &r), 1.5186393661368275330  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (-5.5,0.01, &r), 1.8584945633829063516  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (-5.5,-1.0/8.0, &r), 1.0883319303552135488  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (-5.5,-1.0/256.0, &r), 1.7678268037726177453  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (-5.5,-11.0, &r), 0.09090909090939652475 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_gamma_inc_P_e, (0.001, 0.001), 0.9936876467088602902, TEST_TOL0, GSL_SUCCESS) ;
  TEST_SF(s, gsl_sf_gamma_inc_P_e, (0.001, 1.0), 0.9997803916424144436, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_P_e, (0.001, 10.0), 0.9999999958306921828, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_P_e, (1.0, 0.001), 0.0009995001666250083319, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_P_e, (1.0, 1.01), 0.6357810204284766802, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_P_e, (1.0, 10.0), 0.9999546000702375151, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_P_e, (10.0, 10.01), 0.5433207586693410570, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_P_e, (10.0, 20.0), 0.9950045876916924128, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_P_e, (1000.0, 1000.1), 0.5054666401440661753, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_P_e, (1000.0, 2000.0), 1.0, TEST_TOL0, GSL_SUCCESS);
  
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (0.001, 0.001), 0.006312353291139709793, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (0.001, 1.0), 0.00021960835758555639171, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (0.001, 2.0), 0.00004897691783098147880, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (0.001, 5.0), 1.1509813397308608541e-06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (1.0, 0.001), 0.9990004998333749917, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (1.0, 1.01), 0.3642189795715233198, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (1.0, 10.0), 0.00004539992976248485154, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (10.0, 10.01), 0.4566792413306589430, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (10.0, 100.0), 1.1253473960842733885e-31, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (1000.0, 1000.1), 0.4945333598559338247, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_gamma_inc_Q_e, (1000.0, 2000.0), 6.847349459614753180e-136, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lnbeta_e, (1.0e-8, 1.0e-8, &r),  19.113827924512310617 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (1.0e-8, 0.01, &r),  18.420681743788563403 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (1.0e-8, 1.0, &r),  18.420680743952365472 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (1.0e-8, 10.0, &r),  18.420680715662683009 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (1.0e-8, 1000.0, &r),  18.420680669107656949 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (0.1, 0.1, &r), 2.9813614810376273949 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (0.1, 1.0, &r),  2.3025850929940456840 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (0.1, 100.0, &r),  1.7926462324527931217 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (0.1, 1000, &r),  1.5619821298353164928 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (1.0, 1.0001, &r),  -0.00009999500033330833533 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (1.0, 1.01, &r),  -0.009950330853168082848 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (1.0, 1000.0, &r),  -6.907755278982137052 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (100.0, 100.0, &r),  -139.66525908670663927 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (100.0, 1000.0, &r),  -336.4348576477366051 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_lnbeta_e, (100.0, 1.0e+8, &r),  -1482.9339185256447309 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_e, (1.0,   1.0, &r), 1.0                   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_e, (1.0, 1.001, &r), 0.9990009990009990010 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_e, (1.0,   5.0, &r), 0.2                   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_e, (1.0,  100.0, &r), 0.01                  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_e, (10.0, 100.0, &r), 2.3455339739604649879e-15 , TEST_TOL0, GSL_SUCCESS);

  return s;
}

int check_gegen(void)
{
  gsl_sf_result r;
  double ga[100];
  double y;
  int status = 0;
  int s = 0;

  TEST_SF(s,  gsl_sf_gegenpoly_1_e, (-0.2,   1.0, &r), -0.4 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_1_e, ( 0.0,   1.0, &r), 2.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_1_e, ( 1.0,   1.0, &r), 2.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_1_e, ( 1.0,   0.5, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_1_e, ( 5.0,   1.0, &r), 10.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_1_e, ( 100.0, 0.5, &r), 100.0 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_gegenpoly_2_e, (-0.2,   0.5, &r), 0.12 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_2_e, ( 0.0,   1.0, &r), 1.00 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_2_e, ( 1.0,   1.0, &r), 3.00 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_2_e, ( 1.0,   0.1, &r), -0.96 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_2_e, ( 5.0,   1.0, &r), 55.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_2_e, ( 100.0, 0.5, &r), 4950.0 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_gegenpoly_3_e, (-0.2,   0.5, &r), 0.112 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_3_e, ( 0.0,   1.0, &r), -2.0/3.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_3_e, ( 1.0,   1.0, &r), 4.000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_3_e, ( 1.0,   0.1, &r), -0.392 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_3_e, ( 5.0,   1.0, &r), 220.000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_3_e, ( 100.0, 0.5, &r), 161600.000 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_gegenpoly_n_e, (1,       1.0, 1.0, &r), 2.000		    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_e, (10,      1.0, 1.0, &r), 11.000		    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_e, (10,      1.0, 0.1, &r), -0.4542309376	    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_e, (10,      5.0, 1.0, &r), 9.23780e+4  	    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_e, (10,    100.0, 0.5, &r), 1.5729338392690000e+13  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_e, (1000,  100.0, 1.0, &r), 3.3353666135627322e+232 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_e, (100,  2000.0, 1.0, &r), 5.8753432034937579e+202 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_e, (103,   207.0, 2.0, &r), 1.4210272202235983e+145 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_e, (103,    -0.4, 0.3, &r), -1.64527498094522e-04    , TEST_TOL0, GSL_SUCCESS);

  s = 0;
  gsl_sf_gegenpoly_array_impl(99, 5.0, 1.0, ga);
  s += ( frac_diff( ga[1],     10.0    ) > TOL );
  s += ( frac_diff( ga[10], 9.23780e+4 ) > TOL );
  gsl_test(s, "  gsl_sf_gegenpoly_array_impl");
  status += s;

  return status;
}


int check_hyperg(void)
{
  double y;
  int status = 0;
  int s;


  /* 0F1 */

  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1, 0.5, &r), 1.5660829297563505373 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (5, 0.5, &r), 1.1042674404828684574 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (100, 30, &r), 1.3492598639485110176 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-0.5, 3, &r), -39.29137997543434276  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-100.5, 50, &r), 0.6087930289227538496 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1, -5.0, &r), -0.3268752818235339109 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-0.5, -5.0, &r), -4.581634759005381184  , TEST_TOL0, GSL_SUCCESS);


  /* 1F1 for integer parameters */

  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (1, 1, 0.5, &r), 1.6487212707001281468 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (1, 2, 500.0, &r), 2.8071844357056748215e+214 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (1, 2, -500.0, &r), 0.002 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (8, 1, 0.5, &r), 13.108875178030540372 , TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, 1.0, &r),  131.63017574352619931 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, 10.0, &r), 8.514625476546280796e+09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, 100.0, &r),  1.5671363646800353320e+56 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 20, 1.0, &r),  1.6585618002669675465 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 20, 10.0, &r),  265.26686430340188871 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 20, 100.0, &r), 3.640477355063227129e+34 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 1.0, &r),  1.1056660194025527099 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 10.0, &r),  2.8491063634727594206 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 40.0, &r),  133.85880835831230986 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 80.0, &r),  310361.16228011433406 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 100.0, &r),  8.032171336754168282e+07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 500.0, &r),  7.633961202528731426e+123 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, 1.0, &r),  6.892842729046469965e+07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, 10.0, &r),  2.4175917112200409098e+28 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, 100.0, &r),  1.9303216896309102993e+110 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 1.0, &r),  1.6497469106162459226 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 10.0, &r),  157.93286197349321981 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 100.0, &r),  2.1819577501255075240e+24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 400.0, &r),  3.728975529926573300e+119 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 10.0, &r),  12.473087623658878813 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 100.0, &r),  9.071230376818550241e+11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 150.0, &r),  7.160949515742170775e+18 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 200.0, &r),   2.7406690412731576823e+26 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 300.0, &r),  6.175110613473276193e+43 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 400.0, &r),  1.1807417662711371440e+64 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 600.0, &r),  2.4076076354888886030e+112 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, -1.0, &r),  0.11394854824644542810 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, -10.0, &r),  0.0006715506365396127863 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, -100.0, &r),  -4.208138537480269868e-32 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 50, -1.0, &r),  0.8200061961 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -10.0, &r),   0.38378859043466243 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -100.0, &r),  0.0008460143401464189061 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -500.0, &r),  1.1090822141973655929e-08 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -10000.0, &r),  5.173783508088272292e-21 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (50, 1, -90.0, &r),  -1.6624258547648311554e-21 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (50, 1, -100.0, &r),  4.069661775122048204e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (50, 1, -110.0, &r),  1.0072444993946236025e-25 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 10, -100.0, &r),  -2.7819353611733941962e-37 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, -90.0, &r),  7.501705041159802854e-22 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, -100.0, &r),  6.305128893152291187e-25 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, -110.0, &r),  -7.007122115422439755e-26 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 10, -100.0, &r),  -2.7819353611733941962e-37 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 50, -1.0, &r),  0.016087060191732290813 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 50, -300.0, &r),  -4.294975979706421471e-121 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -1.0, &r),  0.13397521083325179687 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -10.0, &r),  5.835134393749807387e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -100.0, &r),  4.888460453078914804e-74 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -500.0, &r),  -1.4478509059582015053e-195 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, 1, 2.0, &r),  -1.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, -2, 2.0, &r),  2.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-2, -3, 2.0, &r),  3.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, 1, 1.0, &r),  0.4189459325396825397 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, 1, 10.0, &r),  27.984126984126984127 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, 1, 100.0, &r),  9.051283795429571429e+12 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, 20, 1.0, &r),  0.0020203016320697069566 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, 1.0, &r),  1.6379141878548080173 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, 10.0, &r),  78.65202404521289970 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, 100.0, &r),  4.416169713262624315e+08 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, 1.0, &r),  1.1046713999681950919 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, 10.0, &r),  2.6035952191039006838 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, 100.0, &r),  1151.6852040836932392 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, 1.0, &r),  1.6476859702535324743 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, 10.0, &r),  139.38026829540687270 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, 100.0, &r),  1.1669433576237933752e+19 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, -1.0, &r),  0.6025549561148035735 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, -10.0, &r),  0.00357079636732993491 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, -100.0, &r),  1.64284868563391159e-35 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, -1.0, &r),  0.90442397250313899 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, -10.0, &r),  0.35061515251367215, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, -100.0, &r),  8.19512187960476424e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, -1.0, &r),  0.6061497939628952629 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, -10.0, &r),  0.0063278543908877674 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, -100.0, &r),  4.34111795007336552e-25 , TEST_TOL0, GSL_SUCCESS);


  /* 1F1 */

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, 1, &r), 2.0300784692787049755 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, 10, &r),  6172.859561078406855 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, 100, &r),  2.3822817898485692114e+42 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, 500, &r),  5.562895351723513581e+215 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.5, 2.5, 1, &r), 1.8834451238277954398 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.5, 2.5, 10, &r),  3128.7352996840916381 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, 1, &r),  110.17623733873889579 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, 10, &r),  6.146657975268385438e+09 , TEST_TOL0, GSL_SUCCESS);

  s = 0;
  gsl_sf_hyperg_1F1_impl(10, 1.1, 100, &y);
  s += ( frac_diff(y, 9.331833897230312331e+55 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(10, 1.1, 100)");
  status += s;

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, 500, &r),  4.519403368795715843e+235 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 50.1, 2, &r),  1.5001295507968071788 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 50.1, 10, &r),  8.713385849265044908 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 50.1, 100, &r),  5.909423932273380330e+18 , TEST_TOL0, GSL_SUCCESS);

  s = 0;
  gsl_sf_hyperg_1F1_impl(10, 50.1, 500, &y);
  s += ( frac_diff(y, 9.740060618457198900e+165 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(10, 50.1, 500)");
  status += s;

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, 1, &r),  5.183531067116809033e+07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, 10, &r),  1.6032649110096979462e+28 , TEST_TOL0, GSL_SUCCESS);

  s = 0;
  gsl_sf_hyperg_1F1_impl(100, 1.1, 100, &y);
  s += ( frac_diff(y, 1.1045151213192280064e+110 ) > 1.e-11 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(100, 1.1, 100)");
  status += s;


  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 50.1, 1, &r),  7.222953133216603757 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 50.1, 10, &r),  1.0998696410887171538e+08 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 50.1, 100, &r),  7.235304862322283251e+63 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, -1, &r), 0.5380795069127684191 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, -10, &r),  0.05303758099290164485 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, -100, &r), 0.005025384718759852803 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, -500, &r), 0.0010010030151059555322 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.1, -500, &r), 0.00020036137599690208265 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, -1, &r), 0.07227645648935938168 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, -10, &r),  0.0003192415409695588126 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, -100, &r),  -8.293425316123158950e-16 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, -500, &r),  -3.400379216707701408e-23 , TEST_TOL0, GSL_SUCCESS);

  s = 0;
  gsl_sf_hyperg_1F1_impl(50, 1.1, -90.0, &y);
  s += ( frac_diff(y, -7.843129411802921440e-22 ) > 1.e-08 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(50, 1.1, -90.0)");
  printf("%22.28g\n", y);
  status += s;

  s = 0;
  gsl_sf_hyperg_1F1_impl(50, 1.1, -100.0, &y);
  s += ( frac_diff(y, 4.632883869540640460e-24 ) > 1.e-08 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(50, 1.1, -100.0)");
  printf("%22.28g\n", y);
  status += s;

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (50, 1.1, -110.0, &r),  5.642684651305310023e-26 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -1, &r),  0.0811637344096042096 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -10, &r),  0.00025945610092231574387 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -50, &r),  2.4284830988994084452e-13 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -90, &r),  2.4468224638378426461e-22 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -99, &r),  1.0507096272617608461e-23 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -100, &r),  1.8315497474210138602e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -101, &r),  -2.3916306291344452490e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -110, &r),  -4.517581986037732280e-26 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 10.1, -220, &r),  -4.296130300021696573e-64 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -10.1, 10.0, &r),  10959.603204633058116 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -10.1, 1000.0, &r),  2.0942691895502242831e+23 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -100.1, 10.0, &r),  2.6012036337980078062 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1000, -1000.1, 10.0, &r),  22004.341698908631636 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1000, -1000.1, 200.0, &r),  7.066514294896245043e+86 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-8.1, -10.1, -10.0, &r),  0.00018469685276347199258 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-8.1, -1000.1, -10.0, &r),  0.9218280185080036020 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -5.1, 1, &r),  16.936141866089601635 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -5.1, 10, &r),  771534.0349543820541 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -5.1, 100, &r),  2.2733956505084964469e+17 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, -1, &r),  0.13854540373629275583 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, -10, &r),  -9.142260314353376284e+19 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, -100, &r),  -1.7437371339223929259e+87 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, 1, &r),  7.516831748170351173 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, 10, &r),  1.0551632286359671976e+11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, 50, &r),  -7.564755600940346649e+36 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, 100, &r),  4.218776962675977e+55 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10.5, -8.1, 0.1, &r),  1.1387201443786421724 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10.5, -11.1, 1, &r),  2.5682766147138452362 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100.5, -80.1, 10, &r),  355145.4517305220603 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100.5, -102.1, 10, &r),  18678.558725244365016 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100.5, -500.1, 10, &r),  7.342209011101454 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100.5, -500.1, 100, &r),  1.2077443075367177662e+8 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-500.5, -80.1, 2, &r),  774057.8541325341699 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, 1, &r),  -2.1213846338338567395e+12, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, 10, &r),  -6.624849346145112398e+39 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, 100, &r),  -1.2413466759089171904e+129 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, -1, &r),  34456.29405305551691 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, -10, &r),  -7.809224251467710833e+07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, -100, &r),   -5.214065452753988395e-07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, 1, &r),  0.21519810496314438414 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, 10, &r),  8.196123715597869948 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, 100, &r),  -1.4612966715976530293e+20 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 20.1, 1, &r),  0.0021267655527278456412 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 20.1, 10, &r),   2.0908665169032186979e-11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 20.1, 100, &r),  -0.04159447537001340412 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, -1, &r),  2.1214770215694685282e+07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, -10, &r),  1.0258848879387572642e+24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, -100, &r),  1.1811367147091759910e+67 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 50.1, -1, &r),  6.965259317271427390 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 50.1, -10, &r),  1.0690052487716998389e+07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 50.1, -100, &r),  6.889644435777096248e+36 , TEST_TOL0, GSL_SUCCESS);


  /* U for integer parameters */

  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 0.0001, &r),  8.634088070212725330 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 0.01, &r),  4.078511443456425847 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 0.5, &r),  0.9229106324837304688 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 2.0, &r),  0.3613286168882225847 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 100, &r),  0.009901942286733018406 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 1000, &r),  0.0009990019940238807150 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 0.01, &r),  7.272361203006010000e+16 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 1, &r),  1957.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 5, &r),  1.042496 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 8, &r),  0.3207168579101562500 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 50, &r),  0.022660399001600000000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 100, &r),  0.010631236727200000000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 1000, &r),  0.0010060301203607207200 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 20, 1, &r),  1.7403456103284421000e+16 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 20, 20, &r),  0.22597813610531052969 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 50, 1, &r),  3.374452117521520758e+61 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 50, 50, &r),  0.15394136814987651785 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 0.1, &r),  1.0418325171990852858e+253 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 1, &r),  2.5624945006073464385e+154 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 50, &r),  3.0978624160896431391e+07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 100, &r),  0.11323192555773717475 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 200, &r),  0.009715680951406713589 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 1000, &r),  0.0011085142546061528661 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1000, 2000, &r),  0.0009970168547036318206 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -1, 1, &r),  0.29817368116159703717 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -1, 10, &r),  0.07816669698940409380 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 1, &r),  0.08271753756946041959 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 5, &r),  0.06127757419425055261 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 10, &r),  0.04656199948873187212 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 20, &r),  0.031606421847946077709 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 0.01, &r),  0.009900000099999796950 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 1, &r),  0.009802970197050404429 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 10, &r),  0.009001648897173103447 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 20, &r),  0.008253126487166557546 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 50, &r),  0.006607993916432051008 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 90, &r),  0.005222713769726871937 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 110, &r),  0.004727658137692606210 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -1000, 1, &r),  0.0009980029970019970050 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -1000, 1010, &r),  0.0004971408839859245170 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 0.001, &r),  0.0007505359326875706975 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 0.5, &r),  6.449509938973479986e-06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 8, &r),  6.190694573035761284e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 20, &r),  3.647213845460374016e-12 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 8, 1, &r),  0.12289755012652317578 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 8, 10, &r),  5.687710359507564272e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 8, 20, &r),  2.8175404594901039724e-11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 0.01, &r),  1.0099979491941914867e+196 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 0.1, &r),  1.0090713562719862833e+97 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 1, &r),  0.009998990209084729106 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 20, &r),  1.3239363905866130603e-131 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 0.01, &r),  3.274012540759009536e+06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 1, &r),  1.5202710000000000000e+06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 10, &r),  1.0154880000000000000e+08 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 100, &r),  3.284529863685482880e+19 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 10, 1, &r),  1.1043089864100000000e+11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 1, &r),  1.3991152402448957897e+20 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 10, &r),  5.364469916567136000e+19 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 100, &r),  3.909797568000000000e+12 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 500, &r),  8.082625576697984130e+25 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 1, 0.01, &r),  1.6973422555823855798e+64 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 1, 1, &r),  7.086160198304780325e+63 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 1, 10, &r),  5.332862895528712200e+65 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 10, 1, &r),  -7.106713471565790573e+71 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 100, 1, &r),  2.4661377199407186476e+104 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 10, 10, &r),  5.687538583671241287e+68 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 100, 10, &r),  1.7880761664553373445e+102 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1, 0.01, &r),  4.185245354032917715e+137 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1, 0.1, &r),  2.4234043408007841358e+137 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1, 10, &r),  -1.8987677149221888807e+139 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 10, 10, &r),  -5.682999988842066677e+143 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 100, 10, &r),  2.3410029853990624280e+189 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1000, 10, &r),  1.9799451517572225316e+271 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -1, 10, &r),  -9.083195466262584149e+64 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -10, 10, &r),  -1.4418257327071634407e+62 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -100, 0.01, &r),  3.0838993811468983931e+93 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -100, 10, &r),  4.014552630378340665e+95 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-100, -100, 10, &r),  2.0556466922347982030e+162 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-100, -200, 10, &r),  1.1778399522973555582e+219 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-100, -200, 100, &r),  9.861313408898201873e+235 , TEST_TOL0, GSL_SUCCESS);


  /* U */

  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 0.0001, 0.0001, &r), 1.0000576350699863577 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 0.0001, 1.0, &r), 0.9999403679233247536 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 0.0001, 100.0, &r), 0.9995385992657260887 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1, 0.0001, &r), 1.0009210608660065989 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1.0, 1.0, &r), 0.9999999925484179084 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 10, 1, &r), 13.567851006281412726 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 10, 5, &r), 1.0006265020064596364 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 10, 10, &r), 0.9999244381454633265 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 100, 1, &r),  2.5890615708804247881e+150 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 100, 10, &r),  2.3127845417739661466e+55 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 100, 50, &r), 6402.818715083582554 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 100, 98, &r), 0.9998517867411840044 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1000, 300, &r),  2.5389557274938010716e+213 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1000, 999, &r), 0.9997195294193261604 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1000, 1100, &r),  0.9995342990014584713 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.5, 1000, 300, &r), 1.1977955438214207486e+217 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.5, 1000, 800, &r), 9.103916020464797207e+08 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.5, 1000, 998, &r), 0.21970269691801966806 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.5, 0.5, 1.0, &r), 0.7578721561413121060 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 0.0001, 0.0001, &r), 0.9992361337764090785 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 0.0001, 1, &r), 0.4036664068111504538 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 0.0001, 100, &r), 0.009805780851264329587 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 1.2, 2.0, &r), 0.3835044780075602550 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, -0.0001, 1, &r), 0.4036388693605999482 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (8, 10.5, 1, &r),  27.981926466707438538 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (8, 10.5, 10, &r),  2.4370135607662056809e-8 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (8, 10.5, 100, &r),  1.1226567526311488330e-16 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (10, -2.5, 10, &r),  6.734690720346560349e-14 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (10, 2.5, 10, &r),  6.787780794037971638e-13 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (10, 2.5, 50, &r),  2.4098720076596087125e-18 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 1, &r),  -3.990841457734147e+6 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 10, &r),  1.307472052129343e+8 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 50, &r),  3.661978424114088e+16 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 90, &r),  8.09469542130868e+19 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 99, &r),  2.546328328942063e+20 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 100, &r),  2.870463201832814e+20 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 200, &r),  8.05143453769373e+23 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 0.1, &r),  -3.043016255306515e+20 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 1, &r),  -3.194745265896115e+12 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 4, &r),  -6.764203430361954e+07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 10, &r),  -2.067399425480545e+09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 50, &r),  4.661837330822824e+14 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 10, &r),  -6.805460513724838e+66 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 50, &r),  -2.081052558162805e+18 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 80, &r),  2.034113191014443e+14 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 100, &r),  6.85047268436107e+13 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 200, &r),  1.430815706105649e+20 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-19.5, 82.1, 10, &r),  5.464313196201917432e+60 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 10, &r),  -5.5740216266953e+126 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 40, &r),  5.937463786613894e+91 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 50, &r),  -1.631898534447233e+89 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 70, &r),  3.249026971618851e+84 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 100, &r),  1.003401902126641e+85 , TEST_TOL0, GSL_SUCCESS);


  /* 2F1 */

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1, 1, 1, 0.5, &r), 2.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, 8, 1, 0.5, &r), 12451584.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, -8, 1, 0.5, &r), 0.13671875 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, -8.1, 1, 0.5, &r), 0.14147385378899930422 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, -8, 1, -0.5, &r), 4945.136718750000000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, -8, -5.5, 0.5, &y, &r),  -906.6363636363636364 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, -8, -5.5, -0.5, &r), 24565.363636363636364 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, 8, 1, -0.5, &r), -0.006476312098196747669 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, 8, 5, 0.5, &r), 4205.714285714285714 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, 8, 5, -0.5, &r), 0.0028489656290296436616 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, 1, 0.99, &r), 1.2363536673577259280e+38  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -1.5, 0.99, &r), 3.796186436458346579e+46 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -1.5, -0.99, &r), 0.14733409946001025146 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -8.5, 0.99, &r), -1.1301780432998743440e+65 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -8.5, -0.99, &r), -8.856462606575344483 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -21.5, 0.99, &r), 2.0712920991876073253e+95 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -21.5, -0.99, &r), -74.30517015382249216 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -100.5, 0.99, &y, &r),  -3.186778061428268980e+262 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -100.5, -0.99, &y, &r),  2.4454358338375677520 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (25, 25, 1, -0.5, &r), -2.9995530823639545027e-06 , TEST_TOL0, GSL_SUCCESS);


  /* 2F1 conj */

  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (1, 1, 1, 0.5, &r), 3.352857095662929028 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 1, 0.5, &r), 1.7078067538891293983e+09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 5, 0.5, &r), 285767.15696901140627 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 1, -0.5, &r), 0.007248196261471276276 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 5, -0.5, &r), 0.00023301916814505902809 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (25, 25, 1, -0.5, &r), 5.169694409566320627e-06 , TEST_TOL0, GSL_SUCCESS);

  /* FIXME: the "true" values here may not be so good */
  /*
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (0.01, 1.0, -0.02, &r), 0.999803886708565    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (0.1,  0.5, -0.02, &r), 0.999015947934831    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (1,   1, -0.02, &r), 0.980755496569062     , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (8,   8, -0.02, &r), 0.3299059284994299    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (50, 50, -0.02, &r), 2.688995263773233e-13 , TEST_TOL0, GSL_SUCCESS);
  */


  /* 2F1 renorm */

  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (1, 1, 1, 0.5, &r), 2.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, 8, 1, 0.5, &r), 12451584.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, 1, 0.5, &r), 0.13671875 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, 1, -0.5, &r), 4945.13671875 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, -5.5, 0.5, &r), -83081.19167659493609 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, -5.5, -0.5, &r), 2.2510895952730178518e+06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, 8, 5, 0.5, &r), 175.2380952380952381 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (9, 9, -1.5, 0.99, &r), 1.6063266334913066551e+46 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (9, 9, -1.5, -0.99, &r), 0.06234327316254516616 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (5, 5, -1, 0.5, &r), 4949760.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (5, 5, -10, 0.5, &r), 139408493229637632000.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (5, 5, -100, 0.5, &r), 3.0200107544594411315e+206 , TEST_TOL0, GSL_SUCCESS);
  
  /* 2F1 conj renorm */

  TEST_SF(s, gsl_sf_hyperg_2F1_conj_renorm_e, (9, 9, -1.5, 0.99, &r), 5.912269095984229412e+49 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_renorm_e, (9, 9, -1.5, -0.99, &r), 0.10834020229476124874 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5, 5, -1, 0.5, &r), 1.4885106335357933625e+08 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5, 5, -10, 0.5, &r), 7.968479361426355095e+21 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5, 5, -100, 0.5, &r), 3.1113180227052313057e+208 , TEST_TOL0, GSL_SUCCESS);


  return status;
}

int check_jac(void)
{
  double u, m;
  double sn, cn, dn;
  int stat_ej;
  int status = 0;
  int s;
  
  u = 0.5;
  m = 0.5;
  s = 0;
  stat_ej = gsl_sf_elljac_impl(u, m, &sn, &cn, &dn);
  s += ( frac_diff( sn, 0.4707504736556572833 ) > TOL );
  s += ( frac_diff( cn, 0.8822663948904402865 ) > TOL );
  s += ( frac_diff( dn, 0.9429724257773856873 ) > TOL );
  gsl_test(s, "  gsl_sf_elljac(0.5|0.5)");
  status += s;

  u = 2.0;
  m = 0.999999;
  s = 0;
  stat_ej = gsl_sf_elljac_impl(u, m, &sn, &cn, &dn);
  s += ( frac_diff( sn, 0.96402778575700186570 ) > TOL );
  s += ( frac_diff( cn, 0.26580148285600686381 ) > TOL );
  s += ( frac_diff( dn, 0.26580323105264131136 ) > TOL );
  gsl_test(s, "  gsl_sf_elljac(2.0|0.999999)");
  status += s;

  return status;
}


int check_laguerre(void)
{
  int status = 0;
  int s;

  TEST_SF(s,  gsl_sf_laguerre_1_e, (0.5, -1.0, &r), 2.5 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_1_e, (0.5,  1.0, &r), 0.5 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_1_e, (1.0,  1.0, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_laguerre_2_e, (0.5, -1.0, &r), 4.875 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_2_e, (0.5,  1.0, &r), -0.125 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_2_e, (1.0,  1.0, &r), 0.5   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_laguerre_3_e, (0.5, -1.0, &r), 8.479166666666666667   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_3_e, (0.5,  1.0, &r), -0.6041666666666666667  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_3_e, (1.0,  1.0, &r), -0.16666666666666666667 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_3_e, (2.0,  1.0, &r), 2.3333333333333333333  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_laguerre_4_e, (0.5, -1.0, &r), 13.752604166666666667 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_4_e, (0.5,  1.0, &r), -0.8723958333333333333 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_4_e, (1.0,  1.0, &r), -0.7916666666666666667 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_4_e, (2.0,  1.0, &r), 1.5416666666666666667 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_4_e, (2.0,  0.5, &r), 6.752604166666666667  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_laguerre_5_e, (0.5, -1.0, &r), 21.24921875 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_5_e, (0.5,  1.0, &r), -0.9393229166666666667  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_5_e, (1.0,  1.0, &r), -1.2583333333333333333  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_5_e, (2.0,  1.0, &r), 0.28333333333333333333 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_5_e, (2.0,  0.5, &r), 7.45546875             , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_laguerre_n_e, (1, 0.5, 1.0, &r), 0.5 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_e, (2, 1.0, 1.0, &r), 0.5 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_e, (3, 2.0, 1.0, &r), 2.3333333333333333333 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_e, (4, 2.0, 0.5, &r), 6.752604166666666667  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_e, (90, 2.0,  0.5, &r), -48.79047157201507897  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_e, (90, 2.0, -100.0, &r), 2.5295879275042410902e+63 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_e, (90, 2.0,  100.0, &r), -2.0929042259546928670e+20 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_e, (100, 2.0, -0.5, &r), 2.2521795545919391405e+07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_e, (100, 2.0,  0.5, &r), -28.764832945909097418 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_e, (1000, 2.0, -0.5, &r), 2.4399915170947549589e+21 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_e, (1000, 2.0,  0.5, &r), -306.77440254315317525 , TEST_TOL0, GSL_SUCCESS); /**/
  TEST_SF(s,  gsl_sf_laguerre_n_e, (100000, 2.0, 1.0, &r), 5107.73491348319 , TEST_TOL0, GSL_SUCCESS);

  status += s;

  return status;
}


int check_legendre(void)
{
  double L[256];
  double y;
  int status = 0;
  int s;

  TEST_SF(s,  gsl_sf_legendre_P1_e, (-0.5, &r), -0.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P1_e, ( 0.5, &r), 0.5, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_legendre_P2_e, (0.0, &r), -0.5   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P2_e, (0.5, &r), -0.125 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P2_e, (1.0, &r), 1.0   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P2_e, (100.0, &r), 14999.5   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_legendre_P3_e, ( -0.5, &r), 0.4375 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P3_e, (  0.5, &r), -0.4375 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P3_e, (  1.0, &r), 1.0         , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P3_e, (100.0, &r), 2.49985e+06 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_legendre_P4_e, (  0.0, &r), 0.375     , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P4_e, (  0.5, &r), -0.2890625 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P4_e, (  1.0, &r), 1.0       , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P4_e, (100.0, &r), 4.37462500375e+08 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_legendre_P5_e, (  -0.5, &r), -0.08984375 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P5_e, (   0.5, &r), 0.08984375 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P5_e, (   1.0, &r), 1.0        , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P5_e, ( 100.0, &r), 7.87412501875e+10 , TEST_TOL0, GSL_SUCCESS);

  s = 0;
  s += (frac_diff( gsl_sf_legendre_Pl(1, -0.5),    -0.5 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(1,  1.0e-8),  1.0e-08 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(1,  0.5),     0.5 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(1,  1.0),     1.0 ) > TOL );  
  gsl_test(s, "  gsl_sf_legendre_Pl(1)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_Pl(10, -0.5),    -0.1882286071777345     ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(10,  1.0e-8), -0.24609374999999864648 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(10,  0.5),    -0.18822860717773437500 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(10,  1.0),     1.0 ) > TOL );  
  gsl_test(s, "  gsl_sf_legendre_Pl(10)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_Pl(99, -0.5),     0.08300778172138770477   ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(99,  1.0e-8), -7.958923738716563193e-08 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(99,  0.5),    -0.08300778172138770477   ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(99,  0.999),  -0.3317727359254778874    ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(99,  1.0),     1.0 ) > TOL );  
  gsl_test(s, "  gsl_sf_legendre_Pl(99)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_Pl(1000, -0.5),    -0.019168251091650277878 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(1000,  1.0e-8),  0.02522501817709828     ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(1000,  0.5),    -0.019168251091650277878 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(1000,  1.0),     1.0 ) > TOL );  
  gsl_test(s, "  gsl_sf_legendre_Pl(1000)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_Pl(4000, -0.5),    -0.009585404456573080972 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(4000,  0.5),    -0.009585404456573080972 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Pl(4000,  1.0),    1.0 ) > TOL );  
  gsl_test(s, "  gsl_sf_legendre_Pl(4000)");
  status += s;

  s = 0;
  gsl_sf_legendre_Pl_array_impl(100, 0.5, L);
  s += ( frac_diff(L[0],    1.0 ) > TOL );
  s += ( frac_diff(L[10],  -0.18822860717773437500 ) > TOL );
  s += ( frac_diff(L[100], -0.06051802596186118687 ) > TOL );
  gsl_test(s, "  gsl_sf_legendre_Pl_array_impl(100)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_Plm(10, 0, -0.5),	-0.18822860717773437500 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Plm(10, 0,  1.0e-8), -0.24609374999999864648 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Plm(10, 0,  0.5),	-0.18822860717773437500 ) > TOL );
  gsl_test(s, "  gsl_sf_legendre_Plm(10, 0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_Plm(10, 1, -0.5),    -2.0066877394361256516     ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Plm(10, 1,  1.0e-8), -2.7070312499999951725e-07 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Plm(10, 1,  0.5),     2.0066877394361256516     ) > TOL );
  gsl_test(s, "  gsl_sf_legendre_Plm(10, 1)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_Plm(10, 5, -0.5),    -30086.169706116174977 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Plm(10, 5,  1.0e-8), -0.0025337812499999964949 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Plm(10, 5,  0.5),	 30086.169706116174977 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Plm(10, 5,  0.999),	-0.5036411489013270406 ) > TOL );  
  gsl_test(s, "  gsl_sf_legendre_Plm(10, 5)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_Plm(100, 5, -0.5),    -6.617107444248382171e+08 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Plm(100, 5,  1.0e-8),  817.8987598063712851     ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Plm(100, 5,  0.5),     6.617107444248382171e+08 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_Plm(100, 5,  0.999),	-1.9831610803806212189e+09 ) > TOL );  
  gsl_test(s, "  gsl_sf_legendre_Plm(100, 5)");
  status += s;

  s = 0;
  gsl_sf_legendre_Plm_array_impl(100, 5, 0.5, L);
  s += ( frac_diff(L[0],  -460.3466286991656682 ) > TOL );
  s += ( frac_diff(L[10],  38852.51334152290535 ) > TOL );
  s += ( frac_diff(L[95],  6.617107444248382171e+08 ) > TOL );
  gsl_test(s, "  gsl_sf_legendre_Plm_array_impl(100, 5, 0.5)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_sphPlm(10, 0, -0.5),    -0.24332702369300133776 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_sphPlm(10, 0,  0.5),    -0.24332702369300133776 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_sphPlm(10, 0,  0.999),   1.2225754122797385990  ) > 1.0e-12 );  
  gsl_test(s, "  gsl_sf_legendre_sphPlm(10, 5)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_sphPlm(10, 5, -0.5),    -0.3725739049803293972     ) > TOL );
  s += (frac_diff( gsl_sf_legendre_sphPlm(10, 5,  1.0e-8), -3.1377233589376792243e-08 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_sphPlm(10, 5,  0.5),	    0.3725739049803293972     ) > TOL );
  s += (frac_diff( gsl_sf_legendre_sphPlm(10, 5,  0.999),  -6.236870674727370094e-06  ) > 1.0e-12 );  
  gsl_test(s, "  gsl_sf_legendre_sphPlm(10, 5)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_sphPlm(10, 10, -0.5),    0.12876871185785724117 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_sphPlm(10, 10,  0.5),    0.12876871185785724117 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_sphPlm(10, 10,  0.999),  1.7320802307583118647e-14 ) > 1.0e-12 );  
  gsl_test(s, "  gsl_sf_legendre_sphPlm(10, 10)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_legendre_sphPlm(200, 1, -0.5),    0.3302975570099492931 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_sphPlm(200, 1,  0.5),   -0.3302975570099492931 ) > TOL );
  s += (frac_diff( gsl_sf_legendre_sphPlm(200, 1,  0.999), -1.4069792055546256912 ) > 1.0e-12 );  
  gsl_test(s, "  gsl_sf_legendre_sphPlm(200, 1)");
  status += s;

  s = 0;
  gsl_sf_legendre_sphPlm_array_impl(100, 5, 0.5, L);
  s += ( frac_diff(L[0],   -0.22609703187800460722 ) > TOL );
  s += ( frac_diff(L[10],   0.07452710323813558940 ) > TOL );
  s += ( frac_diff(L[95],   0.25865355990880161717 ) > TOL );
  gsl_test(s, "  gsl_sf_legendre_sphPlm_array_impl(100, 5, 0.5)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_half(0.0, -0.5),   0.8573827581049917129  ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_half(0.0,  0.5),   0.8573827581049917129  ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_half(0.0,  2.0),   0.6062611623284649811  ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_half(0.0,  100.0), 0.07979045091636735635 ) > TOL );
  gsl_test(s, "  gsl_sf_conicalP_half(0.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_half(10.0, -0.5),  5.345484922591867188e+08  ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_half(10.0,  0.5),  15137.910380385258370     ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_half(10.0,  2.0),     0.4992680691891618544  ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_half(10.0,  100.0),  -0.07272008163718195685 ) > TOL );
  gsl_test(s, "  gsl_sf_conicalP_half(10.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_half(200.0, -1.0e-3), 1.3347639529084185010e+136 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_half(200.0,  1.0e-8), 1.0928098010940058507e+136 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_half(200.0,  0.5),    3.895546021611205442e+90   ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_half(200.0,  10.0),  -0.04308567180833581268   ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_half(200.0,  100.0), -0.04694669186576399194   ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_half(200.0,  1000.0),  0.023698140704121273277 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_half(200.0,  1.0e+8), -0.00006790983312124277891 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_conicalP_half(200.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_half(1.0e+8,  1.1),    1.1599311133054742944  ) > 1.0e-08 );
  s += (frac_diff( gsl_sf_conicalP_half(1.0e+8,  100.0),  0.07971967557381557875 ) > 1.0e-08 );
  gsl_test(s, "  gsl_sf_conicalP_half(1.0e+8)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_mhalf(0.0, -0.5),   1.7956982494514644808 ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_mhalf(0.0,  0.5),   0.8978491247257322404 ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_mhalf(0.0,  2.0),   0.7984204253272901551 ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_mhalf(0.0,  100.0), 0.4227531369388072584 ) > TOL );
  gsl_test(s, "  gsl_sf_conicalP_mhalf(0.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_mhalf(10.0, -0.5),  5.345484922591867181e+07 ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_mhalf(10.0,  0.5),  1513.7910356104985334    ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_mhalf(10.0,  2.0),  0.03439243987215615642   ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_mhalf(10.0,  100.0), 0.003283756665952609624 ) > TOL );
  gsl_test(s, "  gsl_sf_conicalP_mhalf(10.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_mhalf(200.0, -0.5),    1.7699538115312304280e+179 ) > 1.0e-13 );
  s += (frac_diff( gsl_sf_conicalP_mhalf(200.0,  1.0e-8), 5.464049005470029253e+133  ) > 1.0e-13 );
  s += (frac_diff( gsl_sf_conicalP_mhalf(200.0,  0.5),    1.9477730108056027211e+88  ) > 1.0e-13 );
  s += (frac_diff( gsl_sf_conicalP_mhalf(200.0,  10.0),    0.0012462575917716355362  ) > 1.0e-13 );
  s += (frac_diff( gsl_sf_conicalP_mhalf(200.0,  100.0),  -0.0003225881344802625149  ) > 1.0e-13 );
  s += (frac_diff( gsl_sf_conicalP_mhalf(200.0,  1000.0), -0.00004330652890886567623 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_mhalf(200.0,  1.0e+8),  2.0943091278037078483e-07 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_conicalP_mhalf(200.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_mhalf(1.0e+8,  1.1),    2.092320445620989618e-09 ) > 1.0e-06 );
  s += (frac_diff( gsl_sf_conicalP_mhalf(1.0e+8,  100.0), -3.359967833599016923e-11 ) > 1.0e-06 );  
  gsl_test(s, "  gsl_sf_conicalP_mhalf(1.0e+8)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_0(0.0, -0.5),   1.3728805006183501647  ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_0(0.0,  0.5),   1.0731820071493643751  ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_0(0.0,  2.0),   0.9012862993604472987  ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_0(0.0,  100.0), 0.30091748588199264556 ) > TOL );
  gsl_test(s, "  gsl_sf_conicalP_0(0.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_0(10.0, -0.5),  1.6795592815421804669e+08 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_0(10.0,  0.5),  4826.034132009618240   ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_0(10.0,  2.0),  0.18798468917758716146 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_0(10.0,  100.0), -0.008622130749987962529 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_conicalP_0(10.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_0(200.0, -0.5),    2.502194818646823e+180 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_conicalP_0(200.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_0(1000.0,  100.0),   0.0017908817653497715844  ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_0(1000.0,  1000.0), -0.0006566893804926284301  ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_0(1000.0,  1.0e+8),  2.3167213561756390068e-06 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_conicalP_0(1000.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_1(0.0, -0.5),     0.4939371126656998499  ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_1(0.0,  0.5),     0.14933621085538265636 ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_1(0.0,  2.0),    -0.13666874968871549533 ) > TOL );
  s += (frac_diff( gsl_sf_conicalP_1(0.0,  100.0),  -0.10544528203156629098 ) > TOL );
  gsl_test(s, "  gsl_sf_conicalP_1(0.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_1(10.0, -0.5),    1.7253802958788312520e+09 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_1(10.0,  0.5),    46781.02294059967988 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_1(10.0,  2.0),    0.26613342643657444400 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_1(10.0,  100.0), -0.23281959695501029796 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_conicalP_1(10.0)");
  status += s;

  s = 0;
  /* FIXME: Mathematica gets some brain-damaged numbers for
   * these x < 0 points. I have checked what I am doing in detail,
   * and it must be right because you can do it by summing
   * manifestly positive definite quantities.
   */
  s += (frac_diff( gsl_sf_conicalP_1(200.0, -0.999),  2.71635193070012709e+270 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_1(200.0, -0.9),    4.29524931765857131e+234   ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_1(200.0, -0.5),    5.01159205956053439e+182 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_1(200.0,  0.999),  195733.0396081538 ) > 1.0e-13 );
  s += (frac_diff( gsl_sf_conicalP_1(200.0,  10.0),  -2.9272610662414349553 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_conicalP_1(200.0)");
  status += s;

  s = 0;
  s += (frac_diff( gsl_sf_conicalP_1(1000.0,  100.0),  -1.7783258105862399857 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_1(1000.0,  1000.0),  0.4535161075156427179 ) > 1.0e-12 );
  s += (frac_diff( gsl_sf_conicalP_1(1000.0,  1.0e+8),  0.0009983414549874888478 ) > 1.0e-10 );
  gsl_test(s, "  gsl_sf_conicalP_1(1000.0)");
  status += s;


  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (2, 1.0, -0.5, &r),  1.6406279287008789526 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (10, 1.0, -0.5, &r),  0.000029315266725049129448 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (20, 1.0, -0.5, &r),  7.335769429462034431e-15 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (30, 1.0, -0.5, &r),  1.3235612394267378871e-26 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (10, 1.0, 0.5, &r),  2.7016087199857873954e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (20, 1.0, 0.5, &r),  1.1782569701435933399e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (30, 1.0, 0.5, &r),  3.636240588303797919e-41 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (10, 1.0, 2.0, &r),  2.4934929626284934483e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (20, 1.0, 2.0, &r),  1.1284762488012616191e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (30, 100.0, 100.0, &r),  -1.6757772087159526048e-64 , TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (2, 1.0, -0.5, &r),  2.2048510472375258708 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (10, 1.0, -0.5, &r),  0.00007335034531618655690 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (20, 1.0, -0.5, &r),  2.5419860619212164696e-14 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (30, 1.0, -0.5, &r),  5.579714972260536827e-26 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (10, 1.0, 0.5, &r),  1.1674078819646475282e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (20, 1.0, 0.5, &r),  7.066408031229072207e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (30, 1.0, 0.5, &r),  2.6541973286862588488e-40 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (10, 1.0, 2.0, &r),  1.0736109751890863051e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (20, 1.0, 2.0, &r),  6.760965304863386741e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (30, 100.0, 100.0, &r),  -4.268753482520651007e-63 , TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0e-06, 1.0e-06, &r), 0.9999999999998333333     , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0, 0.0, &r), 1.0                       , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0, 1.0, &r), 0.7160229153604338713     , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0, 100.0, &r), -3.767437313149604566e-44  , TEST_TOL0, GSL_SUCCESS);  
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0, 500.0, &r), -6.665351935878582205e-218 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (100.0, 1.0, &r), -0.004308757035378200029   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (100.0, 10.0, &r), 7.508054627912986427e-07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1000.0, 1.0, &r), 0.0007036067909088818319  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0e+08, 1.0, &r), 7.927485371429105968e-09  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0e+08, 100.0, &r), -3.627118904186918957e-52  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0e-06, 1.0e-06, &r), 3.333333333334222222e-07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0, 1.0e-10, &r), 4.714045207910316829e-11  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0, 1.0, &r), 0.3397013994799344639     , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0, 100.0, &r), -7.200624449531811272e-44  , TEST_TOL0, GSL_SUCCESS);  
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0, 500.0, &r), 4.192260336821728677e-218 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (100.0, 0.01, &r), 0.30117664944267412324    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (100.0, 1.0, &r), -0.007393833425336299309   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (100.0, 10.0, &r), -5.031062029821254982e-07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1000.0, 0.001, &r), 0.30116875865090396421    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1000.0, 1.0, &r), -0.0004776144516074971885  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0e+08, 1.0e-08, &r), 0.30116867893975679722    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0e+08, 1.0, &r), 3.0921097047369081582e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0e+08, 100.0, &r), -6.496142701296286936e-52  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0e-06, 1.0e-06, &r),  1.1544011544013627977e-32 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 1.0e-10, &r),  2.0224912016958766992e-52 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 1.0, &r),  0.011498635037491577728 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 5.0, &r),  0.0020696945662545205776 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 7.0, &r),  -0.0017555303787488993676 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 10.0, &r),  0.00008999979724504887101 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 100.0, &r),  -4.185397793298567945e-44 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 500.0, &r),  1.4235113901091961263e-217 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 100.0, 0.001, &r),  9.642762597222417946e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 100.0, 0.002, &r),  3.0821201254308036109e-08 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 100.0, 0.01, &r),  0.00009281069019005840532 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 100.0, 1.0, &r),  -0.008043100696178624653 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 100.0, 10.0, &r),  -3.927678432813974207e-07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1000.0, 0.001, &r),  0.00009256365284253254503 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1000.0, 0.01, &r),  -0.05553733815473079983 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0e+08, 1.0e-08, &y, &r),   0.00009256115861125841299 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0e+08, 100.0, &y, &r),    -6.496143209092860765e-52  , TEST_TOL0, GSL_SUCCESS);

  s = 0;
  gsl_sf_legendre_H3d_array_impl(100, 1.0, 3.0, L);
  s += ( frac_diff(L[  0], gsl_sf_legendre_H3d(  0, 1.0, 3.0)) > 1.0e-12 );
  s += ( frac_diff(L[  1], gsl_sf_legendre_H3d(  1, 1.0, 3.0)) > 1.0e-12 );
  s += ( frac_diff(L[ 10], gsl_sf_legendre_H3d( 10, 1.0, 3.0)) > 1.0e-12 );
  s += ( frac_diff(L[100], gsl_sf_legendre_H3d(100, 1.0, 3.0)) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_legendre_H3d_array_impl(100, 1.0, 3.0)");
  status += s;


  TEST_SF(s, gsl_sf_legendre_Q0_e, (-0.5, &r), -0.5493061443340548457 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, ( 1.5, &r), 0.8047189562170501873 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Q1_e, (-0.5, &r), -0.7253469278329725772  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, ( 1.5, &r), 0.20707843432557528095 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Ql_e, (10, -0.5, &r), -0.29165813966586752393 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (10,  0.5, &r), 0.29165813966586752393 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (10,  1.5, &r), 0.000014714232718207477406 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Ql_e, (100, -0.5, &r), -0.09492507395207282096 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (100,  0.5, &r), 0.09492507395207282096 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (100,  1.5, &r), 1.1628163435044121988e-43 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Ql_e, (1000, -0.5, &r), -0.030105074974005303500 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (1000,  0.5, &r), 0.030105074974005303500 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (1000,  1.1, &r), 1.0757258447825356443e-194 , TEST_TOL0, GSL_SUCCESS);

  return status;
}


int check_log(void)
{
  double x, y;
  int status = 0;
  int s;

  TEST_SF(s, gsl_sf_log_e, (0.1, &r), -2.3025850929940456840  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_e, (1.1, &r), 0.09531017980432486004 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_e, (1000.0, &r), 6.907755278982137052   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log_abs_e, (-0.1, &r), -2.3025850929940456840  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_abs_e, (-1.1, &r), 0.09531017980432486004 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_abs_e, (-1000.0, &r), 6.907755278982137052   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_abs_e, (0.1, &r), -2.3025850929940456840  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_abs_e, (1.1, &r), 0.09531017980432486004 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_abs_e, (1000.0, &r), 6.907755278982137052   , TEST_TOL0, GSL_SUCCESS);

  s = 0;
  gsl_sf_complex_log_impl(1.0, 1.0, &x, &y);
  s += ( frac_diff( x, 0.3465735902799726547 ) > TOL );
  s += ( frac_diff( y, 0.7853981633974483096 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_log(1 + I)");
  status += s;
  
  s = 0;
  gsl_sf_complex_log_impl(1.0, -1.0, &x, &y);
  s += ( frac_diff( x,  0.3465735902799726547 ) > TOL );
  s += ( frac_diff( y, -0.7853981633974483096 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_log(1 - I)");
  status += s;
  
  s = 0;
  gsl_sf_complex_log_impl(1.0, 100.0, &x, &y);
  s += ( frac_diff( x, 4.605220183488258022 ) > TOL );
  s += ( frac_diff( y, 1.560796660108231381 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_log(1 + 100 I)");
  status += s;

  s = 0;
  gsl_sf_complex_log_impl(-1000.0, -1.0, &x, &y);
  s += ( frac_diff( x,  6.907755778981887052  ) > TOL );
  s += ( frac_diff( y, -3.1405926539231263718 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_log(-1000 - I)");
  status += s;

  s = 0;
  gsl_sf_complex_log_impl(-1.0, 0.0, &x, &y);
  s += ( frac_diff( y, 3.1415926535897932385 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_log(-1)");
  status += s;

  TEST_SF(s,  gsl_sf_log_1plusx_e, (1.0e-10, &r), 9.999999999500000000e-11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_e, (1.0e-8, &r), 9.999999950000000333e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_e, (1.0e-4, &r), 0.00009999500033330833533 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_e, (0.1, &r), 0.09531017980432486004 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_e, (0.49, &r), 0.3987761199573677730 , TEST_TOL0, GSL_SUCCESS);
  
  TEST_SF(s,  gsl_sf_log_1plusx_e, (-0.49, &r), -0.6733445532637655964 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_e, (1.0, &r), M_LN2 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_e, (-0.99, &r), -4.605170185988091368 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_log_1plusx_mx_e, (1.0e-10, &r), -4.999999999666666667e-21 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_e, (1.0e-8, &r), -4.999999966666666917e-17 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_e, (1.0e-4, &r), -4.999666691664666833e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_e, (0.1, &r), -0.004689820195675139956 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_e, (0.49, &r), -0.09122388004263222704 , TEST_TOL0, GSL_SUCCESS);
  
  TEST_SF(s,  gsl_sf_log_1plusx_mx_e, (-0.49, &r), -0.18334455326376559639 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_e, (1.0, &r), M_LN2-1.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_e, (-0.99, &r), -3.615170185988091368 , TEST_TOL0, GSL_SUCCESS);

  return status;
}

int check_poly(void)
{
  int status = 0;
  int s;

  double x;
  double y;
  double c[3]  = { 1.0, 0.5, 0.3 };
  double d[11] = { 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1};

  s = 0;
  x = 0.5;
  y = gsl_sf_poly_eval(c, 3, x);
  s += ( frac_diff(y, 1 + 0.5*x + 0.3*x*x) > TOL );
  gsl_test(s, "  gsl_sf_poly_eval({1, 0.5, 0.3}, 0.5)");
  status += s;
  
  s = 0;
  x = 1.0;
  y = gsl_sf_poly_eval(d, 11, x);
  s += ( frac_diff(y, 1.0) > TOL );
  gsl_test(s, "  gsl_sf_poly_eval({1,-1, 1, -1, 1, -1, 1, -1, 1, -1, 1}, 1.0)");
  status += s;

  return status;
}

int check_pow_int(void)
{
  int status = 0;
  int s;
  
  TEST_SF(s,  gsl_sf_pow_int_e, (2.0, 3, &r), 8.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (-2.0, 3, &r), -8.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (2.0, -3, &r), 1.0/8.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (-2.0, -3, &r), -1.0/8.0 , TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s,  gsl_sf_pow_int_e, (10.0, 4, &r), 1.0e+4 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (10.0, -4, &r), 1.0e-4 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (-10.0, 4, &r), 1.0e+4 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (-10.0, -4, &r), 1.0e-4 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_pow_int_e, (10.0, 40, &r), 1.0e+40 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (8.0, -40, &r), 7.523163845262640051e-37 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (-10.0, 40, &r), 1.0e+40 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (-8.0, -40, &r), 7.523163845262640051e-37 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_pow_int_e, (10.0, 41, &r), 1.0e+41 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (8.0, -41, &r), 9.403954806578300064e-38 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (-10.0, 41, &r), -1.0e+41 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_e, (-8.0, -41, &r), -9.403954806578300064e-38 , TEST_TOL0, GSL_SUCCESS);

  return status;
}

int check_psi(void)
{
  int status = 0;
  int s;
  
  TEST_SF(s, gsl_sf_psi_int_e, (5, &r), 1.5061176684318004727 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_int_e, (100, &r), 4.600161852738087400 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_int_e, (110, &r), 4.695928024251535633 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_int_e, (5000, &r), 8.517093188082904107 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_psi_e, (5.0, &r), 1.5061176684318004727 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_e, (5000.0, &r), 8.517093188082904107 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_psi_e, (-100.5, &r), 4.615124601338064117 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_e, (-1.0e+5-0.5, &r), 11.512935464924395337 , TEST_TOL0, GSL_SUCCESS);
  
  TEST_SF(s, gsl_sf_psi_1piy_e, (0.8, &r), -0.07088340212750589223 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1piy_e, (1.0, &r), 0.09465032062247697727 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1piy_e, (5.0, &r), 1.6127848446157465854  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1piy_e, (100.0, &r), 4.605178519404762003 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1piy_e, (2000.0, &r), 7.600902480375416216 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_psi_1_int_e, (5, &r), 0.22132295573711532536 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1_int_e, (100, &r), 0.010050166663333571395 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1_int_e, (110, &r), 0.009132356622022545705 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1_int_e, (500, &r), 0.0020020013333322666697 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_psi_n_e, (3, 5.0, &r), 0.021427828192755075022 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_n_e, (3, 500.0, &r), 1.6048063999872000683e-08  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_n_e, (10, 5.0, &r), -0.08675107579196581317 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_n_e, (10, 50.0, &r), -4.101091112731268288e-12 , TEST_TOL0, GSL_SUCCESS);

  return status;
}


int check_synch(void)
{
  double y;
  int status = 0;
  int s = 0;

  TEST_SF(s, gsl_sf_synchrotron_1_e, (0.01, &y, &r),  0.444973 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_synchrotron_1_e, (1.0, &y, &r),  0.651423 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_synchrotron_1_e, (10.0, &y, &r),  0.000192238 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_synchrotron_1_e, (100.0, &y, &r),  4.69759e-43 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_synchrotron_2_e, (0.01, &y, &r),  0.23098077342226277732 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_synchrotron_2_e, (1.0, &y, &r),  0.4944750621042082670 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_synchrotron_2_e, (10.0, &y, &r),  0.00018161187569530204281 , TEST_TOL0, GSL_SUCCESS);  
  TEST_SF(s, gsl_sf_synchrotron_2_e, (256.0, &y, &r),  1.3272635474353774058e-110 , TEST_TOL0, GSL_SUCCESS);  /* exp()... not my fault */

  return status;
}


int check_transport(void)
{
  int status = 0;
  int s;

  TEST_SF(s, gsl_sf_transport_2_e, (1.0e-10, &r), 9.9999999999999999999e-11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_2_e, (1.0, &r), 0.97303256135517012845 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_2_e, (3.0, &r), 2.41105004901695346199 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_2_e, (10.0, &r), 3.28432911449795173575 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_2_e, (100.0, &r), 3.28986813369645287294 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_2_e, (1.0e+05, &r), 3.28986813369645287294 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_transport_3_e, (1.0e-10, &r), 4.999999999999999999997e-21 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_e, (1.0, &r), 0.479841006572417499939 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_e, (3.0, &r), 3.210604662942246772338 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_e, (5.0, &r), 5.614386613842273228585 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_e, (10.0, &r), 7.150322712008592975030 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_e, (30.0, &r), 7.212341416160946511930 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_e, (100.0, &r), 7.212341418957565712398 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_e, (1.0e+05, &r), 7.212341418957565712398 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_transport_4_e, (1.0e-10, &r), 3.33333333333333333333e-31 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_e, (1.0e-07, &r), 3.33333333333333166666e-22 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_e, (1.0e-04, &r), 3.33333333166666666726e-13 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_e, (0.1, &r), 0.000333166726172109903824 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_e, (1.0, &r), 0.31724404523442648241 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_e, (3.0, &r), 5.96482239737147652446 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_e, (5.0, &r), 15.3597843168821829816 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_e, (10.0, &r), 25.2736676770304417334 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_e, (30.0, &r), 25.9757575220840937469 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_e, (100.0, &r), 25.9757576090673165963 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_e, (1.0e+05, &r), 25.9757576090673165963 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_transport_5_e, (1.0e-10, &r), 2.49999999999999999999e-41 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_e, (1.0e-07, &r), 2.49999999999999861111e-29 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_e, (1.0e-04, &r), 2.49999999861111111163e-17 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_e, (0.1, &r), 0.000024986116317791487410 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_e, (1.0, &r), 0.236615879239094789259153 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_e, (3.0, &r), 12.77055769104415951115760 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_e, (5.0, &r), 50.26309221817518778543615 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_e, (10.0, &r), 116.3807454024207107698556 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_e, (30.0, &r), 124.4313279083858954839911 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_e, (100.0, &r), 124.4313306172043911597639 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_e, (1.0e+05, &r), 124.43133061720439115976   , TEST_TOL0, GSL_SUCCESS);

  return status;
}


int check_trig(void)
{
  double zr, zi;
  double yr, yi;
  double x, y;
  double theta;
  int status = 0;
  int s;

  s = 0;
  y = gsl_sf_sin_pi_x(1000.5);
  s += ( frac_diff(y, 1.0 ) > TOL );
  gsl_test(s, "  gsl_sf_sin_pi_x(1000.5)");
  status += s;

  s = 0;
  y = gsl_sf_sin_pi_x(10000.0 + 1.0/65536.0);
  s += ( frac_diff(y, 0.00004793689960306688455 ) > TOL );
  gsl_test(s, "  gsl_sf_sin_pi_x(10000.0 + 1.0/65536.0)");
  status += s;

  s = 0;
  y = gsl_sf_sin_pi_x(1099511627776.0 + 1 + 0.125);
  s += ( frac_diff(y, -0.3826834323650897717 ) > TOL );
  gsl_test(s, "  gsl_sf_sin_pi_x(2^40 + 1 + 0.125)");
  status += s;


  yr = 1.0;
  yi = 5.0;
  gsl_sf_complex_sin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr, 62.44551846769653403 ) > TOL );
  s += ( frac_diff( zi, 40.09216577799840254 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_sin_impl(1 + 5 I)");
  status += s;
  
  yr = 1.0;
  yi = 5.0;
  gsl_sf_complex_cos_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr,  40.09580630629882573 ) > TOL );
  s += ( frac_diff( zi, -62.43984868079963017 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_cos_impl(1 + 5 I)");
  status += s;

  yr =   1.0;
  yi = 100.0;
  gsl_sf_complex_logsin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr, 99.3068528194400546900 ) > TOL );
  s += ( frac_diff( zi,  0.5707963267948966192 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_logsin_impl(1 + 100 I)");
  status += s;
  
  yr =    1.0;
  yi = -100.0;
  gsl_sf_complex_logsin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr,  99.3068528194400546900 ) > TOL );
  s += ( frac_diff( zi,  -0.5707963267948966192 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_logsin_impl(1 - 100 I)");
  status += s;

  yr = 5.0;
  yi = 5.0;
  gsl_sf_complex_logsin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr, 4.3068909128079757420 ) > TOL );
  s += ( frac_diff( zi, 2.8540063315538773952 ) > TOL );
  gsl_test(s, "  gsl_sf_complex_logsin_impl(5 + 5 I)");
  status += s;

  TEST_SF(s,  gsl_sf_lnsinh_e, (0.1, &r), -2.3009189815304652235  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnsinh_e, (1.0, &r), 0.16143936157119563361 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnsinh_e, (5.0, &r), 4.306807418479684201   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnsinh_e, (100.0, &r), 99.30685281944005469  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_lncosh_e, (0.125, &r), 0.007792239318898252791 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lncosh_e, (1.0, &r), 0.4337808304830271870   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lncosh_e, (5.0, &r), 4.306898218339271555    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lncosh_e, (100.0, &r), 99.30685281944005469    , TEST_TOL0, GSL_SUCCESS);

  gsl_sf_polar_to_rect_impl(10.0, M_PI/6.0, &x, &y);
  s = 0;
  s += ( frac_diff( x, 10.0 * sqrt(3) / 2.0 ) > TOL );
  s += ( frac_diff( y, 10.0 * 0.5           ) > TOL );
  gsl_test(s, "  gsl_sf_polar_to_rect_impl(10, Pi/6)");
  status += s;

  gsl_sf_polar_to_rect_impl(10.0, -2.0/3.0*M_PI, &x, &y);
  s = 0;
  s += ( frac_diff( x, 10.0 * (-0.5)           ) > TOL );
  s += ( frac_diff( y, 10.0 * (-sqrt(3) / 2.0) ) > TOL );
  gsl_test(s, "  gsl_sf_polar_to_rect_impl(10, -2/3 Pi)");
  status += s;

  theta = 5.0*M_PI + M_PI/2.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, 3.0/2.0*M_PI ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_pos_impl: theta =  11/2 Pi");
  status += s;

  theta = -5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_pos_impl: theta = -11/2 Pi");
  status += s;

  theta = 50000.0 + 1.0/65536.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, 4.6945260308194656055 ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_pos_impl: theta = 50000.0 + 1.0/65536.0");
  status += s;

  theta = 5000000.0 + 1.0/65536.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, 4.49537973053997376 ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_pos_impl: theta = 5000000.0 + 1.0/65536.0");
  status += s;

  theta = 140737488355328.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, 3.20652300406795792638 ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_pos_impl: theta = 2^47");
  status += s;

  theta = 5.0*M_PI + M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, -M_PI/2.0 ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta =  11/2 Pi");
  status += s;

  theta = -5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta = -11/2 Pi");
  status += s;

  theta =  5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta = -9/2 Pi");
  status += s;

  theta =  3.0/2.0*M_PI;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, -M_PI/2.0 ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta =  3/2 Pi");
  status += s;

  theta = -3.0/2.0*M_PI;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta = -3/2 Pi");
  status += s;

  theta = 50000.0 + 1.0/65536.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, -1.5886592763601208714 ) > TOL );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta = 50000.0 + 1.0/65536.0");
  status += s;

  return status;
}


int check_zeta(void)
{
  int status = 0;
  int s;

  TEST_SF(s, gsl_sf_zeta_int_e, (-61, &r), -3.30660898765775767257e+34, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_int_e, (-5, &r), -0.003968253968253968253968, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_int_e, (5, &r), 1.0369277551433699263313655, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_int_e, (31, &r), 1.0000000004656629065033784, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_zeta_e, (-151, &r), 8.195215221831378294e+143  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_e, (-51, &r), 9.68995788746359406565e+24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_e, (-5, &r), -0.003968253968253968253968 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_e, (-0.5, &r), -0.207886224977354566017307 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_e, (0.5, &r), -1.460354508809586812889499 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_e, (1.0-1.0/1024.0, &r), -1023.4228554489429787      , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_e, (1.0+1.0/1048576, &r), 1.0485765772157343441e+06  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_e, (5, &r), 1.036927755143369926331365 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_e, (25.5, &r), 1.000000021074106110269959 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hzeta_e, (2,  1.0, &r), 1.6449340668482264365   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_e, (2, 10.0, &r), 0.1051663356816857461   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_e, (5,  1.0, &r), 1.0369277551433699263   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_e, (5, 10.0, &r), 0.000030413798676470276 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_e, (9,  0.1, &r), 1.0000000004253980e+09  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_e, (30, 0.5, &r), 1.0737418240000053e+09  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_e, (30, 0.9, &r), 2.3589824880264765e+01  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_e, (75, 0.25, &r), 1.4272476927059599e+45  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_eta_int_e, (-91, &r), -4.945598888750002040e+94 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_e, (-51, &r), -4.363969073121683116e+40 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_e, (-5, &r), 0.25 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_e, (-1, &r), 0.25 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_e, ( 0, &r), 0.5  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_e, ( 5, &r), 0.9721197704469093059 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_e, ( 6, &r), 0.9855510912974351041 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_e, ( 20, &r), 0.9999990466115815221 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_e, ( 1000, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_eta_e, (-51.5, &r), -1.2524184036924703656e+41 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_e, (-5, &r), 0.25 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_e, (0.5, &r), 0.6048986434216303702 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_e, (0.999, &r), 0.6929872789683383574 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_e, (1.0, &r), 0.6931471805599453094 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_e, (1.0+1.0e-10, &r), 0.6931471805759321998 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_e, ( 5, &r), 0.9721197704469093059 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_e, ( 5.2, &r), 0.9755278712546684682 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_e, ( 6, &r), 0.9855510912974351041 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_e, ( 20, &r), 0.9999990466115815221 , TEST_TOL0, GSL_SUCCESS);  status += s;

  return status;
}


int main(int argc, char * argv[])
{
  gsl_test(check_airy(),       "Airy Functions");
  gsl_test(check_bessel(),     "Bessel Functions");
  gsl_test(check_cheb(),       "Chebyshev Evaluation");
  gsl_test(check_clausen(),    "Clausen Integral");
  gsl_test(check_coulomb(),    "Coulomb Wave Functions");
  gsl_test(check_coupling(),   "Coupling Coefficients");
  gsl_test(check_dawson(),     "Dawson Integral");
  gsl_test(check_debye(),      "Debye Functions");
  gsl_test(check_dilog(),      "Dilogarithm");
  gsl_test(check_elementary(), "Elementary Functions (Misc)");
  gsl_test(check_ellint(),     "Elliptic Integrals");
  gsl_test(check_jac(),        "Elliptic Functions (Jacobi)");
  gsl_test(check_erf(),        "Error Functions");
  gsl_test(check_exp(),        "Exponential Functions");
  gsl_test(check_expint(),     "Exponential/Sine/Cosine Integrals");
  gsl_test(check_fermidirac(), "Fermi-Dirac Functions");
  gsl_test(check_gamma(),      "Gamma Functions");
  gsl_test(check_gegen(),      "Gegenbauer Polynomials");
  gsl_test(check_hyperg(),     "Hypergeometric Functions");
  gsl_test(check_laguerre(),   "Laguerre Polynomials");
  gsl_test(check_legendre(),   "Legendre Functions");
  gsl_test(check_log(),        "Logarithm");
  gsl_test(check_poly(),       "Polynomial Evaluation");
  gsl_test(check_pow_int(),    "Integer Powers");
  gsl_test(check_psi(),        "Psi Functions");
  gsl_test(check_synch(),      "Synchrotron Functions");
  gsl_test(check_transport(),  "Transport Functions");
  gsl_test(check_trig(),       "Trigonometric and Related Functions");
  gsl_test(check_zeta(),       "Zeta Functions");

  return gsl_test_summary();
}
