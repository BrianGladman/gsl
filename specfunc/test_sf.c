/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_test.h>
#include <gsl_sf.h>

#include "test_sf.h"


double
test_sf_frac_diff(double x1, double x2)
{
  if(x1 == 0.0 && x2 == 0.0) {
    return 0.0;
  }
  else if(x1 <= DBL_MAX && x2 <= DBL_MAX && (x1 + x2 != 0.0)) {
    return fabs((x1-x2)/(x1+x2));
  }
  else {
    return 1.0;
  }
}


/* Check a result against a given value and tolerance.
 */
int
test_sf_check_result(char * message_buff, gsl_sf_result r, double val, double tol)
{
  int    s = 0;
  double f = test_sf_frac_diff(val, r.val);

  if(fabs(val - r.val) > 2.0*r.err) s |= TEST_SF_INCONS;
  if(r.err < 0.0)                   s |= TEST_SF_ERRNEG;
  if(f > tol)                       s |= TEST_SF_TOLBAD;

  if(s != 0) {
    char buff[2048];
    sprintf(buff, "  expected: %20.16g\n", val);
    sprintf(buff, "  obtained: %20.16g   %20.16g  %g\n", r.val, r.err, r.err/(fabs(r.val) + r.err));
    sprintf(buff, "  fracdiff: %20.16g\n", f);
    strcat(message_buff, buff);
  }

  if(s & TEST_SF_INCONS) {
    strcat(message_buff, "  value/expected not consistent within reported error\n");
  }
  if(s & TEST_SF_ERRNEG) {
    strcat(message_buff, "  reported error negative\n");
  }
  if(s & TEST_SF_TOLBAD) {
    strcat(message_buff, "  value not within tolerance of expected value\n");
  }

  return s;
}


/* Check a return value.
 */
int
test_sf_check_return(char * message_buff, int val_return, int expected_return)
{
  if(val_return != expected_return) {
    char buff[256];
    sprintf(buff, "  unexpected return code: %d\n", val_return);
    strcat(message_buff, buff);
    return TEST_SF_RETBAD;
  }
  else {
    return 0;
  }
}



int test_cheb(void)
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
  s += ( f > 200.0 * TEST_TOL0 );
  gsl_test(s, "  gsl_sf_cheb_eval_impl()");
  status += s;
  
  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    gsl_sf_cheb_eval_n_impl(cs, 25, x, &r);
    f += fabs(r.val - sin(x));
  }
  s = 0;
  s += ( f > 200.0 * TEST_TOL0 );
  gsl_test(s, "  gsl_sf_cheb_eval_n_impl()");
  status += s;

  gsl_sf_cheb_calc_impl(cs, sin);
  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    gsl_sf_cheb_eval_impl(cs, x, &r);
    f += fabs(r.val - sin(x));
  }
  s = 0;
  s += ( f > 200.0 * TEST_TOL0 );
  gsl_test(s, "  gsl_sf_cheb_calc_impl()");
  status += s;

  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    gsl_sf_cheb_eval_deriv_impl(cs, x, &r);
    f += fabs(r.val - cos(x));
  }
  s = 0;
  s += ( f > 200.0 * 10.0 * TEST_TOL0 );
  gsl_test(s, "  gsl_sf_cheb_eval_deriv_impl()");
  status += s;

  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    gsl_sf_cheb_eval_integ_impl(cs, x, &r);
    f += fabs(r.val + (1.0 + cos(x)));
  }
  s = 0;
  s += ( f > 200.0 * TEST_TOL0 );
  gsl_test(s, "  gsl_sf_cheb_eval_integ_impl()");
  status += s;

  gsl_sf_cheb_free(cs);

  return status;
}


int test_clausen(void)
{
  gsl_sf_result r;
  int s = 0;

  TEST_SF(s,  gsl_sf_clausen_impl, (M_PI/20.0, &r), 0.4478882448133546, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_clausen_impl, (M_PI/6.0, &r), 0.8643791310538927, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_clausen_impl, (M_PI/3.0, &r), 1.0149416064096535, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_clausen_impl, (  2.0*M_PI + M_PI/3.0, &r), 1.0149416064096535, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_clausen_impl, (100.0*M_PI + M_PI/3.0, &r), 1.0149416064096535, TEST_TOL0, GSL_SUCCESS);

  return s;
}


#define PRINT(n) printf("%22.18g  %22.18g  %22.18g  %22.18g\n", F[n], Fp[n], G[n], Gp[n])

int test_coulomb(void)
{
  gsl_sf_result r;
  int status = 0;
  int s = 0;
  
  const int kmax = 20;
  double F[kmax+1], Fp[kmax+1], G[kmax+1], Gp[kmax+1];
  double Fe, Ge;
  double lam_min;
  double lam_F;
  double eta, x;
  int k_G;

  TEST_SF(s, gsl_sf_hydrogenicR_1_impl, (3.0, 2.0, &r), 0.025759948256148471036, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_1_impl, (3.0, 10.0, &r), 9.724727052062819704e-13, TEST_TOL0, GSL_SUCCESS);  status += s;

  TEST_SF(s, gsl_sf_hydrogenicR_impl, (4, 0, 3.0, 2.0, &r), -0.03623182256981820062, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_impl, (4, 1, 3.0, 2.0, &r), -0.028065049083129581005, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_impl, (4, 2, 3.0, 2.0, &r), 0.14583027278668431009, TEST_TOL0, GSL_SUCCESS);  status += s;

  TEST_SF(s, gsl_sf_hydrogenicR_impl, (100,  0, 3.0, 2.0, &r), -0.00007938950980052281367, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_impl, (100, 10, 3.0, 2.0, &r), 7.112823375353605977e-12, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_impl, (100, 90, 3.0, 2.0, &r), 5.845231751418131548e-245, TEST_TOL0, GSL_SUCCESS);

#if 0
  lam_F = 0.0;
  k_G   = 0;
  eta = 1.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( test_sf_frac_diff(  F[0],  0.6849374120059439677 ) > 1.e-10 );
  s += ( test_sf_frac_diff( Fp[0], -0.7236423862556063963 ) > 1.e-10 );
  s += ( test_sf_frac_diff(  G[0], -0.8984143590920205487 ) > 1.e-10 );
  s += ( test_sf_frac_diff( Gp[0], -0.5108047585190350106 ) > 1.e-10 );
  gsl_test(s,"  gsl_sf_coulomb_wave_FG_impl(1.0, 5.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 10.0;
  k_G   = 2;
  eta = 1.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( test_sf_frac_diff(  F[0],  0.0006423773354915823698 ) > 1.e-10 );
  s += ( test_sf_frac_diff( Fp[0],  0.0013299570958719702545 ) > 1.e-10 );
  s += ( test_sf_frac_diff(  G[0],  33.27615734455096130     ) > 1.e-10 );
  s += ( test_sf_frac_diff( Gp[0], -45.49180102261540580     ) > 1.e-10 );
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

int test_coupling(void)
{
  gsl_sf_result r;
  int s = 0;

  TEST_SF(s, gsl_sf_coupling_3j_impl, (0, 1, 1, 0, 1, -1, &r), sqrt(1.0/2.0), TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_coupling_3j_impl, (1, 1, 2, 1, -1, 0, &r), sqrt(1.0/6.0), TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_coupling_3j_impl, (2, 4, 6, 0, 2, -2, &r), sqrt(8.0/105.0), TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_coupling_6j_impl, (2, 2, 4, 2, 2, 2, &r), 1.0/6.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_coupling_9j_impl, (4, 2, 4, 3, 3, 2, 1, 1, 2, &r), -0.040824829046386, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_coupling_9j_impl, (8, 4, 10, 7, 3, 8, 1, 1, 2, &r), 0.025458753860866, TEST_TOL0, GSL_SUCCESS);

  return s;
}

int test_dawson(void)
{
  gsl_sf_result r;
  int status = 0;
  int s;

  TEST_SF(s,  gsl_sf_dawson_impl, (1.0e-15, &r), 1.0e-15, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dawson_impl, (0.5, &r), 0.4244363835020222959, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dawson_impl, (2.0, &r), 0.30134038892379196603, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dawson_impl, (1000.0, &r), 0.0005000002500003750009, TEST_TOL0, GSL_SUCCESS);
  
  return status;
}

int test_debye(void)
{
  gsl_sf_result r;
  int s = 0;

  TEST_SF(s, gsl_sf_debye_1_impl, (0.1, &r),  0.975277750004723276, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_1_impl, (1.0, &r),  0.777504634112248239, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_1_impl, (10.0, &r), 0.164443465679946027, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_debye_2_impl, (0.1, &r),  0.967083287045302664,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_2_impl, (1.0, &r),  0.70787847562782924,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_2_impl, (10.0, &r), 0.0479714980201218708, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_debye_3_impl, (0.1, &r),  0.962999940487211048,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_3_impl, (1.0, &r),  0.674415564077814667,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_3_impl, (10.0, &r), 0.0192957656903454886, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_debye_4_impl, (0.1, &r),  0.960555486124335944,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_4_impl, (1.0, &r),  0.654874068886737049,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_debye_4_impl, (10.0, &r), 0.00967367556027115896, TEST_TOL0, GSL_SUCCESS);

  return s;
}


int test_dilog(void)
{
  gsl_sf_result r;
  gsl_sf_result r1, r2;
  int s = 0;

  /* real dilog */

  TEST_SF(s,  gsl_sf_dilog_impl, (-3.0, &r), -1.9393754207667089531, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (-0.5, &r), -0.4484142069236462024, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (-0.001, &r), -0.0009997501110486510834, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (0.1, &r), 0.1026177910993911, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (0.7, &r), 0.8893776242860387386, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (1.0, &r), 1.6449340668482260, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (1.5, &r), 2.3743952702724802007, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (2.0, &r), 2.4674011002723397, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, ( 5.0, &r), 1.7837191612666306277, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, ( 11.0, &r), 0.3218540439999117111, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (12.59, &r), 0.0010060918167266208634, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (12.595, &r), 0.00003314826006436236810, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (13.0, &r), -0.07806971248458575855, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (20.0, &r), -1.2479770861745251168, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (150.0, &r), -9.270042702348657270, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_dilog_impl, (1100.0, &r), -21.232504073931749553, TEST_TOL0, GSL_SUCCESS);


  /* complex dilog */
  /* FIXME: probably need more tests here... 
   * also need to work on accuracy for r->1; need to
   * adjust the switch-over point I suppose.
   */

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (1.00001, M_PI/2.0, &r1, &r2),
            -0.20562022409960237363, TEST_TOL0,
             0.91597344814458309320, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (0.99999, M_PI/2.0, &r1, &r2),
            -0.20561329262779687646, TEST_TOL0,
             0.91595774018131512060, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (0.991, M_PI/2.0, &r1, &r2),
            -0.20250384721077806127, TEST_TOL0,
             0.90888544355846447810, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (0.98, M_PI/2.0, &r1, &r2),
            -0.19871638377785918403, TEST_TOL0,
             0.90020045882981847610, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (0.95, M_PI/2.0, &r1, &r2),
            -0.18848636456893572091, TEST_TOL0,
             0.87633754133420277830, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (0.8, M_PI/2.0, &r1, &r2),
            -0.13980800855429037810, TEST_TOL0,
             0.75310609092419884460, TEST_TOL0,
	     GSL_SUCCESS);


  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (0.5, M_PI/2.0, &r1, &r2),
            -0.05897507442156586346, TEST_TOL0,
             0.48722235829452235710, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (0.01, M_PI/2.0, &r1, &r2),
            -0.000024999375027776215378, TEST_TOL0,
             0.009999888892888684820, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (10.0, M_PI/2.0, &r1, &r2),
            -3.0596887943287347304, TEST_TOL0,
             3.7167814930680685900, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (100.0, M_PI/2.0, &r1, &r2),
            -11.015004738293824854, TEST_TOL0,
             7.2437843013083534970, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (0.99, M_PI/8.0, &r1, &r2),
            1.0571539648820244720, TEST_TOL0,
            0.7469145254610851318, TEST_TOL0,
	    GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (0.99, M_PI/64.0, &r1, &r2),
            1.5381800285902999666, TEST_TOL0,
            0.1825271634987756651, TEST_TOL0,
	    GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_impl, (0.9, 3.0*M_PI/4.0, &r1, &r2),
            -0.6062840301356530985, TEST_TOL0,
             0.4836632833122775721, TEST_TOL0,
	    GSL_SUCCESS);

  return s;
}


int test_elementary(void)
{
  gsl_sf_result r;
  double x = 0.2*DBL_MAX;
  int s = 0;

  TEST_SF(s,  gsl_sf_multiply_impl, (-3.0,2.0, &r), -6.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_multiply_impl, (x, 1.0/x, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_multiply_impl, (x, 0.2, &r), 0.2, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_multiply_impl, (x, 2.0, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  s += ( gsl_sf_multiply_impl(DBL_MAX, 1.1, &r) != GSL_EOVRFLW);
  s += ( gsl_sf_multiply_impl(DBL_MIN, 0.9, &r) != GSL_EUNDRFLW);

  return s;
}


int test_ellint(void)
{
  gsl_sf_result r;
  gsl_mode_t mode = GSL_MODE_DEFAULT;
  int s = 0;
  
  TEST_SF(s,  gsl_sf_ellint_Kcomp_impl, ( 0.99, mode, &r), 3.3566005233611923760, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_Kcomp_impl, ( 0.50, mode, &r), 1.6857503548125960429, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_Kcomp_impl, (0.010, mode, &r), 1.5708355989121522360, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_Ecomp_impl, (0.99, mode, &r), 1.0284758090288040010, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_Ecomp_impl, (0.50, mode, &r), 1.4674622093394271555, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_Ecomp_impl, (0.01, mode, &r), 1.5707570561503852873, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_F_impl, (M_PI/3.0, 0.99, mode, &r), 1.3065333392738766762, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_F_impl, (M_PI/3.0, 0.50, mode, &r), 1.0895506700518854093, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_F_impl, (M_PI/3.0, 0.01, mode, &r), 1.0472129063770918952, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_E_impl, (M_PI/3.0, 0.99, mode, &r), 0.8704819220377943536, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_E_impl, (M_PI/3.0, 0.50, mode, &r), 1.0075555551444720293, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_E_impl, (M_PI/3.0, 0.01, mode, &r), 1.0471821963889481104, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_P_impl, (M_PI/3.0, 0.99, 0.5, mode, &r), 1.1288726598764099882, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_P_impl, (M_PI/3.0, 0.50, 0.5, mode, &r), 0.9570574331323584890, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_ellint_P_impl, (M_PI/3.0, 0.01, 0.5, mode, &r), 0.9228868127118118465, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RF_impl, (5.0e-11, 1.0e-10, 1.0, mode, &r), 12.36441982979439, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RF_impl, (1.0, 2.0, 3.0, mode, &r), 0.7269459354689082, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RD_impl, (5.0e-11, 1.0e-10, 1.0, mode, &r), 34.0932594919337362, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RD_impl, (1.0, 2.0, 3.0, mode, &r), 0.2904602810289906, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RC_impl, (1.0, 2.0, mode, &r), 0.7853981633974482, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_ellint_RJ_impl, (2.0, 3.0, 4.0, 5.0, mode, &r), 0.1429757966715675, TEST_TOL0, GSL_SUCCESS);


  return s;
}


int test_erf(void)
{
  gsl_sf_result r;
  int s = 0;

  TEST_SF(s, gsl_sf_erfc_impl, (-10.0, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_impl, (-1.0, &r), 1.8427007929497148693, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_impl, (-0.5, &r), 1.5204998778130465377, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_impl, (1.0, &r), 0.15729920705028513066, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_impl, (3.0, &r), 0.000022090496998585441373, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_impl, (7.0, &r), 4.183825607779414399e-23, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erfc_impl, (10.0, &r), 2.0884875837625447570e-45, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log_erfc_impl, (-1.0, &r), log(1.842700792949714869) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_erfc_impl, (1.0, &r), log(0.15729920705028513066) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_erfc_impl, (10.0, &r), log(2.0884875837625447570e-45) , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_erf_impl, (-10.0, &r), -1.0000000000000000000, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erf_impl, (0.5, &r), 0.5204998778130465377, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erf_impl, (1.0, &r), 0.8427007929497148693, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erf_impl, (10.0, &r), 1.0000000000000000000, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_erf_Z_impl, (1.0, &r), 0.24197072451914334980, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_erf_Q_impl, (10.0, &r), 7.619853024160526066e-24 , TEST_TOL0, GSL_SUCCESS);

  return s;
}


int test_exp(void)
{
  gsl_sf_result r;
  double x;
  int s = 0;

  TEST_SF(s,  gsl_sf_exp_impl, (-10.0, &r), exp(-10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exp_impl, ( 10.0, &r), exp( 10.0) , TEST_TOL0, GSL_SUCCESS);

  x = 0.8*GSL_LOG_DBL_MAX;
  TEST_SF(s, gsl_sf_exp_mult_impl, (-10.0,  1.0e-06, &r), 1.0e-06*exp(-10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_impl, (-10.0,  2.0, &r), 2.0*exp(-10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_impl, (-10.0, -2.0, &r), -2.0*exp(-10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_impl, ( 10.0,  1.0e-06, &r), 1.0e-06*exp( 10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_impl, ( 10.0, -2.0, &r), -2.0*exp( 10.0) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_impl, (x, 1.00001, &r), 1.00001*exp(x) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_impl, (x, 1.000001, &r), 1.000001*exp(x) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_impl, (x, 1.000000001, &r), 1.000000001*exp(x) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_impl, (x, 100.0, &r), 100.0*exp(x) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_impl, (x, 1.0e+20, &r), 1.0e+20*exp(x) , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_impl, (x, exp(-x)*exp(M_LN2), &r),  2.0, 1.0e-13, GSL_SUCCESS );

  TEST_SF(s,  gsl_sf_expm1_impl, (-10.0, &r), exp(-10.0)-1.0          , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expm1_impl, (-0.001, &r), -0.00099950016662500845, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expm1_impl, (-1.0e-8, &r), -1.0e-08 + 0.5e-16      , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expm1_impl, ( 1.0e-8, &r), 1.0e-08 + 0.5e-16      , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expm1_impl, ( 0.001, &r), 0.0010005001667083417  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expm1_impl, ( 10.0, &r), exp(10.0)-1.0           , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exprel_impl, (-10.0, &r), 0.0999954600070237515 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_impl, (-0.001, &r), 0.9995001666250084    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_impl, (-1.0e-8, &r), 1.0 - 0.5e-08       , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_impl, ( 1.0e-8, &r), 1.0 + 0.5e-08       , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_impl, ( 0.001, &r), 1.0005001667083417   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_impl, ( 10.0, &r), 2202.5465794806716517 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exprel_2_impl, (-10.0, &r), 0.18000090799859524970 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_2_impl, (-0.001, &r), 0.9996667499833361107  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_2_impl, (-1.0e-8, &r), 0.9999999966666666750  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_2_impl, ( 1.0e-8, &r), 1.0000000033333333417  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_2_impl, ( 0.001, &r), 1.0003334166833361115  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_2_impl, ( 10.0, &r), 440.3093158961343303   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exprel_n_impl, (3, -1000.0, &r), 0.00299400600000000000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3, -100.0, &r), 0.02940600000000000000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3, -10.0, &r), 0.24599972760042142509 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3, -3.0, &r), 0.5444917625849191238  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3, -0.001, &r), 0.9997500499916678570  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3, -1.0e-8, &r), 0.9999999975000000050  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3,  1.0e-8, &r), 1.0000000025000000050  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3,  0.001, &r), 1.0002500500083345240  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3,  3.0, &r), 2.5745637607083706091  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3,  3.1, &r), 2.6772417068460206247  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3,  10.0, &r), 131.79279476884029910  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (3,  100.0, &r), 1.6128702850896812690e+38 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exprel_n_impl, (50, -1000.0, &r), 0.04766231609253975959 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50, -100.0, &r), 0.3348247572345889317  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50, -10.0, &r), 0.8356287051853286482  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50, -3.0, &r), 0.9443881609152163615  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50, -1.0, &r), 0.980762245565660617   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50, -1.0e-8, &r), 1.0 -1.0e-8/51.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50,  1.0e-8, &r), 1.0 +1.0e-8/51.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50,  1.0, &r), 1.01999216583666790   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50,  3.0, &r), 1.0624205757460368307 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50,  48.0, &r), 7.499573876877194416  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50,  50.1, &r), 9.311803306230992272  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50,  100.0, &r), 8.175664432485807634e+07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (50,  500.0, &r), 4.806352370663185330e+146 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_exprel_n_impl, (500, -1000.0, &r), 0.3334815803127619256 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500, -100.0, &r), 0.8335646217536183909 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500, -10.0, &r), 0.9804297803131823066 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500, -3.0, &r), 0.9940475488850672997 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500, -1.0, &r), 0.9980079602383488808 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500, -1.0e-8, &r), 1.0 -1.0e-8/501.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500,  1.0e-8, &r), 1.0 +1.0e-8/501.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500,  1.0, &r), 1.0019999920160634252 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500,  3.0, &r), 1.0060240236632444934 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500,  48.0, &r), 1.1059355517981272174 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500,  100.0, &r), 1.2492221464878287204 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500,  500.0, &r), 28.363019877927630858 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500,  1000.0, &r), 2.4037563160335300322e+68 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_exprel_n_impl, (500,  1600.0, &r), 7.899293535320607403e+226 , TEST_TOL0, GSL_SUCCESS);

  return s;
}


int test_expint(void)
{
  gsl_sf_result r;
  int s;

  TEST_SF(s,  gsl_sf_expint_E1_impl, (-1.0, &r), -1.8951178163559367555  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_impl, (1.0e-10, &r), 22.448635265138923980  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_impl, (1.0e-05, &r), 10.935719800043695615  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_impl, (0.1, &r), 1.82292395841939066610 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_impl, (1.0, &r), 0.21938393439552027368 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_impl, (10.0, &r), 4.156968929685324277e-06  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_impl, (50.0, &r), 3.783264029550459019e-24  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E1_impl, (300.0, &r), 1.710384276804510115e-133 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_expint_E2_impl, (-1.0, &r), 0.8231640121031084799  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_impl, (1.0/4294967296.0, &r), 0.9999999947372139168  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_impl, (1.0/65536.0, &r), 0.9998243233207178845  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_impl, (0.1, &r), 0.7225450221940205066  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_impl, (1.0, &r), 0.14849550677592204792 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_impl, (10.0, &r), 3.830240465631608762e-06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_impl, (50.0, &r), 3.711783318868827367e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_E2_impl, (300.0, &r), 1.7047391998483433998e-133 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_expint_Ei_impl, (-1.0, &r), -0.21938393439552027368 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_Ei_impl, (1.0/4294967296.0, &r), -21.603494112783886397  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_Ei_impl, (1.0, &r), 1.8951178163559367555  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_Shi_impl, (-1.0, &r), -1.0572508753757285146     , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_impl, (1.0/4294967296.0, &r), 2.3283064365386962891e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_impl, (1.0/65536.0, &r), 0.00001525878906269737298 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_impl, (0.1, &r), 0.1000555722250569955 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_impl, (1.0, &r), 1.0572508753757285146 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_impl, (10.0, &r), 1246.1144901994233444 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_impl, (50.0, &r), 5.292818448565845482e+19  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Shi_impl, (300.0, &r), 3.248241254044332895e+127 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_Chi_impl, (-1.0, &r), 0.8378669409802082409 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_impl, (1.0/4294967296.0, &r), -21.603494113016717041 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_impl, (1.0/65536.0, &r), -10.513139223999384429 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_impl, (1.0/8.0, &r), -1.4983170827635760646 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_impl, (1.0, &r), 0.8378669409802082409 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_impl, (10.0, &r), 1246.1144860424544147 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_impl, (50.0, &r), 5.292818448565845482e+19  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Chi_impl, (300.0, &r), 3.248241254044332895e+127 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_expint_3_impl, (1.0e-10, &r), 1.0e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_impl, (1.0e-05, &r), 9.9999999999999975e-06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_impl, (0.1, &r), 0.09997500714119079665122 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_impl, (0.5, &r), 0.48491714311363971332427 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_impl, (1.0, &r), 0.80751118213967145285833 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_impl, (2.0, &r), 0.89295351429387631138208 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_impl, (5.0, &r), 0.89297951156924921121856 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_impl, (10.0, &r), 0.89297951156924921121856 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_expint_3_impl, (100.0, &r), 0.89297951156924921121856 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_Si_impl, (-1.0, &r), -0.9460830703671830149 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_impl, (1.0e-10, &r), 1.0e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_impl, (1.0e-05, &r), 9.999999999944444444e-06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_impl, (0.1, &r), 0.09994446110827695016   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_impl, (1.0, &r), 0.9460830703671830149    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_impl, (10.0, &r), 1.6583475942188740493    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_impl, (50.0, &r), 1.5516170724859358947 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_impl, (300.0, &r), 1.5708810882137495193 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Si_impl, (1.0e+20, &r), 1.5707963267948966192 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_Ci_impl, (1.0/4294967296.0, &r), -21.603494113016717041   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_impl, (1.0/65536.0, &r), -10.513139224115799751   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_impl, (1.0/8.0, &r), -1.5061295845296396649   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_impl, (1.0, &r), 0.3374039229009681347   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_impl, (10.0, &r), -0.04545643300445537263  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_impl, (50.0, &r), -0.005628386324116305440 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_impl, (300.0, &r), -0.003332199918592111780 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_impl, (65536.0, &r), 0.000010560248837656279453 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_impl, (4294967296.0, &r), -1.0756463261957757485e-10  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_Ci_impl, (1099511627776.0, &r), -3.689865584710764214e-13   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_atanint_impl, (1.0e-10, &r), 1.0e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_impl, (1.0e-05, &r), 9.99999999988888888889e-06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_impl, (0.1, &r), 0.09988928686033618404 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_impl, (1.0, &r), 0.91596559417721901505 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_impl, (2.0, &r), 1.57601540344632342236 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_impl, (10.0, &r), 3.71678149306806859029 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_impl, (50.0, &r), 6.16499047850274874222 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_impl, (300.0, &r), 8.96281388924518959990 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_atanint_impl, (1.0e+5, &r), 18.084471031038661920  , TEST_TOL0, GSL_SUCCESS);


  return s;
}


int test_fermidirac(void)
{
  gsl_sf_result r;
  int s = 0;

  TEST_SF(s, gsl_sf_fermi_dirac_m1_impl, (-10.0, &r), 0.00004539786870243439450 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_m1_impl, ( -1.0, &r), 0.26894142136999512075 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_m1_impl, (  1.0, &r), 0.7310585786300048793  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_m1_impl, ( 10.0, &r), 0.9999546021312975656  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_0_impl, (-10.0, &r), 0.00004539889921686464677 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_0_impl, ( -1.0, &r), 0.31326168751822283405 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_0_impl, (  1.0, &r), 1.3132616875182228340  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_0_impl, ( 10.0, &r), 10.000045398899216865  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, (-10.0, &r), 0.00004539941448447633524 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, ( -2.0, &r), 0.13101248471442377127 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, ( -1.0, &r), 0.3386479964034521798  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, ( -0.4, &r), 0.5825520806897909028  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, (  0.4, &r), 1.1423819861584355337  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, (  1.0, &r), 1.8062860704447742567  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, (  1.5, &r), 2.5581520872227806402  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, (  2.5, &r), 4.689474797599761667   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, ( 10.0, &r), 51.64488866743374196   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, ( 12.0, &r), 73.64492792264531092   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, ( 20.0, &r), 201.64493406478707282  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_1_impl, ( 50.0, &r), 1251.6449340668482264  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, (-10.0, &r), 0.00004539967212174776662 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, ( -2.0, &r), 0.13313272938565030508 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, ( -1.0, &r), 0.3525648792978077590  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, ( -0.4, &r), 0.6229402647001272120  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, (  0.4, &r), 1.2915805581060844533  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, (  1.0, &r), 2.1641656128127008622  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, (  1.5, &r), 3.247184513920792475   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, (  2.5, &r), 6.797764392735056317   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, ( 10.0, &r), 183.11605273482105278  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, ( 12.0, &r), 307.73921494638635166  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, ( 20.0, &r), 1366.2320146723590157  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, ( 50.0, &r), 20915.580036675744655  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_2_impl, (200.0, &r), 1.3336623201467029786e+06 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, (-10.0, &r), 0.00004539847236080549532 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, ( -2.0, &r), 0.12366562180120994266 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, ( -1.0, &r), 0.29402761761145122022 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, ( -0.4, &r), 0.4631755336886027800 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, (  0.4, &r), 0.7654084737661656915 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, (  1.0, &r), 1.0270571254743506890 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, (  1.5, &r), 1.2493233478527122008 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, (  2.5, &r), 1.6663128834358313625 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, ( 10.0, &r), 3.552779239536617160 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, ( 12.0, &r), 3.897268231925439359 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, ( 20.0, &r), 5.041018507535328603 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_mhalf_impl, ( 50.0, &r), 7.977530858581869960 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, (-10.0, &r), 0.00004539920105264132755 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, ( -2.0, &r), 0.12929851332007559106 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, ( -1.0, &r), 0.3277951592607115477 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, ( -0.4, &r), 0.5522452153690688947 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, (  0.4, &r), 1.0386797503389389277 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, (  1.0, &r), 1.5756407761513002308 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, (  1.5, &r), 2.1448608775831140360 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, (  2.5, &r), 3.606975377950373251  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, ( 10.0, &r), 24.084656964637653615 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, ( 12.0, &r), 31.540203287044242593 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, ( 20.0, &r), 67.49151222165892049  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_half_impl, ( 50.0, &r), 266.09281252136259343 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, (-10.0, &r), 0.00004539956540456176333 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, ( -2.0, &r), 0.13224678225177236685 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, ( -1.0, &r), 0.3466747947990574170  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, ( -0.4, &r), 0.6056120213305040910  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, (  0.4, &r), 1.2258236403963668282  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, (  1.0, &r), 2.0022581487784644573  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, (  1.5, &r), 2.9277494127932173068  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, (  2.5, &r), 5.768879312210516582   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, ( 10.0, &r), 101.00510084332600020  , TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, ( 12.0, &r), 156.51518642795728036  , TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, ( 20.0, &r), 546.5630100657601959   , TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_3half_impl, ( 50.0, &r), 5332.353566687145552   , TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (3,  -2.0, &r), 0.1342199155038680215 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (3,   0.0, &r), 0.9470328294972459176 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (3,   0.1, &r), 1.0414170610956165759 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (3,   1.0, &r), 2.3982260822489407070 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (3,   3.0, &r), 12.621635313399690724 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (3, 100.0, &r), 4.174893231066566793e+06 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (3, 500.0, &r), 2.604372285319088354e+09 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (5,  -2.0, &r), 0.13505242246823676478 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (5,   0.0, &r), 0.9855510912974351041  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (5,   0.1, &r), 1.0876519750101492782  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (5,   1.0, &r), 2.6222337848692390539  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (5,   3.0, &r), 17.008801618012113022  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (5, 100.0, &r), 1.3957522531334869874e+09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (5, 500.0, &r), 2.1705672808114817955e+13 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (7,  -2.0, &r), 0.1352641105671255851 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (7,   0.0, &r), 0.9962330018526478992 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (7,   0.1, &r), 1.1005861815180315485 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (7,   1.0, &r), 2.6918878172003129203 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (7,   3.0, &r), 19.033338976999367642 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (7,  10.0, &r), 5654.530932873610014  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (7,  50.0, &r), 1.005005069985066278e+09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (7, 500.0, &r), 9.691690268341569514e+16 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (9,  -2.0, &r), 0.1353174385330242691 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (9,   0.0, &r), 0.9990395075982715656 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (9,   0.1, &r), 1.1039997234712941212 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (9,   1.0, &r), 2.7113648898129249947 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (9,   3.0, &r), 19.768544008138602223 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (9,  10.0, &r), 10388.990167312912478 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (9,  50.0, &r), 2.85466960802601649e+10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (9, 500.0, &r), 2.69273849842695876e+20 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (10,  -2.0, &r), 0.13532635396712288092 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (10,   0.0, &r), 0.9995171434980607541 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (10,   0.1, &r), 1.1045818238852612296 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (10,   1.0, &r), 2.7147765350346120647 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (10,   3.0, &r), 19.917151938411675171 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (10,  10.0, &r), 12790.918595516495955 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (10,  50.0, &r), 1.3147703201869657654e+11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (10, 500.0, &r), 1.2241331244469204398e+22 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (11,  -2.0, &r), 0.1353308162894847149 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (11,   0.0, &r), 0.9997576851438581909 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (11,   0.1, &r), 1.1048751811565850418 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (11,   1.0, &r), 2.7165128749007313436 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (11,   3.0, &r), 19.997483022044603065 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (11,  10.0, &r), 14987.996005901818036 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (11,  50.0, &r), 5.558322924078990628e+11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (11, 500.0, &r), 5.101293089606198280e+23 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (20,  -2.0, &r), 0.13533527450327238373 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (20,   0.0, &r), 0.9999995232582155428  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (20,   0.1, &r), 1.1051703357941368203  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (20,   1.0, &r), 2.7182783069905721654  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (20,   3.0, &r), 20.085345296028242734  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (20,  10.0, &r), 21898.072920149606475  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (20,  50.0, &r), 1.236873256595717618e+16 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_fermi_dirac_int_impl, (20, 500.0, &r), 9.358938204369557277e+36 , TEST_TOL0, GSL_SUCCESS);


  return s;
}


int test_gegen(void)
{
  gsl_sf_result r;
  double ga[100];
  int s = 0;
  int sa;

  TEST_SF(s,  gsl_sf_gegenpoly_1_impl, (-0.2,   1.0, &r), -0.4 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_1_impl, ( 0.0,   1.0, &r), 2.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_1_impl, ( 1.0,   1.0, &r), 2.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_1_impl, ( 1.0,   0.5, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_1_impl, ( 5.0,   1.0, &r), 10.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_1_impl, ( 100.0, 0.5, &r), 100.0 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_gegenpoly_2_impl, (-0.2,   0.5, &r), 0.12 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_2_impl, ( 0.0,   1.0, &r), 1.00 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_2_impl, ( 1.0,   1.0, &r), 3.00 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_2_impl, ( 1.0,   0.1, &r), -0.96 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_2_impl, ( 5.0,   1.0, &r), 55.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_2_impl, ( 100.0, 0.5, &r), 4950.0 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_gegenpoly_3_impl, (-0.2,   0.5, &r), 0.112 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_3_impl, ( 0.0,   1.0, &r), -2.0/3.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_3_impl, ( 1.0,   1.0, &r), 4.000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_3_impl, ( 1.0,   0.1, &r), -0.392 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_3_impl, ( 5.0,   1.0, &r), 220.000 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_3_impl, ( 100.0, 0.5, &r), 161600.000 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_gegenpoly_n_impl, (1,       1.0, 1.0, &r), 2.000		    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_impl, (10,      1.0, 1.0, &r), 11.000		    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_impl, (10,      1.0, 0.1, &r), -0.4542309376	    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_impl, (10,      5.0, 1.0, &r), 9.23780e+4  	    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_impl, (10,    100.0, 0.5, &r), 1.5729338392690000e+13  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_impl, (1000,  100.0, 1.0, &r), 3.3353666135627322e+232 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_impl, (100,  2000.0, 1.0, &r), 5.8753432034937579e+202 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_impl, (103,   207.0, 2.0, &r), 1.4210272202235983e+145 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_gegenpoly_n_impl, (103,    -0.4, 0.3, &r), -1.64527498094522e-04    , TEST_TOL0, GSL_SUCCESS);

  sa = 0;
  gsl_sf_gegenpoly_array_impl(99, 5.0, 1.0, ga);
  sa += ( test_sf_frac_diff( ga[1],     10.0    ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff( ga[10], 9.23780e+4 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_sf_gegenpoly_array_impl");
  s += sa;

  return s;
}


int test_jac(void)
{
  double u, m;
  double sn, cn, dn;
  int stat_ej;
  int s = 0;
  int sa;

  u = 0.5;
  m = 0.5;
  sa = 0;
  stat_ej = gsl_sf_elljac_impl(u, m, &sn, &cn, &dn);
  sa += ( test_sf_frac_diff( sn, 0.4707504736556572833 ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff( cn, 0.8822663948904402865 ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff( dn, 0.9429724257773856873 ) > TEST_TOL0 );
  gsl_test(s, "  gsl_sf_elljac_impl(0.5|0.5)");
  s += sa;

  u = 2.0;
  m = 0.999999;
  sa = 0;
  stat_ej = gsl_sf_elljac_impl(u, m, &sn, &cn, &dn);
  sa += ( test_sf_frac_diff( sn, 0.96402778575700186570 ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff( cn, 0.26580148285600686381 ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff( dn, 0.26580323105264131136 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_sf_elljac_impl(2.0|0.999999)");
  s += sa;

  return s;
}


int test_laguerre(void)
{
  gsl_sf_result r;
  int s = 0;

  TEST_SF(s,  gsl_sf_laguerre_1_impl, (0.5, -1.0, &r), 2.5 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_1_impl, (0.5,  1.0, &r), 0.5 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_1_impl, (1.0,  1.0, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_laguerre_2_impl, (0.5, -1.0, &r), 4.875 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_2_impl, (0.5,  1.0, &r), -0.125 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_2_impl, (1.0,  1.0, &r), 0.5   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_laguerre_3_impl, (0.5, -1.0, &r), 8.479166666666666667   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_3_impl, (0.5,  1.0, &r), -0.6041666666666666667  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_3_impl, (1.0,  1.0, &r), -0.16666666666666666667 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_3_impl, (2.0,  1.0, &r), 2.3333333333333333333  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_laguerre_n_impl, (1, 0.5, 1.0, &r), 0.5 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (2, 1.0, 1.0, &r), 0.5 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (3, 2.0, 1.0, &r), 2.3333333333333333333 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (4, 2.0, 0.5, &r), 6.752604166666666667  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (90, 2.0,  0.5, &r), -48.79047157201507897  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (90, 2.0, -100.0, &r), 2.5295879275042410902e+63 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (90, 2.0,  100.0, &r), -2.0929042259546928670e+20 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (100, 2.0, -0.5, &r), 2.2521795545919391405e+07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (100, 2.0,  0.5, &r), -28.764832945909097418 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (1000, 2.0, -0.5, &r), 2.4399915170947549589e+21 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (1000, 2.0,  0.5, &r), -306.77440254315317525 , TEST_TOL0, GSL_SUCCESS); /**/
  TEST_SF(s,  gsl_sf_laguerre_n_impl, (100000, 2.0, 1.0, &r), 5107.73491348319 , TEST_TOL0, GSL_SUCCESS);

  return s;
}


int test_legendre(void)
{
  gsl_sf_result r;
  double L[256];
  int s = 0;
  int sa;

  TEST_SF(s,  gsl_sf_legendre_P1_impl, (-0.5, &r), -0.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P1_impl, ( 0.5, &r), 0.5, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_legendre_P2_impl, (0.0, &r), -0.5   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P2_impl, (0.5, &r), -0.125 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P2_impl, (1.0, &r), 1.0   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P2_impl, (100.0, &r), 14999.5   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_legendre_P3_impl, ( -0.5, &r), 0.4375 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P3_impl, (  0.5, &r), -0.4375 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P3_impl, (  1.0, &r), 1.0         , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P3_impl, (100.0, &r), 2.49985e+06 , TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_legendre_Pl_impl, (1, -0.5, &r), -0.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (1,  1.0e-8, &r), 1.0e-08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (1,  0.5, &r), 0.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (1,  1.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
 
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (10, -0.5, &r), -0.1882286071777345, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (10,  1.0e-8, &r), -0.24609374999999864648, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (10,  0.5, &r), -0.18822860717773437500, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (10,  1.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Pl_impl, (99, -0.5, &r), 0.08300778172138770477, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (99,  1.0e-8, &r), -7.958923738716563193e-08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (99,  0.5, &r), -0.08300778172138770477, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (99,  0.999, &r), -0.3317727359254778874, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (99,  1.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Pl_impl, (1000, -0.5, &r), -0.019168251091650277878, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (1000,  1.0e-8, &r), 0.02522501817709828, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (1000,  0.5, &r), -0.019168251091650277878, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (1000,  1.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Pl_impl, (4000, -0.5, &r), -0.009585404456573080972, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (4000,  0.5, &r), -0.009585404456573080972, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_impl, (4000,  1.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  sa = 0;
  gsl_sf_legendre_Pl_array_impl(100, 0.5, L);
  sa += ( test_sf_frac_diff(L[0],    1.0 ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff(L[10],  -0.18822860717773437500 ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff(L[100], -0.06051802596186118687 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_sf_legendre_Pl_array_impl(100)");
  s += sa;

  TEST_SF(s, gsl_sf_legendre_Plm_impl, (10, 0, -0.5, &r), -0.18822860717773437500, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_impl, (10, 0, 1.0e-08, &r), -0.24609374999999864648, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_impl, (10, 0, 0.5, &r), -0.18822860717773437500, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Plm_impl, (10, 1, -0.5, &r), -2.0066877394361256516, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_impl, (10, 1, 1.0e-08, &r), -2.7070312499999951725e-07, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_impl, (10, 1, 0.5, &r), 2.0066877394361256516, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Plm_impl, (10, 5, -0.5, &r), -30086.169706116174977, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_impl, (10, 5, 1.0e-08, &r), -0.0025337812499999964949, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_impl, (10, 5, 0.5, &r), 30086.169706116174977, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_impl, (10, 5, 0.999, &r), -0.5036411489013270406, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Plm_impl, (100, 5, -0.5, &r), -6.617107444248382171e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_impl, (100, 5, 1.0e-08, &r), 817.8987598063712851, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_impl, (100, 5, 0.5, &r), 6.617107444248382171e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_impl, (100, 5, 0.999, &r), -1.9831610803806212189e+09, TEST_TOL0, GSL_SUCCESS);

  sa = 0;
  gsl_sf_legendre_Plm_array_impl(100, 5, 0.5, L);
  sa += ( test_sf_frac_diff(L[0],  -460.3466286991656682 ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff(L[10],  38852.51334152290535 ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff(L[95],  6.617107444248382171e+08 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_sf_legendre_Plm_array_impl(100, 5, 0.5)");
  s += sa;

  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (10, 0, -0.5, &r), -0.24332702369300133776, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (10, 0, 0.5, &r), -0.24332702369300133776, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (10, 0, 0.999, &r), 1.2225754122797385990, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (10, 5, -0.5, &r), -0.3725739049803293972, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (10, 5, 1.0e-08, &r), -3.1377233589376792243e-08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (10, 5, 0.5, &r), 0.3725739049803293972, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (10, 5, 0.999, &r), -6.236870674727370094e-06, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (10, 10, -0.5, &r), 0.12876871185785724117, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (10, 10, 0.5, &r), 0.12876871185785724117, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (10, 10, 0.999, &r), 1.7320802307583118647e-14, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (200, 1, -0.5, &r), 0.3302975570099492931, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (200, 1, 0.5, &r), -0.3302975570099492931, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_impl, (200, 1, 0.999, &r), -1.4069792055546256912, TEST_TOL0, GSL_SUCCESS);

  sa = 0;
  gsl_sf_legendre_sphPlm_array_impl(100, 5, 0.5, L);
  sa += ( test_sf_frac_diff(L[0],   -0.22609703187800460722 ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff(L[10],   0.07452710323813558940 ) > TEST_TOL0 );
  sa += ( test_sf_frac_diff(L[95],   0.25865355990880161717 ) > TEST_TOL0 );
  gsl_test(s, "  gsl_sf_legendre_sphPlm_array_impl(100, 5, 0.5)");
  s += sa;

  TEST_SF(s, gsl_sf_conicalP_half_impl, (0.0, -0.5, &r),   0.8573827581049917129, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (0.0,  0.5, &r),   0.8573827581049917129, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (0.0,  2.0, &r),   0.6062611623284649811, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (0.0,  100.0, &r), 0.07979045091636735635, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_half_impl, (10.0, -0.5, &r),    5.345484922591867188e+08 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (10.0,  0.5, &r),    15137.910380385258370 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (10.0,  2.0, &r),    0.4992680691891618544 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (10.0,  100.0, &r), -0.07272008163718195685, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_half_impl, (200.0, -1.0e-3, &r),    1.3347639529084185010e+136, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (200.0,  1.0e-8, &r),    1.0928098010940058507e+136, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (200.0,  0.5, &r),    3.895546021611205442e+90, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (200.0,  10.0, &r),  -0.04308567180833581268, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (200.0,  100.0, &r),   -0.04694669186576399194, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (200.0,  1000.0, &r),   0.023698140704121273277, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (200.0,  1.0e+8, &r),  -0.00006790983312124277891, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_half_impl, (1.0e+8,  1.1, &r),   1.1599311133054742944, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_impl, (1.0e+8,  100.0, &r), 0.07971967557381557875, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (0.0, -0.5, &r),  1.7956982494514644808, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (0.0,  0.5, &r),  0.8978491247257322404, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (0.0,  2.0, &r),  0.7984204253272901551, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (0.0,  100.0, &r),  0.4227531369388072584, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (10.0, -0.5, &r),  5.345484922591867181e+07, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (10.0,  0.5, &r),  1513.7910356104985334, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (10.0,  2.0, &r),  0.03439243987215615642, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (10.0,  100.0, &r),  0.003283756665952609624, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (200.0, -0.5, &r),  1.7699538115312304280e+179, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (200.0,  1.0e-8, &r),  5.464049005470029253e+133, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (200.0,  0.5, &r),  1.9477730108056027211e+88, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (200.0,  10.0, &r),  0.0012462575917716355362, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (200.0,  100.0, &r),  -0.0003225881344802625149, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (200.0,  1000.0, &r), -0.00004330652890886567623 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (200.0,  1.0e+8, &r),  2.0943091278037078483e-07, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (1.0e+8,  1.1, &r), 2.092320445620989618e-09, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_impl, (1.0e+8,  100.0, &r),  -3.359967833599016923e-11, TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_conicalP_0_impl, (0.0, -0.5, &r),  1.3728805006183501647, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_impl, (0.0,  0.5, &r),  1.0731820071493643751, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_impl, (0.0,  2.0, &r),  0.9012862993604472987, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_impl, (0.0,  100.0, &r),  0.30091748588199264556, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_0_impl, (10.0, -0.5, &r),  1.6795592815421804669e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_impl, (10.0,  0.5, &r),  4826.034132009618240, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_impl, (10.0,  2.0, &r),  0.18798468917758716146, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_impl, (10.0,  100.0, &r), -0.008622130749987962529 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_0_impl, (200.0,  -0.5, &r), 2.502194818646823e+180, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_0_impl, (1000.0,  100.0, &r),   0.0017908817653497715844, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_impl, (1000.0,  1000.0, &r), -0.0006566893804926284301, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_impl, (1000.0,  1.0e+8, &r),  2.3167213561756390068e-06, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_1_impl, (0.0, -0.5, &r),    0.4939371126656998499, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (0.0,  0.5, &r),    0.14933621085538265636, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (0.0,  2.0, &r),   -0.13666874968871549533, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (0.0,  100.0, &r), -0.10544528203156629098, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_1_impl, (10.0, -0.5, &r),   1.7253802958788312520e+09, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (10.0,  0.5, &r),   46781.02294059967988, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (10.0,  2.0, &r),   0.26613342643657444400, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (10.0,  100.0, &r), -0.23281959695501029796, TEST_TOL0, GSL_SUCCESS);


  /* FIXME: Mathematica gets some brain-damaged numbers for
   * these x < 0 points. I have checked what I am doing in detail,
   * and it must be right because you can do it by summing
   * manifestly positive definite quantities.
   */
  TEST_SF(s, gsl_sf_conicalP_1_impl, (200.0, -0.999, &r),  2.71635193070012709e+270 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (200.0, -0.9, &r),   4.29524931765857131e+234, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (200.0, -0.5, &r),   5.01159205956053439e+182, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (200.0,  0.999, &r),  195733.0396081538, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (200.0,  10.0, &r), -2.9272610662414349553, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_1_impl, (1000.0,  100.0, &r),   -1.7783258105862399857, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (1000.0,  1000.0, &r),   0.4535161075156427179, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_impl, (1000.0,  1.0e+8, &r),   0.0009983414549874888478, TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_conicalP_sph_reg_impl, (2, 1.0, -0.5, &r),  1.6406279287008789526 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_impl, (10, 1.0, -0.5, &r),  0.000029315266725049129448 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_impl, (20, 1.0, -0.5, &r),  7.335769429462034431e-15 , TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_impl, (30, 1.0, -0.5, &r),  1.3235612394267378871e-26 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_impl, (10, 1.0, 0.5, &r),  2.7016087199857873954e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_impl, (20, 1.0, 0.5, &r),  1.1782569701435933399e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_impl, (30, 1.0, 0.5, &r),  3.636240588303797919e-41 , TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_impl, (10, 1.0, 2.0, &r),  2.4934929626284934483e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_impl, (20, 1.0, 2.0, &r),  1.1284762488012616191e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_impl, (30, 100.0, 100.0, &r),  -1.6757772087159526048e-64 , TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_conicalP_cyl_reg_impl, (2, 1.0, -0.5, &r),  2.2048510472375258708 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_impl, (10, 1.0, -0.5, &r),  0.00007335034531618655690 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_impl, (20, 1.0, -0.5, &r),  2.5419860619212164696e-14 , TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_impl, (30, 1.0, -0.5, &r),  5.579714972260536827e-26 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_impl, (10, 1.0, 0.5, &r),  1.1674078819646475282e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_impl, (20, 1.0, 0.5, &r),  7.066408031229072207e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_impl, (30, 1.0, 0.5, &r),  2.6541973286862588488e-40 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_impl, (10, 1.0, 2.0, &r),  1.0736109751890863051e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_impl, (20, 1.0, 2.0, &r),  6.760965304863386741e-24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_impl, (30, 100.0, 100.0, &r),  -4.268753482520651007e-63 , TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_legendre_H3d_0_impl, (1.0e-06, 1.0e-06, &r), 0.9999999999998333333     , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_impl, (1.0, 0.0, &r), 1.0                       , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_impl, (1.0, 1.0, &r), 0.7160229153604338713     , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_impl, (1.0, 100.0, &r), -3.767437313149604566e-44  , TEST_TOL0, GSL_SUCCESS);  
  TEST_SF(s, gsl_sf_legendre_H3d_0_impl, (1.0, 500.0, &r), -6.665351935878582205e-218 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_impl, (100.0, 1.0, &r), -0.004308757035378200029   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_impl, (100.0, 10.0, &r), 7.508054627912986427e-07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_impl, (1000.0, 1.0, &r), 0.0007036067909088818319  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_impl, (1.0e+08, 1.0, &r), 7.927485371429105968e-09  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_impl, (1.0e+08, 100.0, &r), -3.627118904186918957e-52  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (1.0e-06, 1.0e-06, &r), 3.333333333334222222e-07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (1.0, 1.0e-10, &r), 4.714045207910316829e-11  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (1.0, 1.0, &r), 0.3397013994799344639     , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (1.0, 100.0, &r), -7.200624449531811272e-44  , TEST_TOL0, GSL_SUCCESS);  
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (1.0, 500.0, &r), 4.192260336821728677e-218 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (100.0, 0.01, &r), 0.30117664944267412324    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (100.0, 1.0, &r), -0.007393833425336299309   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (100.0, 10.0, &r), -5.031062029821254982e-07  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (1000.0, 0.001, &r), 0.30116875865090396421    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (1000.0, 1.0, &r), -0.0004776144516074971885  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (1.0e+08, 1.0e-08, &r), 0.30116867893975679722    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (1.0e+08, 1.0, &r), 3.0921097047369081582e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_impl, (1.0e+08, 100.0, &r), -6.496142701296286936e-52  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1.0e-06, 1.0e-06, &r),  1.1544011544013627977e-32 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1.0, 1.0e-10, &r),  2.0224912016958766992e-52 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1.0, 1.0, &r),  0.011498635037491577728 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1.0, 5.0, &r),  0.0020696945662545205776 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1.0, 7.0, &r),  -0.0017555303787488993676 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1.0, 10.0, &r),  0.00008999979724504887101 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1.0, 100.0, &r),  -4.185397793298567945e-44 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1.0, 500.0, &r),  1.4235113901091961263e-217 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 100.0, 0.001, &r),  9.642762597222417946e-10 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 100.0, 0.002, &r),  3.0821201254308036109e-08 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 100.0, 0.01, &r),  0.00009281069019005840532 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 100.0, 1.0, &r),  -0.008043100696178624653 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 100.0, 10.0, &r),  -3.927678432813974207e-07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1000.0, 0.001, &r),  0.00009256365284253254503 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1000.0, 0.01, &r),  -0.05553733815473079983 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1.0e+08, 1.0e-08, &r),   0.00009256115861125841299 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_impl, (5, 1.0e+08, 100.0, &r),    -6.496143209092860765e-52  , TEST_TOL0, GSL_SUCCESS);

#if 0
  sa = 0;
  gsl_sf_legendre_H3d_array_impl(100, 1.0, 3.0, L);
  sa += ( test_sf_frac_diff(L[  0], gsl_sf_legendre_H3d(  0, 1.0, 3.0)) > 1.0e-12 );
  sa += ( test_sf_frac_diff(L[  1], gsl_sf_legendre_H3d(  1, 1.0, 3.0)) > 1.0e-12 );
  sa += ( test_sf_frac_diff(L[ 10], gsl_sf_legendre_H3d( 10, 1.0, 3.0)) > 1.0e-12 );
  sa += ( test_sf_frac_diff(L[100], gsl_sf_legendre_H3d(100, 1.0, 3.0)) > 1.0e-12 );
  gsl_test(sa, "  gsl_sf_legendre_H3d_array_impl(100, 1.0, 3.0)");
  s += sa;
#endif

  TEST_SF(s, gsl_sf_legendre_Q0_impl, (-0.5, &r), -0.5493061443340548457 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_impl, ( 1.5, &r), 0.8047189562170501873 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Q1_impl, (-0.5, &r), -0.7253469278329725772  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_impl, ( 1.5, &r), 0.20707843432557528095 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Ql_impl, (10, -0.5, &r), -0.29165813966586752393 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_impl, (10,  0.5, &r), 0.29165813966586752393 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_impl, (10,  1.5, &r), 0.000014714232718207477406 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Ql_impl, (100, -0.5, &r), -0.09492507395207282096 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_impl, (100,  0.5, &r), 0.09492507395207282096 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_impl, (100,  1.5, &r), 1.1628163435044121988e-43 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Ql_impl, (1000, -0.5, &r), -0.030105074974005303500 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_impl, (1000,  0.5, &r), 0.030105074974005303500 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_impl, (1000,  1.1, &r), 1.0757258447825356443e-194 , TEST_TOL0, GSL_SUCCESS);

  return s;
}


int test_log(void)
{
  gsl_sf_result r;
  gsl_sf_result r1, r2;
  int s = 0;

  TEST_SF(s, gsl_sf_log_impl, (0.1, &r), -2.3025850929940456840  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_impl, (1.1, &r), 0.09531017980432486004 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_impl, (1000.0, &r), 6.907755278982137052   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log_abs_impl, (-0.1, &r), -2.3025850929940456840  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_abs_impl, (-1.1, &r), 0.09531017980432486004 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_abs_impl, (-1000.0, &r), 6.907755278982137052   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_abs_impl, (0.1, &r), -2.3025850929940456840  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_abs_impl, (1.1, &r), 0.09531017980432486004 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log_abs_impl, (1000.0, &r), 6.907755278982137052   , TEST_TOL0, GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_log_impl, (1.0, 1.0, &r1, &r2),
            0.3465735902799726547, TEST_TOL0,
	    0.7853981633974483096, TEST_TOL0,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_log_impl, (1.0, -1.0, &r1, &r2),
             0.3465735902799726547, TEST_TOL0,
	    -0.7853981633974483096, TEST_TOL0,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_log_impl, (1.0, 100.0, &r1, &r2),
            4.605220183488258022, TEST_TOL0,
	    1.560796660108231381, TEST_TOL0,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_log_impl, (-1000.0, -1.0, &r1, &r2),
             6.907755778981887052, TEST_TOL0,
	    -3.1405926539231263718, TEST_TOL0,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_log_impl, (-1.0, 0.0, &r1, &r2),
             0.0, TEST_TOL0,
	     3.1415926535897932385, TEST_TOL0,
             GSL_SUCCESS);


  TEST_SF(s,  gsl_sf_log_1plusx_impl, (1.0e-10, &r), 9.999999999500000000e-11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_impl, (1.0e-8, &r), 9.999999950000000333e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_impl, (1.0e-4, &r), 0.00009999500033330833533 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_impl, (0.1, &r), 0.09531017980432486004 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_impl, (0.49, &r), 0.3987761199573677730 , TEST_TOL0, GSL_SUCCESS);
  
  TEST_SF(s,  gsl_sf_log_1plusx_impl, (-0.49, &r), -0.6733445532637655964 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_impl, (1.0, &r), M_LN2 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_impl, (-0.99, &r), -4.605170185988091368 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_log_1plusx_mx_impl, (1.0e-10, &r), -4.999999999666666667e-21 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_impl, (1.0e-8, &r), -4.999999966666666917e-17 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_impl, (1.0e-4, &r), -4.999666691664666833e-09 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_impl, (0.1, &r), -0.004689820195675139956 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_impl, (0.49, &r), -0.09122388004263222704 , TEST_TOL0, GSL_SUCCESS);
  
  TEST_SF(s,  gsl_sf_log_1plusx_mx_impl, (-0.49, &r), -0.18334455326376559639 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_impl, (1.0, &r), M_LN2-1.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_log_1plusx_mx_impl, (-0.99, &r), -3.615170185988091368 , TEST_TOL0, GSL_SUCCESS);

  return s;
}

int test_poly(void)
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
  s += ( test_sf_frac_diff(y, 1 + 0.5*x + 0.3*x*x) > TEST_TOL0 );
  gsl_test(s, "  gsl_sf_poly_eval({1, 0.5, 0.3}, 0.5)");
  status += s;
  
  s = 0;
  x = 1.0;
  y = gsl_sf_poly_eval(d, 11, x);
  s += ( test_sf_frac_diff(y, 1.0) > TEST_TOL0 );
  gsl_test(s, "  gsl_sf_poly_eval({1,-1, 1, -1, 1, -1, 1, -1, 1, -1, 1}, 1.0)");
  status += s;

  return status;
}

int test_pow_int(void)
{
  gsl_sf_result r;
  int status = 0;
  int s = 0;
  
  TEST_SF(s,  gsl_sf_pow_int_impl, (2.0, 3, &r), 8.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (-2.0, 3, &r), -8.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (2.0, -3, &r), 1.0/8.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (-2.0, -3, &r), -1.0/8.0 , TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s,  gsl_sf_pow_int_impl, (10.0, 4, &r), 1.0e+4 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (10.0, -4, &r), 1.0e-4 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (-10.0, 4, &r), 1.0e+4 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (-10.0, -4, &r), 1.0e-4 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_pow_int_impl, (10.0, 40, &r), 1.0e+40 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (8.0, -40, &r), 7.523163845262640051e-37 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (-10.0, 40, &r), 1.0e+40 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (-8.0, -40, &r), 7.523163845262640051e-37 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_pow_int_impl, (10.0, 41, &r), 1.0e+41 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (8.0, -41, &r), 9.403954806578300064e-38 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (-10.0, 41, &r), -1.0e+41 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pow_int_impl, (-8.0, -41, &r), -9.403954806578300064e-38 , TEST_TOL0, GSL_SUCCESS);

  return status;
}

int test_psi(void)
{
  gsl_sf_result r;
  int s = 0;
  
  TEST_SF(s, gsl_sf_psi_int_impl, (5, &r), 1.5061176684318004727 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_int_impl, (100, &r), 4.600161852738087400 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_int_impl, (110, &r), 4.695928024251535633 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_int_impl, (5000, &r), 8.517093188082904107 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_psi_impl, (5.0, &r), 1.5061176684318004727 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_impl, (5000.0, &r), 8.517093188082904107 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_psi_impl, (-100.5, &r), 4.615124601338064117 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_impl, (-1.0e+5-0.5, &r), 11.512935464924395337 , TEST_TOL0, GSL_SUCCESS);
  
  TEST_SF(s, gsl_sf_psi_1piy_impl, (0.8, &r), -0.07088340212750589223 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1piy_impl, (1.0, &r), 0.09465032062247697727 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1piy_impl, (5.0, &r), 1.6127848446157465854  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1piy_impl, (100.0, &r), 4.605178519404762003 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1piy_impl, (2000.0, &r), 7.600902480375416216 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_psi_1_int_impl, (5, &r), 0.22132295573711532536 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1_int_impl, (100, &r), 0.010050166663333571395 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1_int_impl, (110, &r), 0.009132356622022545705 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_1_int_impl, (500, &r), 0.0020020013333322666697 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_psi_n_impl, (3, 5.0, &r), 0.021427828192755075022 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_n_impl, (3, 500.0, &r), 1.6048063999872000683e-08  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_n_impl, (10, 5.0, &r), -0.08675107579196581317 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_psi_n_impl, (10, 50.0, &r), -4.101091112731268288e-12 , TEST_TOL0, GSL_SUCCESS);

  return s;
}


int test_synch(void)
{
  gsl_sf_result r;
  int s = 0;

  TEST_SF(s, gsl_sf_synchrotron_1_impl, (0.01, &r),  0.444973 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_synchrotron_1_impl, (1.0, &r),  0.651423 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_synchrotron_1_impl, (10.0, &r),  0.000192238 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_synchrotron_1_impl, (100.0, &r),  4.69759e-43 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_synchrotron_2_impl, (0.01, &r),  0.23098077342226277732 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_synchrotron_2_impl, (1.0, &r),  0.4944750621042082670 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_synchrotron_2_impl, (10.0, &r),  0.00018161187569530204281 , TEST_TOL0, GSL_SUCCESS);  
  TEST_SF(s, gsl_sf_synchrotron_2_impl, (256.0, &r),  1.3272635474353774058e-110 , TEST_TOL0, GSL_SUCCESS);  /* exp()... not my fault */

  return s;
}


int test_transport(void)
{
  gsl_sf_result r;
  int s;

  TEST_SF(s, gsl_sf_transport_2_impl, (1.0e-10, &r), 9.9999999999999999999e-11 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_2_impl, (1.0, &r), 0.97303256135517012845 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_2_impl, (3.0, &r), 2.41105004901695346199 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_2_impl, (10.0, &r), 3.28432911449795173575 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_2_impl, (100.0, &r), 3.28986813369645287294 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_2_impl, (1.0e+05, &r), 3.28986813369645287294 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_transport_3_impl, (1.0e-10, &r), 4.999999999999999999997e-21 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_impl, (1.0, &r), 0.479841006572417499939 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_impl, (3.0, &r), 3.210604662942246772338 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_impl, (5.0, &r), 5.614386613842273228585 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_impl, (10.0, &r), 7.150322712008592975030 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_impl, (30.0, &r), 7.212341416160946511930 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_impl, (100.0, &r), 7.212341418957565712398 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_3_impl, (1.0e+05, &r), 7.212341418957565712398 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_transport_4_impl, (1.0e-10, &r), 3.33333333333333333333e-31 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_impl, (1.0e-07, &r), 3.33333333333333166666e-22 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_impl, (1.0e-04, &r), 3.33333333166666666726e-13 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_impl, (0.1, &r), 0.000333166726172109903824 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_impl, (1.0, &r), 0.31724404523442648241 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_impl, (3.0, &r), 5.96482239737147652446 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_impl, (5.0, &r), 15.3597843168821829816 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_impl, (10.0, &r), 25.2736676770304417334 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_impl, (30.0, &r), 25.9757575220840937469 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_impl, (100.0, &r), 25.9757576090673165963 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_4_impl, (1.0e+05, &r), 25.9757576090673165963 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_transport_5_impl, (1.0e-10, &r), 2.49999999999999999999e-41 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_impl, (1.0e-07, &r), 2.49999999999999861111e-29 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_impl, (1.0e-04, &r), 2.49999999861111111163e-17 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_impl, (0.1, &r), 0.000024986116317791487410 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_impl, (1.0, &r), 0.236615879239094789259153 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_impl, (3.0, &r), 12.77055769104415951115760 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_impl, (5.0, &r), 50.26309221817518778543615 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_impl, (10.0, &r), 116.3807454024207107698556 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_impl, (30.0, &r), 124.4313279083858954839911 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_impl, (100.0, &r), 124.4313306172043911597639 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_transport_5_impl, (1.0e+05, &r), 124.43133061720439115976   , TEST_TOL0, GSL_SUCCESS);

  return s;
}


int test_trig(void)
{
  gsl_sf_result r;
  gsl_sf_result r1, r2;
  double theta;
  int s = 0;
  int sa;

  TEST_SF(s, gsl_sf_sin_pi_x_impl, (1000.5, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_pi_x_impl, (10000.0 + 1.0/65536.0, &r), 0.00004793689960306688455, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_pi_x_impl, (1099511627776.0 + 1 + 0.125, &r), -0.3826834323650897717, TEST_TOL0, GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_sin_impl, (1.0, 5.0, &r1, &r2),
            62.44551846769653403, TEST_TOL0,
            40.09216577799840254, TEST_TOL0,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_cos_impl, (1.0, 5.0, &r1, &r2),
             40.09580630629882573, TEST_TOL0,
            -62.43984868079963017, TEST_TOL0,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_logsin_impl, (1.0, 100.0, &r1, &r2),
            99.3068528194400546900, TEST_TOL0,
            0.5707963267948966192, TEST_TOL0,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_logsin_impl, (1.0, -100.0, &r1, &r2),
             99.3068528194400546900, TEST_TOL0,
            -0.5707963267948966192, TEST_TOL0,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_logsin_impl, (5.0, 5.0, &r1, &r2),
            4.3068909128079757420, TEST_TOL0,
            2.8540063315538773952, TEST_TOL0,
            GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_lnsinh_impl, (0.1, &r),  -2.3009189815304652235,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnsinh_impl, (1.0, &r),   0.16143936157119563361, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnsinh_impl, (5.0, &r),   4.306807418479684201,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnsinh_impl, (100.0, &r), 99.30685281944005469,   TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_lncosh_impl, (0.125, &r), 0.007792239318898252791, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lncosh_impl, (1.0, &r),   0.4337808304830271870,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lncosh_impl, (5.0, &r),   4.306898218339271555, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lncosh_impl, (100.0, &r), 99.30685281944005469, TEST_TOL0, GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_polar_to_rect_impl, (10.0, M_PI/6.0, &r1, &r2),
            (10.0 * sqrt(3) / 2.0), TEST_TOL0,
	    (10.0 * 0.5), TEST_TOL0,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_polar_to_rect_impl, (10.0, -2.0/3.0*M_PI, &r1, &r2),
            (10.0 * (-0.5)), TEST_TOL0,
	    (10.0 * (-sqrt(3) / 2.0)), TEST_TOL0,
            GSL_SUCCESS);


  theta = 5.0*M_PI + M_PI/2.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, 3.0/2.0*M_PI ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_pos_impl: theta =  11/2 Pi");
  s += sa;

  theta = -5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, M_PI/2.0 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_pos_impl: theta = -11/2 Pi");
  s += sa;

  theta = 50000.0 + 1.0/65536.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, 4.6945260308194656055 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_pos_impl: theta = 50000.0 + 1.0/65536.0");
  s += sa;

  theta = 5000000.0 + 1.0/65536.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, 4.49537973053997376 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_pos_impl: theta = 5000000.0 + 1.0/65536.0");
  s += sa;

  theta = 140737488355328.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, 3.20652300406795792638 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_pos_impl: theta = 2^47");
  s += sa;

  theta = 5.0*M_PI + M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, -M_PI/2.0 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_symm_impl: theta =  11/2 Pi");
  s += sa;

  theta = -5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, M_PI/2.0 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_symm_impl: theta = -11/2 Pi");
  s += sa;

  theta =  5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, M_PI/2.0 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_symm_impl: theta = -9/2 Pi");
  s += sa;

  theta =  3.0/2.0*M_PI;
  gsl_sf_angle_restrict_symm_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, -M_PI/2.0 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_symm_impl: theta =  3/2 Pi");
  s += sa;

  theta = -3.0/2.0*M_PI;
  gsl_sf_angle_restrict_symm_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, M_PI/2.0 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_symm_impl: theta = -3/2 Pi");
  s += sa;

  theta = 50000.0 + 1.0/65536.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  sa = 0;
  sa += ( test_sf_frac_diff( theta, -1.5886592763601208714 ) > TEST_TOL0 );
  gsl_test(sa, "  gsl_angle_restrict_symm_impl: theta = 50000.0 + 1.0/65536.0");
  s += sa;

  return s;
}


int test_zeta(void)
{
  gsl_sf_result r;
  int s = 0;

  TEST_SF(s, gsl_sf_zeta_int_impl, (-61, &r), -3.30660898765775767257e+34, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_int_impl, (-5, &r), -0.003968253968253968253968, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_int_impl, (5, &r), 1.0369277551433699263313655, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_int_impl, (31, &r), 1.0000000004656629065033784, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_zeta_impl, (-151, &r), 8.195215221831378294e+143  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_impl, (-51, &r), 9.68995788746359406565e+24 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_impl, (-5, &r), -0.003968253968253968253968 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_impl, (-0.5, &r), -0.207886224977354566017307 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_impl, (0.5, &r), -1.460354508809586812889499 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_impl, (1.0-1.0/1024.0, &r), -1023.4228554489429787      , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_impl, (1.0+1.0/1048576, &r), 1.0485765772157343441e+06  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_impl, (5, &r), 1.036927755143369926331365 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_zeta_impl, (25.5, &r), 1.000000021074106110269959 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hzeta_impl, (2,  1.0, &r), 1.6449340668482264365   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_impl, (2, 10.0, &r), 0.1051663356816857461   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_impl, (5,  1.0, &r), 1.0369277551433699263   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_impl, (5, 10.0, &r), 0.000030413798676470276 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_impl, (9,  0.1, &r), 1.0000000004253980e+09  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_impl, (30, 0.5, &r), 1.0737418240000053e+09  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_impl, (30, 0.9, &r), 2.3589824880264765e+01  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hzeta_impl, (75, 0.25, &r), 1.4272476927059599e+45  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_eta_int_impl, (-91, &r), -4.945598888750002040e+94 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_impl, (-51, &r), -4.363969073121683116e+40 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_impl, (-5, &r), 0.25 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_impl, (-1, &r), 0.25 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_impl, ( 0, &r), 0.5  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_impl, ( 5, &r), 0.9721197704469093059 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_impl, ( 6, &r), 0.9855510912974351041 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_impl, ( 20, &r), 0.9999990466115815221 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_int_impl, ( 1000, &r), 1.0 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_eta_impl, (-51.5, &r), -1.2524184036924703656e+41 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_impl, (-5, &r), 0.25 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_impl, (0.5, &r), 0.6048986434216303702 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_impl, (0.999, &r), 0.6929872789683383574 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_impl, (1.0, &r), 0.6931471805599453094 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_impl, (1.0+1.0e-10, &r), 0.6931471805759321998 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_impl, ( 5, &r), 0.9721197704469093059 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_impl, ( 5.2, &r), 0.9755278712546684682 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_impl, ( 6, &r), 0.9855510912974351041 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_eta_impl, ( 20, &r), 0.9999990466115815221 , TEST_TOL0, GSL_SUCCESS);

  return s;
}


int main(int argc, char * argv[])
{
  gsl_test(test_airy(),       "Airy Functions");
  gsl_test(test_bessel(),     "Bessel Functions");
  gsl_test(test_cheb(),       "Chebyshev Evaluation");
  gsl_test(test_clausen(),    "Clausen Integral");
  gsl_test(test_coulomb(),    "Coulomb Wave Functions");
  gsl_test(test_coupling(),   "Coupling Coefficients");
  gsl_test(test_dawson(),     "Dawson Integral");
  gsl_test(test_debye(),      "Debye Functions");
  gsl_test(test_dilog(),      "Dilogarithm");
  gsl_test(test_elementary(), "Elementary Functions (Misc)");
  gsl_test(test_ellint(),     "Elliptic Integrals");
  gsl_test(test_jac(),        "Elliptic Functions (Jacobi)");
  gsl_test(test_erf(),        "Error Functions");
  gsl_test(test_exp(),        "Exponential Functions");
  gsl_test(test_expint(),     "Exponential/Sine/Cosine Integrals");
  gsl_test(test_fermidirac(), "Fermi-Dirac Functions");
  gsl_test(test_gamma(),      "Gamma Functions");
  gsl_test(test_gegen(),      "Gegenbauer Polynomials");
  gsl_test(test_hyperg(),     "Hypergeometric Functions");
  gsl_test(test_laguerre(),   "Laguerre Polynomials");
  gsl_test(test_legendre(),   "Legendre Functions");
  gsl_test(test_log(),        "Logarithm");
  gsl_test(test_poly(),       "Polynomial Evaluation");
  gsl_test(test_pow_int(),    "Integer Powers");
  gsl_test(test_psi(),        "Psi Functions");
  gsl_test(test_synch(),      "Synchrotron Functions");
  gsl_test(test_transport(),  "Transport Functions");
  gsl_test(test_trig(),       "Trigonometric and Related Functions");
  gsl_test(test_zeta(),       "Zeta Functions");

  return gsl_test_summary();
}
