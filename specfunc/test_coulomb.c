/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_test.h>
#include <gsl_sf.h>
#include "test_sf.h"

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

  TEST_SF(s, gsl_sf_hydrogenicR_1_impl, (3.0, 2.0, &r),  0.025759948256148471036,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_1_impl, (3.0, 10.0, &r), 9.724727052062819704e-13, TEST_TOL0, GSL_SUCCESS);  status += s;

  TEST_SF(s, gsl_sf_hydrogenicR_impl, (4, 0, 3.0, 2.0, &r), -0.03623182256981820062,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_impl, (4, 1, 3.0, 2.0, &r), -0.028065049083129581005, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_impl, (4, 2, 3.0, 2.0, &r),  0.14583027278668431009,  TEST_TOL0, GSL_SUCCESS);  status += s;

  TEST_SF(s, gsl_sf_hydrogenicR_impl, (100,  0, 3.0, 2.0, &r), -0.00007938950980052281367, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_impl, (100, 10, 3.0, 2.0, &r),  7.112823375353605977e-12,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_impl, (100, 90, 3.0, 2.0, &r),  5.845231751418131548e-245, TEST_TOL2, GSL_SUCCESS);

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

  return s;
}
