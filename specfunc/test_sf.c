/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_test.h>
#include <gsl_sf.h>


double frac_diff(double x1, double x2)
{
  if(x1 <= DBL_MAX && x2 <= DBL_MAX)
    return fabs((x1-x2)/(x1+x2));
  else
    return 1.0;
}


int check_airy(void)
{
  int status = 0;
  int s;

  /** functions */

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Ai(-5),                   0.3507610090241142    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(-0.3000000000000094),  0.4309030952855831    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(0.6999999999999907),   0.1891624003981519    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(1.649999999999991),    0.05831058618720882   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(2.54999999999999),     0.01446149513295428   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(3.499999999999987),    0.002584098786989702  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(5.39999999999998),     4.272986169411866e-05 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Ai");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(-5),            gsl_sf_airy_Ai(-5)         ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(-0.5),          gsl_sf_airy_Ai(-0.5)       ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(0.6999999999999907),   0.2795125667681217  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(1.649999999999991),    0.2395493001442741  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(2.54999999999999),     0.2183658595899388  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(3.499999999999987),    0.2032920808163519  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(5.39999999999998),     0.1836050093282229  ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Ai_scaled");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_airy_Bi(-5),                   -0.1383691349016005    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(-0.3000000000000094),   0.4779778401098885    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(0.6999999999999907),    0.9733286558781599    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(1.649999999999991),     2.196407956850028     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(2.54999999999999),      6.973628612493443     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(3.499999999999987),     33.05550675461069     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(5.39999999999998),      1604.476078241272     ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Bi");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(-5),            gsl_sf_airy_Bi(-5)        ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(-0.5),          gsl_sf_airy_Bi(-0.5)      ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(0.6999999999999907),  0.6587080754582302  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(1.649999999999991),   0.5346449995597539  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(2.54999999999999),    0.461835455542297   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(3.499999999999987),   0.4201771882353061  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(5.39999999999998),    0.3734050675720473  ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Bi_scaled");
  status += s;
  
  /** derivatives */

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(-5),                   0.3271928185544435)  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(-0.5500000000000094), -0.1914604987143629)  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(0.4999999999999906),  -0.2249105326646850)  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(1.899999999999992),   -0.06043678178575718) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(3.249999999999988),   -0.007792687926790889) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(5.199999999999981),   -0.0001589434526459543) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Ai_deriv");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(-5),          gsl_sf_airy_Ai_deriv(-5)   )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(-0.5),        gsl_sf_airy_Ai_deriv(-0.5) )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(0.5499999999999906), -0.2874057279170166 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(1.499999999999991),  -0.3314199796863637 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(2.49999999999999),   -0.366108938475162  )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(3.649999999999986),  -0.3974033831453963 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(6.299999999999977),  -0.4508799189585947 )  > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Ai_deriv_scaled");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(-5),                   0.778411773001899 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(-0.5500000000000094),  0.5155785358765014)  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(0.4999999999999906),   0.5445725641405883)  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(1.899999999999992),    3.495165862891568 ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(3.249999999999988),    36.55485149250338 ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(5.199999999999981),    2279.748293583233 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Bi_deriv");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(-5),          gsl_sf_airy_Bi_deriv(-5)   )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(-0.5),        gsl_sf_airy_Bi_deriv(-0.5) )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(0.5499999999999906),  0.4322811281817566 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(1.499999999999991),   0.5542307563918037 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(2.49999999999999),    0.6755384441644985 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(3.649999999999986),   0.7613959373000228 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(6.299999999999977),   0.8852064139737571 )  > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Bi_deriv_scaled");
  status += s;
  return status;
}


int check_bessel(void)
{
  const int nmax = 100;
  double J[nmax+1], Jp[nmax+1];
  double Y[nmax+1], Yp[nmax+1];
  double I[nmax+1], Ip[nmax+1];
  double K[nmax+1], Kp[nmax+1];
  int status = 0;
  int s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_J0(0.1),     0.99750156206604003230  ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_J0(2.0),     0.22389077914123566805  ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_J0(100.0),   0.019985850304223122424 ) > 1.0e-13 );
  s += ( frac_diff(gsl_sf_bessel_J0(1.0e+10), 2.1755917502468917269e-06 ) > 1.0e-06 );
  gsl_test(s, "  gsl_sf_bessel_J0");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_bessel_J1(0.1),      0.04993752603624199756 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_J1(2.0),      0.57672480775687338720 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_J1(100.0),   -0.07714535201411215803 ) > 1.0e-13 );
  s += ( frac_diff(gsl_sf_bessel_J1(1.0e+10), -7.676508175684157103e-06 ) > 1.0e-06 );
  gsl_test(s, "  gsl_sf_bessel_J1");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Jn(   4,     0.1), 2.6028648545684032338e-07 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Jn(   5,     2.0), 0.007039629755871685484   ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Jn( 100,   100.0), 0.09636667329586155967    ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Jn(1000, 1.0e+10), 2.1759755729360371995e-06 ) > 1.0e-04 );
  gsl_test(s, "  gsl_sf_bessel_Jn");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Y0(0.1),     -1.5342386513503668441  ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Y0(2.0),      0.5103756726497451196  ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Y0(100.0),   -0.07724431336508315225 ) > 1.0e-13 );
  s += ( frac_diff(gsl_sf_bessel_Y0(1.0e+10), -7.676508175792936690e-06 ) > 1.0e-06 );
  gsl_test(s, "  gsl_sf_bessel_Y0");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Y1(0.1),     -6.45895109470202698800  ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Y1(2.0),     -0.10703243154093754689  ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Y1(100.0),   -0.020372312002759793305 ) > 1.0e-13 );
  s += ( frac_diff(gsl_sf_bessel_Y1(1.0e+10), -2.1755917506307171357e-06 ) > 1.0e-05 );
  gsl_test(s, "  gsl_sf_bessel_Y1");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Yn(   4,     0.1), -305832.29793353160319    ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Yn(   5,     2.0), -9.935989128481974981     ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Yn( 100,   100.0), -0.16692141141757650654   ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Yn(1000, 1.0e+10), -7.676399386609853644e-06 ) > 1.0e-06 );
  gsl_test(s, "  gsl_sf_bessel_Yn");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_I0(0.1),   1.0025015629340956014      ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_I0(2.0),   2.2795853023360672674      ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_I0(100.0), 1.0737517071310738235e+42  ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_bessel_I0");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_bessel_I1(0.1),   0.05006252604709269211    ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_I1(2.0),   1.59063685463732906340    ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_I1(100.0), 1.0683693903381624812e+42 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_bessel_I1");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_In(   4,    0.1), 2.6054690212996573677e-07  ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_In(   5,    2.0), 0.009825679323131702321    ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_In( 100,  100.0), 4.641534941616199114e+21   ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_bessel_In");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_K0(0.1),   2.4270690247020166125    ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_K0(2.0),   0.11389387274953343565   ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_K0(100.0), 4.656628229175902019e-45 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_bessel_K0");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_bessel_K1(0.1),   9.853844780870606135     ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_K1(2.0),   0.13986588181652242728   ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_K1(100.0), 4.679853735636909287e-45 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_bessel_K1");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Kn(   4,    0.1),  479600.2497925682849     ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Kn(   5,    2.0),  9.431049100596467443     ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_bessel_Kn( 100,  100.0),  7.617129630494085416e-25 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_bessel_Kn");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Jnu(0.0001,10.0), -0.2459270166445205        ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Jnu( 1.0, 0.001),  0.0004999999375000026     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_bessel_Jnu( 1.0,   1.0),  0.4400505857449335160     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_bessel_Jnu(30.0,   1.0),  3.482869794251482902e-42  ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Jnu(30.0, 100.0),  0.08146012958117222297    ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Jnu(10.0,   1.0),  2.6306151236874532070e-10 ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Jnu(10.0, 100.0), -0.05473217693547201474    ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Jnu(10.2, 100.0), -0.03548919161046526864  ) > 1.e-13 );
  gsl_test(s, "  gsl_sf_bessel_Jnu");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Ynu(0.0001,10.0),  0.05570979797521875261    ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Ynu( 1.0, 0.001), -636.6221672311394281      ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_bessel_Ynu( 1.0,   1.0), -0.7812128213002887165     ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Ynu(30.0,   1.0), -3.0481287832256432162e+39 ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Ynu(30.0, 100.0),  0.006138839212010033452   ) > 1.e-10 );
  s += ( frac_diff(gsl_sf_bessel_Ynu(10.0,   1.0), -1.2161801427868918929e+08 ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Ynu(10.0, 100.0),  0.05833157423641492875    ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Ynu(10.2, 100.0),  0.07169383985546287091    ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_bessel_Ynu");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Inu_scaled(0.0001,10.0), 0.12783333709581669672    ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Inu_scaled( 1.0, 0.001), 0.0004995003123542213370  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_bessel_Inu_scaled( 1.0,   1.0), 0.20791041534970844887    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_bessel_Inu_scaled(30.0,   1.0), 1.3021094983785914437e-42 ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Inu_scaled(30.0, 100.0), 0.0004486987756920986146  ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Inu_scaled(10.0,   1.0), 1.0127529864692066036e-10 ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Inu_scaled(10.0, 100.0), 0.024176682718258828365   ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Inu_scaled(10.2, 100.0), 0.023691628843913810043   ) > 1.e-13 );
  gsl_test(s, "  gsl_sf_bessel_Inu_scaled");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Inu(0.0001,10.0), 2815.7166269770030352     ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Inu( 1.0, 0.001), 0.0005000000625000026042  ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Inu( 1.0,   1.0), 0.5651591039924850272     ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Inu(30.0,   1.0), 3.539500588106447747e-42  ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Inu(30.0, 100.0), 1.2061548704498434006e+40 ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Inu(10.0,   1.0), 2.7529480398368736252e-10 ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Inu(10.0, 100.0), 6.498975524720147799e+41  ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Inu(10.2, 100.0), 6.368587361287030443e+41  ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_bessel_Inu");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Knu_scaled(0.0001,10.0), 0.3916319346235421817     ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Knu_scaled( 1.0, 0.001), 1000.9967345590684524     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_bessel_Knu_scaled( 1.0,   1.0), 1.6361534862632582465     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_bessel_Knu_scaled(30.0,   1.0), 1.2792629867539753925e+40 ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Knu_scaled(30.0, 100.0), 10.673443449954850040     ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Knu_scaled(10.0,   1.0), 4.912296520990198599e+08  ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Knu_scaled(10.0, 100.0), 0.20578687173955779807    ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Knu_scaled(10.2, 100.0), 0.20995808355244385075    ) > 1.e-13 );
  gsl_test(s, "  gsl_sf_bessel_Knu_scaled");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Knu(0.0001,10.0), 0.000017780062324654874306 ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Knu( 1.0, 0.001), 999.9962381560855743       ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Knu( 1.0,   1.0), 0.6019072301972345747      ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Knu(30.0,   1.0), 4.706145526783626883e+39   ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Knu(30.0, 100.0), 3.970602055959398739e-43   ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Knu(10.0,   1.0), 1.8071328990102945469e+08  ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Knu(10.0, 100.0), 7.655427977388100611e-45   ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Knu(10.2, 100.0), 7.810600225948217841e-45   ) > 1.e-13 );
  gsl_test(s, "  gsl_sf_bessel_Knu");
  status += s;

  return status;
}

int check_cheb(void)
{
  double x;
  double f;
  int status = 0;
  int s;
  
  struct gsl_sf_cheb_series * cs = gsl_sf_cheb_new(sin, -M_PI, M_PI, 30);

  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    f += fabs(gsl_sf_cheb_eval(cs, x) - sin(x));
  }
  s = 0;
  s += ( f > 100.0 * 1.0e-14 );
  gsl_test(s, "  gsl_sf_cheb_eval()");
  status += s;
  
  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    f += fabs(gsl_sf_cheb_eval_n(cs, 25, x) - sin(x));
  }
  s = 0;
  s += ( f > 100.0 * 1.0e-14 );
  gsl_test(s, "  gsl_sf_cheb_eval_n()");
  status += s;
  
  gsl_sf_cheb_calc_impl(cs, sin);
  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    f += fabs(gsl_sf_cheb_eval(cs, x) - sin(x));
  }
  s = 0;
  s += ( f > 100.0 * 1.0e-14 );
  gsl_test(s, "  gsl_sf_cheb_calc()");
  status += s;

  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    f += fabs(gsl_sf_cheb_eval_deriv(cs, x) - cos(x));
  }
  s = 0;
  s += ( f > 100.0 * 10.0 * 1.0e-14 );
  gsl_test(s, "  gsl_sf_cheb_eval_deriv()");
  status += s;

  f = 0.0;
  for(x=-M_PI; x<M_PI; x += M_PI/100.0) {
    f += fabs(gsl_sf_cheb_eval_integ(cs, x) + (1.0 + cos(x)));
  }
  s = 0;
  s += ( f > 100.0 * 1.0e-14 );
  gsl_test(s, "  gsl_sf_cheb_eval_integ()");
  status += s;

  gsl_sf_cheb_free(cs);

  return status;
}


int check_clausen(void)
{
  double y;
  int status = 0;
  int s;

  s = 0;
  s += ( frac_diff( gsl_sf_clausen(M_PI/20.0), 0.4478882448133546 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_clausen(M_PI/6.0),  0.8643791310538927 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_clausen(M_PI/3.0),  1.0149416064096535 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_clausen(  2.0*M_PI + M_PI/3.0),  1.0149416064096535 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_clausen(100.0*M_PI + M_PI/3.0),  1.0149416064096535 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_clausen");
  status += s;
  
  s = 0;
  s += ( gsl_sf_clausen_impl(1.0e+10*M_PI + M_PI/3.0, &y) != GSL_ELOSS);
  gsl_test(s, "  gsl_sf_clausen: trap accuracy loss from large argument");
  status += s;

  return status;
}


#define PRINT(n) printf("%20.16g  %20.16g  %20.16g  %20.16g    %20.16g\n", F[n], Fp[n], G[n], Gp[n], Fp[n]*G[n]-Gp[n]*F[n])

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
  s += ( frac_diff(  F[0],  1.7207454091787930614e-06 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  3.0975994706405458046e-06 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  167637.56609459967623     ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -279370.76655361803075     ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(10.0, 5.0, lam_F=0, lam_G=0)");
  PRINT(0);
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 25.0;
  x = 10.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  1.5451274501076114315e-16 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  3.1390869393378630928e-16 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.6177129008336318136e+15 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -3.1854062013149740860e+15 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(25.0, 10.0, lam_F=0, lam_G=0)");
  PRINT(0);
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
  s += ( frac_diff(  F[0],  0.0016262711250135878249 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  0.0017060476320792806014 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  307.87321661090837987    ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -291.92772380826822871    ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(10.0, 10.0, lam_F=0, lam_G=0)");
  PRINT(0);
  status += s;

  lam_F = 0.0;
  eta = 100.0;
  x = 1.0;
  gsl_sf_coulomb_wave_FG_impl(eta, x, lam_F, k_G, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  8.999367996930662705e-126 ) > 1.e-3 );
  s += ( frac_diff( Fp[0],  1.292746745757069321e-124 ) > 1.e-3 );
  s += ( frac_diff(  G[0],  3.936654148133683610e+123 ) > 1.e-3 );
  s += ( frac_diff( Gp[0], -5.456942268061526371e+124 ) > 1.e-3 );
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_impl(100.0, 1.0, lam_F=0, lam_G=0)");
  status += s;

  return status;
}

int check_coupling(void)
{
  double y;
  int status = 0;
  int s;

  gsl_sf_coupling_3j_impl(0, 1, 1, 0, 1, -1, &y);
  s = 0;
  s += ( frac_diff( y, sqrt(1.0/2.0) ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_coupling_3j(ja=0    jb=1/2  j=1/2  ma=0    mb=1/2  m=-1/2)");
  status += s;

  gsl_sf_coupling_3j_impl(1, 1, 2, 1, -1, 0, &y);
  s = 0;
  s += ( frac_diff( y, sqrt(1.0/6.0) ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_coupling_3j(ja=1/2  jb=1/2  j=1    ma=1/2  mb=-1/2 m=0)");
  status += s;

  gsl_sf_coupling_3j_impl(2, 4, 6, 0, 2, -2, &y);
  s = 0;
  s += ( frac_diff( y, sqrt(8.0/105.0) ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_coupling_3j(ja=1    jb=2    j=3    ma=0    mb=1    m=-1)");
  status += s;

  gsl_sf_coupling_6j_impl(2, 2, 4, 2, 2, 2, &y);
  s = 0;
  s += ( frac_diff( y, 1.0/6.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_coupling_6j(ja=1   jb=1   jc=2    jd=1     je=1     jf=1)");
  status += s;

  gsl_sf_coupling_9j_impl(4, 2, 4, 3, 3, 2, 1, 1, 2, &y);
  s = 0;
  s += ( frac_diff( y, -0.040824829046386 ) > 1.0e-13 );
  gsl_test(s, "  gsl_sf_coupling_9j(ja=2   jb=1   jc=2    jd=3/2   je=3/2   jf=1   jg=1/2  jh=1/2  ji=1)");
  status += s;

  gsl_sf_coupling_9j_impl(8, 4, 10, 7, 3, 8, 1, 1, 2, &y);
  s = 0;
  s += ( frac_diff( y, 0.025458753860866 ) > 1.0e-13 );
  gsl_test(s, "  gsl_sf_coupling_9j(ja=4   jb=2   jc=5    jd=7/2   je=3/2   jf=4   jg=1/2  jh=1/2  ji=1)");
  status += s;

  return status;
}

int check_dawson(void)
{
  int status = 0;
  int s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_dawson(1.0e-15), 1.0e-15) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dawson(1.0e-15)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_dawson(0.5), 0.4244363835020222959 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dawson(0.5)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_dawson(2.0), 0.30134038892379196603 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dawson(2.0)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_dawson(1000), 0.0005000002500003750009 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dawson(1000)");
  status += s;
  
  return status;
}

int check_debye(void)
{
  double y;
  int status = 0;
  int s;

  s = 0;
  gsl_sf_debye_1_impl(0.1, &y);
  s += ( frac_diff( y, 0.975278 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_1(0.1)");
  status += s;

  s = 0;
  gsl_sf_debye_1_impl(1.0, &y);
  s += ( frac_diff( y, 0.777505 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_1(1.0)");
  status += s;

  s = 0;
  gsl_sf_debye_1_impl(10.0, &y);
  s += ( frac_diff( y, 0.164443 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_1(10.0)");
  status += s;

  s = 0;
  gsl_sf_debye_2_impl(0.1, &y);
  s += ( frac_diff( y, 0.967083 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_2(0.1)");
  status += s;

  s = 0;
  gsl_sf_debye_2_impl(1.0, &y);
  s += ( frac_diff( y, 0.707878 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_2(1.0)");
  status += s;

  s = 0;
  gsl_sf_debye_2_impl(10.0, &y);
  s += ( frac_diff( y, 0.047971 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_2(10.0)");
  status += s;

  s = 0;
  gsl_sf_debye_3_impl(0.1, &y);
  s += ( frac_diff( y, 0.963000 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_3(0.1)");
  status += s;

  s = 0;
  gsl_sf_debye_3_impl(1.0, &y);
  s += ( frac_diff( y, 0.674416 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_3(1.0)");
  status += s;

  s = 0;
  gsl_sf_debye_3_impl(10.0, &y);
  s += ( frac_diff( y, 0.019296 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_3(10.0)");
  status += s;

  s = 0;
  gsl_sf_debye_4_impl(0.1, &y);
  s += ( frac_diff( y, 0.960555 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_4(0.1)");
  status += s;

  s = 0;
  gsl_sf_debye_4_impl(1.0, &y);
  s += ( frac_diff( y, 0.654874 ) > 1.0e-5 );
  gsl_test(s, "  gsl_sf_debye_4(1.0)");
  status += s;

  s = 0;
  gsl_sf_debye_4_impl(10.0, &y);
  s += ( frac_diff( y, 0.009674 ) > 1.0e-4 );
  gsl_test(s, "  gsl_sf_debye_4(10.0)");
  status += s;

  return status;
}

/* FIXME: probably need more tests here... 
 * also need to work on accuracy for r->1; need to
 * adjust the switch-over point I suppose.
 */
int check_dilog(void)
{
  double x, y;
  int status = 0;
  int s;

  /* real dilog */

  s = 0;
  s += ( frac_diff( gsl_sf_dilog(-3.0), -1.9393754207667089531 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dilog(-3.0)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_dilog(-0.5), -0.4484142069236462024 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dilog(-0.5)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_dilog(-0.001), -0.0009997501110486510834 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dilog(-0.001");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_dilog(0.1),  0.1026177910993911 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dilog(0.1)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_dilog(0.7),  0.8893776242860387386 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dilog(0.7)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_dilog(1.0),  1.6449340668482260 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dilog(1.0)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_dilog(1.5),  2.3743952702724802007 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dilog(1.5)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_dilog(2.0),  2.4674011002723397 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dilog(2.0)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_dilog( 5.0),  1.7837191612666306277 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dilog(5.0)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_dilog(12.59), 0.0010060918167266208634  ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_dilog(12.59)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_dilog(20.0), -1.2479770861745251168 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_dilog(20.0)");
  status += s;


  /* complex dilog */

  s = 0;
  gsl_sf_complex_dilog_impl(1.00001, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.20562022409960237363 ) > 1.0e-14 );
  s += ( frac_diff( y,  0.91597344814458309320 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_dilog(1.00001, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.99999, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.20561329262779687646 ) > 1.0e-14 );
  s += ( frac_diff( y,  0.91595774018131512060 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_dilog(0.99999, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.991, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.20250384721077806127 ) > 1.0e-06 );
  s += ( frac_diff( y,  0.90888544355846447810 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_dilog(0.991, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.98, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.19871638377785918403 ) > 1.0e-05 );
  s += ( frac_diff( y,  0.90020045882981847610 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_dilog(0.98, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.95, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.18848636456893572091 ) > 1.0e-14 );
  s += ( frac_diff( y,  0.87633754133420277830 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_dilog(0.95, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.8, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.13980800855429037810 ) > 1.0e-14 );
  s += ( frac_diff( y,  0.75310609092419884460 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_dilog(0.8, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.5, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.05897507442156586346 ) > 1.0e-14 );
  s += ( frac_diff( y,  0.48722235829452235710 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_dilog(0.5, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(0.01, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -0.000024999375027776215378 ) > 1.0e-12 );
  s += ( frac_diff( y,  0.009999888892888684820    ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_dilog(0.01, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(10.0, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -3.0596887943287347304 ) > 1.0e-14 );
  s += ( frac_diff( y,  3.7167814930680685900 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_dilog(10.0, Pi/2)");
  status += s;

  s = 0;
  gsl_sf_complex_dilog_impl(100.0, M_PI/2.0, &x, &y);
  s += ( frac_diff( x, -11.015004738293824854 ) > 1.0e-14 );
  s += ( frac_diff( y,  7.2437843013083534970 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_dilog(100.0, Pi/2)");
  status += s;

  return status;
}

int check_ellint(void)
{
  int status = 0;
  int s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_ellint_Kcomp( 0.99, 1.0e-3), 3.3566005233611923760 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_ellint_Kcomp( 0.50, 1.0e-3), 1.6857503548125960429 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_ellint_Kcomp(0.010, 1.0e-3), 1.5708355989121522360 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_ellint_Kcomp");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_ellint_Ecomp(0.99, 1.0e-3), 1.0284758090288040010 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_ellint_Ecomp(0.50, 1.0e-3), 1.4674622093394271555 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_ellint_Ecomp(0.01, 1.0e-3), 1.5707570561503852873 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_ellint_Ecomp");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_ellint_F(M_PI/3.0, 0.99, 1.0e-3), 1.3065333392738766762 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_ellint_F(M_PI/3.0, 0.50, 1.0e-3), 1.0895506700518854093 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_ellint_F(M_PI/3.0, 0.01, 1.0e-3), 1.0472129063770918952 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_ellint_F");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_ellint_E(M_PI/3.0, 0.99, 1.0e-3), 0.8704819220377943536 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_ellint_E(M_PI/3.0, 0.50, 1.0e-3), 1.0075555551444720293 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_ellint_E(M_PI/3.0, 0.01, 1.0e-3), 1.0471821963889481104 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_ellint_E");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_ellint_P(M_PI/3.0, 0.99, 0.5, 1.0e-3), 1.1288726598764099882 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_ellint_P(M_PI/3.0, 0.50, 0.5, 1.0e-3), 0.9570574331323584890 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_ellint_P(M_PI/3.0, 0.01, 0.5, 1.0e-3), 0.9228868127118118465 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_ellint_P");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_ellint_RF(5.0e-11, 1.0e-10, 1.0, 1.0e-3), 12.36441982 ) > 1.0e-09 );
  gsl_test(s, "  gsl_sf_ellint_RF(5.0e-11, 1.0e-10, 1.0, 1.0e-3)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_ellint_RF(1.0, 2.0, 3.0, 1.0e-3), 0.726946 ) > 1.0e-05 );
  gsl_test(s, "  gsl_sf_ellint_RF(1.0, 2.0, 3.0, 1.0e-3)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_ellint_RD(5.0e-11, 1.0e-10, 1.0, 1.0e-3), 34.09325948 ) > 1.0e-09 );
  gsl_test(s, "  gsl_sf_ellint_RD(5.0e-11, 1.0e-10, 1.0, 1.0e-3)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_ellint_RD(1.0, 2.0, 3.0, 1.0e-3), 0.29046 ) > 1.0e-04 );
  gsl_test(s, "  gsl_sf_ellint_RD(1.0, 2.0, 3.0, 1.0e-3)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_ellint_RC(1.0, 2.0, 1.0e-3), 0.7853981630 ) > 1.0e-09 );
  gsl_test(s, "  gsl_sf_ellint_RC(1.0, 2.0, 1.0e-3)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_ellint_RJ(2.0, 3.0, 4.0, 5.0, 1.0e-3), 0.1429757967 ) > 1.0e-09 );
  gsl_test(s, "  gsl_sf_ellint_RJ(2.0, 3.0, 4.0, 5.0, 1.0e-3)");
  status += s;

  return status;
}

int check_erf(void)
{
  int status = 0;
  int s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_erfc(-10.0), 2.0000000000000000000 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_erfc(-10.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_erfc(-1.0), 1.8427007929497148693 ) > 1.0e-9 );
  gsl_test(s, "  gsl_sf_erfc(-1.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_erfc(1.0), 0.15729920705028513066 ) > 1.0e-8 );
  gsl_test(s, "  gsl_sf_erfc(1.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_erfc(10.0), 2.0884875837625447570e-45 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_erfc(10.0)");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_log_erfc(-10.0), log(2.0000000000000000000) ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_log_erfc(-10.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_log_erfc(1.0), log(0.15729920705028513066) ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_log_erfc(1.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_log_erfc(10.0), log(2.0884875837625447570e-45) ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_log_erfc(10.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_erf(-10.0), -1.0000000000000000000 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_erf(-10.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_erf(1.0), 0.8427007929497148693 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_erf(1.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_erf(10.0), 1.0000000000000000000 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_erf(10.0)");
  status += s;

  return status;
}


int check_exp(void)
{
  int status = 0;
  int s;
 
  s = 0;
  s += ( frac_diff( gsl_sf_exp(-10.0), exp(-10.0) ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_exp( 10.0), exp( 10.0) ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_exp");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_expm1(-10.0),   exp(-10.0)-1.0  ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_expm1(-0.001),  -0.00099950016662500845 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_expm1(-1.0e-8), -1.0e-08 + 0.5e-16 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_expm1( 1.0e-8),  1.0e-08 + 0.5e-16 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_expm1( 0.001),   0.0010005001667083417  )> 1.0e-14 );
  s += ( frac_diff( gsl_sf_expm1( 10.0),   exp(10.0)-1.0  ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_expm1");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_exprel(-10.0),   (exp( -10.0)-1.0)/(-10.0)  ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_exprel(-0.001),   0.9995001666250084 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_exprel(-1.0e-8),   1.0 - 0.5e-08 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_exprel( 1.0e-8),   1.0 + 0.5e-08 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_exprel( 0.001),   1.0005001667083417 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_exprel( 10.0),   (exp(  10.0)-1.0)/(10.0)  ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_exprel");
  status += s;

  return status;
}


int check_gamma(void)
{
  double zr, zi, lg_r, lg_i, lg, sgn;
  int status = 0;
  int s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_lngamma(0.1),    2.252712651734205 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_lngamma(0.1)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_lngamma(100.0),  359.1342053695753 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_lngamma(100)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_lngamma(-0.1),    2.368961332728788655 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_lngamma(-0.1)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_lngamma(-1.000001), 13.815510135181261473 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_lngamma(-1.000001)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_lngamma(-1.0-1.0e-8), 18.420680739724522253 ) > 1.e-08 );
  gsl_test(s, "  gsl_sf_lngamma(-1.0-1.0e-8)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_lngamma(-1.0-1.0e-10), 23.025850929898178407 ) > 1.0e-06 );
  gsl_test(s, "  gsl_sf_lngamma(-1.0-1.0e-10)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_lngamma(-100.5), -364.9009683094273518 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_lngamma(-100.5)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(0.7, &lg, &sgn);
  s += ( frac_diff( lg,  0.26086724653166651439 ) > 1.e-14 );
  s += ( frac_diff( sgn, 1.0) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(0.7)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(0.1, &lg, &sgn);
  s += ( frac_diff( lg,  2.2527126517342059599 ) > 1.e-14 );
  s += ( frac_diff( sgn, 1.0) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(0.1)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-0.1, &lg, &sgn);
  s += ( frac_diff( lg,   2.368961332728788655 ) > 1.e-14 );
  s += ( frac_diff( sgn, -1.0) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-0.1)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-1.0-1.0e-6, &lg, &sgn);
  s += ( frac_diff( lg,  13.8155101351812614727 ) > 1.e-10 );
  s += ( frac_diff( sgn, 1.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-1.0-1.0e-6)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-2.0-1.0e-6, &lg, &sgn);
  s += ( frac_diff( lg,  13.1223624546214411633 ) > 1.e-10 );
  s += ( frac_diff( sgn, -1.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-2.0-1.0e-6)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-3.1, &lg, &sgn);
  s += ( frac_diff( lg,  0.4003116967039859208 ) > 1.e-14 );
  s += ( frac_diff( sgn, 1.0) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-3.1)");
  status += s;

  s = 0;
  gsl_sf_lngamma_sgn_impl(-100.5, &lg, &sgn);
  s += ( frac_diff( lg,  -364.9009683094273518   ) > 1.e-14 );
  s += ( frac_diff( sgn, -1.0) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_sgn_impl(-100.5)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_gamma(10.0),   362880.0 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_gamma(-10.5), -2.640121820547716316e-07 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_gamma(-11.2),  8.200835193555368332e-08 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_gamma(-1.000001),  999999.5772170767414 ) > 1.0e-08 );
  s += ( frac_diff( gsl_sf_gamma(-1.0-1.0e-10), 9.999999999577215665e+09 ) > 1.0e-06 );
  gsl_test(s, "  gsl_sf_gamma");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_gammainv(10.0),   1.0/362880.0 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_gammainv(-10.5), -1.0/2.640121820547716316e-07 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_gammainv(-11.2),  1.0/8.200835193555368332e-08 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_gammainv(-1.000001),  1.0/999999.5772170767414 ) > 1.0e-08 );
  s += ( frac_diff( gsl_sf_gammainv(-1.0-1.0e-10),  1.0/9.999999999577215665e+09 ) > 1.0e-06 );
  s += ( fabs(gsl_sf_gammainv(-100.0)) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_gammainv");
  status += s;

  s = 0;
  zr = 5.0;
  zi = 2.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, 2.7487017561338026749 ) > 1.0e-14 );
  s += ( frac_diff( lg_i, 3.0738434100497007915 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(5 + 2 I)");
  status += s;
  
  s = 0;
  zr = 100.0;
  zi = 100.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, 315.07804459949331323 ) > 1.0e-14 );
  s += ( frac_diff( lg_i, 2.0821801804113110099 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(100 + 100 I)");
  status += s;
  
  s = 0;
  zr =   100.0;
  zi = -1000.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, -882.3920483010362817000 ) > 1.0e-14 );
  s += ( frac_diff( lg_i,   -2.1169293725678813270 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(100 - 1000 I)");
  status += s;

  s = 0;
  zr = -100.0;
  zi =   -1.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, -365.0362469529239516000 ) > 1.0e-14 );
  s += ( frac_diff( lg_i,   -3.0393820262864361140 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(-1000 - I)");
  status += s;

  s = 0;
  s += ( frac_diff(  gsl_sf_fact(7), 5040.0 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_fact(33), 8.683317618811886496e+36 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_fact");
  status += s;

  s = 0;
  s += ( frac_diff(  gsl_sf_doublefact(7), 105.0 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_doublefact(33), 6.332659870762850625e+18 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_doublefact");
  status += s;

  s = 0;
  s += ( frac_diff(  gsl_sf_lnfact(7), 8.525161361065414300 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_lnfact(33), 85.05446701758151741 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lnfact");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_lnchoose(7,3), 3.555348061489413680 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_lnchoose(5,2), 2.302585092994045684 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lnchoose");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_choose(7,3), 35.0 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_choose(5,2), 10.0 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_choose(500,200), 5.054949849935532221e+144 ) > 1.0e-11 );
  gsl_test(s, "  gsl_sf_choose");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_lnpoch(7,3), 6.222576268071368616 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_lnpoch(5,2), 3.401197381662155375 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_lnpoch(5,0.01), 0.015072234709399379065 ) > 1.0e-10 );
  gsl_test(s, "  gsl_sf_lnpoch");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_poch(7,3), 504.0 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_poch(5,2),  30.0) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_poch(5,0.01), 1.0151863936613682753 ) > 1.0e-10 );
  gsl_test(s, "  gsl_sf_poch");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_pochrel(7,3), 503.0/3.0 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_pochrel(5,2),  29.0/2.0 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_pochrel(5,0.01),  1.5186393661368275330 ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_pochrel(-5.5,0.01),  1.8584945633829063516  ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_pochrel(-5.5,-0.01), 1.7292598184510393273  ) > 1.0e-14 );
  s += ( frac_diff( gsl_sf_pochrel(-5.5,-11.0), 0.09090909090939652475 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_pochrel");
  status += s;

  return status;
}

int check_gegen(void)
{
  double y;
  int status = 0;
  int s;

  s = 0;
  s += ( frac_diff( gsl_sf_gegenpoly_1(-0.2,   1.0),  -0.4 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_1( 0.0,   1.0),   2.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_1( 1.0,   1.0),   2.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_1( 1.0,   0.5),   1.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_1( 5.0,   1.0),  10.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_1( 100.0, 0.5), 100.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_gegenpoly_1");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_gegenpoly_2(-0.2,   0.5),   0.12 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_2( 0.0,   1.0),   1.00 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_2( 1.0,   1.0),   3.00 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_2( 1.0,   0.1),  -0.96 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_2( 5.0,   1.0),   55.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_2( 100.0, 0.5), 4950.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_gegenpoly_2");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_gegenpoly_3(-0.2,   0.5),      0.112 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_3( 0.0,   1.0),   -2.0/3.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_3( 1.0,   1.0),      4.000 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_3( 1.0,   0.1),     -0.392 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_3( 5.0,   1.0),    220.000 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_3( 100.0, 0.5), 161600.000 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_gegenpoly_3");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_gegenpoly_n(1,       1.0, 1.0),  2.000		    ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(10,      1.0, 1.0), 11.000		    ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(10,      1.0, 0.1), -0.4542309376	    ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(10,      5.0, 1.0),  9.23780e+4  	    ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(10,    100.0, 0.5),  1.5729338392690000e+13  ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(1000,  100.0, 1.0),  3.3353666135627322e+232 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(100,  2000.0, 1.0),  5.8753432034937579e+202 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(103,   207.0, 2.0),  1.4210272202235983e+145 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(103,    -0.4, 0.3), -1.64527498094522e-04    ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_gegenpoly_n");
  status += s;
  
  s += 0;
  s += (gsl_sf_gegenpoly_n_impl(103, -0.5,  1.0, &y) != GSL_EDOM);
  s += (gsl_sf_gegenpoly_n_impl(103, -0.51, 1.0, &y) != GSL_EDOM);
  gsl_test(s, "  gsl_sf_gegenpoly_n_impl: trap lambda <= -1/2");
  status += s;

  s += 0;
  s += (gsl_sf_gegenpoly_n_impl(-2, 0.0,  1.0, &y) != GSL_EDOM);
  gsl_test(s, "  gsl_sf_gegenpoly_n_impl: trap n < 0");
  status += s;
  
  return status;
}


int check_hyperg(void)
{
  double y;
  int status = 0;
  int s;


  /* 0F1 */

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_0F1(1, 0.5),      1.5660829297563505373 ) > 1.e-11 );
  s += ( frac_diff(gsl_sf_hyperg_0F1(5, 0.5),      1.1042674404828684574 ) > 1.e-11 );
  s += ( frac_diff(gsl_sf_hyperg_0F1(100, 30),     1.3492598639485110176 ) > 1.e-11 );
  s += ( frac_diff(gsl_sf_hyperg_0F1(-0.5, 3),    -39.29137997543434276  ) > 1.e-11 );
  s += ( frac_diff(gsl_sf_hyperg_0F1(-100.5, 50),  0.6087930289227538496 ) > 1.e-11 );
  s += ( frac_diff(gsl_sf_hyperg_0F1(1, -5.0),    -0.3268752818235339109 ) > 1.e-11 );
  s += ( frac_diff(gsl_sf_hyperg_0F1(-0.5, -5.0), -4.581634759005381184  ) > 1.e-11 );
  gsl_test(s, "  gsl_sf_hyperg_0F1");
  status += s;


  /* 1F1 */

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_1F1(1, 1, 0.5), 1.6487212707001281468 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(1, 1, 0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_1F1(8, 1, 0.5), 13.108875178030540372 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(8, 1, 0.5)");

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_1F1(8, 1, 8), 5.481422453671217135e+7 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(8, 1, 8)");

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_1F1(1, 8, 8), 4.918996932358889820 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(1, 8, 8)");
  status += s;

  s = 0;
  gsl_sf_hyperg_1F1_impl(10, 2.5, -2.0, &y);
  s += ( frac_diff(y , 0.008423875120475747307 ) > 1.e-11 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(10, 2.5, -2.0)");
  status += s;

  s = 0;
  gsl_sf_hyperg_1F1_impl(10, 2.5, -8.0, &y);
  s += ( frac_diff(y, 0.00019157621858903819077 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(10, 2.5, -8.0)");
  status += s;

  s = 0;
  gsl_sf_hyperg_1F1_impl(100, 2.5, -10.0, &y);
  s += ( frac_diff(y, -5.011573844831150295e-06 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_1F1(100, 2.5, -10.0)");
  status += s;

  
  /* U */

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 0.0001, 0.0001), 1.0000576350699863577 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 0.0001, 0.0001)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 0.0001, 1.0), 0.9999403679233247536 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 0.0001, 1.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 0.0001, 100.0), 0.9995385992657260887 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 0.0001, 100.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 1, 0.0001), 1.0009210608660065989 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 1, 0.0001)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 1.0, 1.0), 0.9999999925484179084 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 1.0, 1.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 10, 1), 13.567851006281412726 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 10, 1)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 10, 5), 1.0006265020064596364 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 10, 5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 10, 10), 0.9999244381454633265 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 10, 10)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 100, 1), 2.5890615708804247881e+150 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 100, 1)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 100, 10), 2.3127845417739661466e+55 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 100, 10)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 100, 50), 6402.818715083582554 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 100, 50)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 100, 98), 0.9998517867411840044 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 100, 98)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 1000, 300), 2.5389557274938010716e+213 ) > 1.e-06 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 1000, 300)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 1000, 998), 0.9997236131433337327 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 1000, 998)");
  printf("%22.18g\n", gsl_sf_hyperg_U(0.0001, 1000, 998));
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.0001, 1000, 999), 0.9997195294193261604 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.0001, 1000, 999)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.5, 1000, 300), 1.1977955438214207486e+217 ) > 1.e-08 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.5, 1000, 300)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.5, 1000, 800), 9.103916020464797207e+08 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.5, 1000, 800)");
  printf("%22.18g\n", gsl_sf_hyperg_U(0.5, 1000, 800));
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.5, 1000, 998), 0.21970269691801966806 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.5, 1000, 998)");
  status += s;



  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(0.5, 0.5, 1.0), 0.7578721561413121060 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(0.5, 0.5, 1.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 0.0001, 0.0001), 0.9992361337764090785 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 0.0001, 0.0001)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 0.0001, 1), 0.4036664068111504538 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 0.0001, 1)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 0.0001, 100), 0.009805780851264329587 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 0.0001, 100)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 1, 0.0001), 8.634088070212725330 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 1, 0.0001)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 1, 0.01), 4.078511443456425847 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 1, 0.01)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 1, 0.5), 0.9229106324837304688 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 1, 0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 1, 2.0), 0.3613286168882225847 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 1, 2.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 1, 100), 0.009901942286733018406 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 1, 100)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 1, 1000), 0.0009990019940238807150 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 1, 1000)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 1.2, 2.0), 0.3835044780075602550 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 1.2, 2.0)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 8, 1), 1957.0000000000000000 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 8, 1)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 8, 5), 1.0424960000000000000 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 8, 5)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 8, 8), 0.3207168579101562500 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 8, 8)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 8, 50), 0.022660399001600000000 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 8, 50)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 8, 100), 0.010631236727200000000 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 8, 100)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 8, 1000), 0.0010060301203607207200 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 8, 1000)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 20, 1), 1.7403456103284421000e+16 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 20, 1)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 20, 20), 0.22597813610531052969 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 20, 20)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 50, 1), 3.374452117521520758e+61 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 50, 1)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 50, 50), 0.15394136814987651785 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 50, 50)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 100, 1), 2.5624945006073464385e+154 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 100, 1)");
  status += s;

   s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 100, 50), 3.0978624160896431391e+07 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 100, 50)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 100, 100), 0.11323192555773717475 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 100, 100)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 100, 200), 0.009715680951406713589 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 100, 200)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 100, 1000), 0.0011085142546061528661 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 100, 1000)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, 1000, 2000), 0.0009970168547036318206 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, 1000, 2000)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -0.0001, 1), 0.4036388693605999482 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -0.0001, 1)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -1, 1), 0.29817368116159703717 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -1, 1)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -1, 10), 0.07816669698940409380 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -1, 10)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -10, 1), 0.08271753756946041959 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -10, 1)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -10, 5), 0.06127757419425055261 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -10, 5)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -10, 10), 0.04656199948873187212 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -10, 10)");
  printf("%22.18g\n", gsl_sf_hyperg_U(1, -10, 10));
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -10, 20), 0.031606421847946077709 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -10, 20)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -100, 1), 0.009802970197050404429 ) > 1.e-08 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -100, 1)");
  printf("%22.18g\n", gsl_sf_hyperg_U(1, -100, 1));
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -100, 10), 0.009001648897173103447 ) > 1.e-06 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -100, 10)");
  printf("%22.18g\n", gsl_sf_hyperg_U(1, -100, 10));
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -100, 20), 0.008253126487166557546 ) > 1.e-06 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -100, 20)");
  printf("%22.18g\n", gsl_sf_hyperg_U(1, -100, 20));
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -100, 50), 0.006607993916432051008 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -100, 50)");
  printf("%22.18g\n", gsl_sf_hyperg_U(1, -100, 50));
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -100, 90), 0.005222713769726871937 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -100, 90)");
  printf("%22.18g\n", gsl_sf_hyperg_U(1, -100, 90));
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -100, 110), 0.004727658137692606210 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -100, 110)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -1000, 1), 0.0009980029970019970050 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -1000, 1)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(1, -1000, 1010), 0.0004971408839859245170 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_U(1, -1000, 1010)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(8, 1, 0.5), 6.449509938973479986e-06 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_U(8, 1, 0.5)");
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(8, 1, 8), 6.190694573035761284e-10 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(8, 1, 8)");
  printf("%22.18g\n", gsl_sf_hyperg_U(8, 1, 8));
  status += s;
  
  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(8, 8, 1), 0.12289755012652317578 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(8, 8, 1)");
  printf("%22.18g\n", gsl_sf_hyperg_U(8, 8, 1));
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(8, 8, 10), 5.687710359507564272e-09 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(8, 8, 10)");
  printf("%22.18g\n", gsl_sf_hyperg_U(8, 8, 10));
  status += s;

  s = 0;
  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(8, 10.5, 1), 27.981926466707438538 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(8, 10.5, 1)");
  printf("%22.18g\n", gsl_sf_hyperg_U(8, 10.5, 1));
  status += s;

  status += s;
  s += ( frac_diff(gsl_sf_hyperg_U(100, 100, 1), 0.009998990209084729106 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(100, 100, 1)");
  printf("%22.18g\n", gsl_sf_hyperg_U(100, 100, 1));
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(10, -2.5, 10), 6.734690720346560349e-14 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(10, -2.5, 10)");
  printf("%22.18g\n", gsl_sf_hyperg_U(10, -2.5, 10));
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_U(10, 2.5, 50), 2.4098720076596087125e-18 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_U(10, 2.5, 50)");
  printf("%22.18g\n", gsl_sf_hyperg_U(10, 2.5, 50));
  status += s;




  /* 2F1 */

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(1, 1, 1, 0.5), 2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(1,1,1,0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(8, 8, 1, 0.5), 12451584.0 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(8, 8, 1, 0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(8, -8, 1, 0.5), 0.13671875 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(8, -8, 1, 0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(8, -8.1, 1, 0.5), 0.14147385378899930422 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(8, -8.1, 1, 0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(8, -8, 1, -0.5), 4945.136718750000000 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(8, -8, 1, -0.5)");
  status += s;

  s = 0;
  gsl_sf_hyperg_2F1_impl(8, -8, -5.5, 0.5, &y);
  s += ( frac_diff(y, -906.6363636363636364 ) > 1.e-11 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(8, -8, -5.5, 0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(8, -8, -5.5, -0.5), 24565.363636363636364 ) > 1.e-12 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(8, -8, -5.5, -0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(8, 8, 1, -0.5), -0.006476312098196747669 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(8, 8, 1, -0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(8, 8, 5, 0.5), 4205.714285714285714 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(8, 8, 5, 0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(8, 8, 5, -0.5), 0.0028489656290296436616 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(8, 8, 5, -0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(9, 9, 1, 0.99),  1.2363536673577259280e+38  ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(9, 9, 1, 0.99)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(9, 9, -1.5, 0.99), 3.796186436458346579e+46 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(9, 9, -1.5, 0.99)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(9, 9, -1.5, -0.99), 0.14733409946001025146 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(9, 9, -1.5, -0.99)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(9, 9, -8.5, 0.99), -1.1301780432998743440e+65 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(9, 9, -8.5, 0.99)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(9, 9, -8.5, -0.99), -8.856462606575344483 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(9, 9, -8.5, -0.99)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(9, 9, -21.5, 0.99), 2.0712920991876073253e+95 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(9, 9, -21.5, 0.99)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(9, 9, -21.5, -0.99), -74.30517015382249216 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(9, 9, -21.5, -0.99)");
  status += s;

  s = 0;
  gsl_sf_hyperg_2F1_impl(9, 9, -100.5, 0.99, &y);
  s += ( frac_diff(y, -3.186778061428268980e+262 ) > 1.e-09 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(9, 9, -100.5, 0.99)");
  status += s;

  s = 0;
  gsl_sf_hyperg_2F1_impl(9, 9, -100.5, -0.99, &y);
  s += ( frac_diff(y, 2.4454358338375677520 ) > 1.e-09 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(9, 9, -100.5, -0.99)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(25, 25, 1, -0.5), -2.9995530823639545027e-06 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1(25, 25, 1, -0.5)");
  printf("%22.18g\n", gsl_sf_hyperg_2F1(25, 25, 1, -0.5));
  status += s;


  /* 2F1 conj */

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1_conj(1, 1, 1, 0.5), 3.352857095662929028 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1_conj(1,1,1,0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1_conj(8, 8, 1, 0.5), 1.7078067538891293983e+09 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1_conj(8, 8, 1, 0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1_conj(8, 8, 5, 0.5), 285767.15696901140627 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1_conj(8, 8, 5, 0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1_conj(8, 8, 1, -0.5), 0.007248196261471276276 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1_conj(8, 8, 1, -0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1_conj(8, 8, 5, -0.5), 0.00023301916814505902809 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1_conj(8, 8, 5, -0.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1_conj(25, 25, 1, -0.5), 5.169694409566320627e-06 ) > 1.e-10 );
  gsl_test(s, "  gsl_sf_hyperg_2F1_conj(25, 25, 1, -0.5)");
  printf("%22.18g\n", gsl_sf_hyperg_2F1_conj(25, 25, 1, -0.5));
  status += s;


  /* FIXME: the "true" values here may not be so good */
  /*
  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F0(0.01, 1.0, -0.02), 0.999803886708565    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_2F0(0.1,  0.5, -0.02), 0.999015947934831    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_2F0(1,   1, -0.02),   0.980755496569062     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_2F0(8,   8, -0.02),   0.3299059284994299    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_2F0(50, 50, -0.02),   2.688995263773233e-13 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_2F0");
  status += s;
  */

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
  s += ( frac_diff( sn, 0.4707504736556572833 ) > 1.0e-14 );
  s += ( frac_diff( cn, 0.8822663948904402865 ) > 1.0e-14 );
  s += ( frac_diff( dn, 0.9429724257773856873 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_elljac(0.5|0.5)");
  status += s;

  u = 2.0;
  m = 0.999999;
  s = 0;
  stat_ej = gsl_sf_elljac_impl(u, m, &sn, &cn, &dn);
  s += ( frac_diff( sn, 0.96402778575700186570 ) > 1.0e-14 );
  s += ( frac_diff( cn, 0.26580148285600686381 ) > 1.0e-14 );
  s += ( frac_diff( dn, 0.26580323105264131136 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_elljac(2.0|0.999999)");
  status += s;

  return status;
}

int check_log(void)
{
  double x, y;
  int status = 0;
  int s;
  
  s = 0;
  gsl_sf_complex_log_impl(1.0, 1.0, &x, &y);
  s += ( frac_diff( x, 0.3465735902799726547 ) > 1.0e-14 );
  s += ( frac_diff( y, 0.7853981633974483096 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_log(1 + I)");
  status += s;
  
  s = 0;
  gsl_sf_complex_log_impl(1.0, -1.0, &x, &y);
  s += ( frac_diff( x,  0.3465735902799726547 ) > 1.0e-14 );
  s += ( frac_diff( y, -0.7853981633974483096 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_log(1 - I)");
  status += s;
  
  s = 0;
  gsl_sf_complex_log_impl(1.0, 100.0, &x, &y);
  s += ( frac_diff( x, 4.605220183488258022 ) > 1.0e-14 );
  s += ( frac_diff( y, 1.560796660108231381 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_log(1 + 100 I)");
  status += s;

  s = 0;
  gsl_sf_complex_log_impl(-1000.0, -1.0, &x, &y);
  s += ( frac_diff( x,  6.907755778981887052  ) > 1.0e-14 );
  s += ( frac_diff( y, -3.1405926539231263718 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_log(-1000 - I)");
  status += s;

  s = 0;
  gsl_sf_complex_log_impl(-1.0, 0.0, &x, &y);
  s += ( frac_diff( y, 3.1415926535897932385 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_log(-1)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_log_1plusx(0.1), 0.09531017980432486004 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_log_1plusx(0.1)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_log_1plusx(0.49), 0.3987761199573677730 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_log_1plusx(0.49)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_log_1plusx(-0.49), -0.6733445532637655964 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_log_1plusx(-0.49)");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_log_1plusx(1.0e-10), 9.999999999500000000e-11 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_log_1plusx(1.0e-10)");
  status += s;

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
  s += ( frac_diff(y, 1 + 0.5*x + 0.3*x*x) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_poly_eval({1, 0.5, 0.3}, 0.5)");
  status += s;
  
  s = 0;
  x = 1.0;
  y = gsl_sf_poly_eval(d, 11, x);
  s += ( frac_diff(y, 1.0) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_poly_eval({1,-1, 1, -1, 1, -1, 1, -1, 1, -1, 1}, 1.0)");
  status += s;

  return status;
}

int check_pow_int(void)
{
  int status = 0;
  int s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_pow_int(2.0, 3), 8.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_pow_int(2.0, 3)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_pow_int(-2.0, 3), -8.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_pow_int(-2.0, 3)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_pow_int(2.0, -3), 1.0/8.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_pow_int(2.0, -3)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_pow_int(-2.0, -3), -1.0/8.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_pow_int(-2.0, -3)");
  status += s;
  
  return status;
}

int check_psi(void)
{
  int status = 0;
  int s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_psi_int(5), 1.5061176684318004727 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi_int(5)");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_psi_int(5000), 8.517093188082904107 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi_int(5000)");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_psi(5.0), 1.5061176684318004727 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi(5.0)");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_psi(5000.0), 8.517093188082904107 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi(5000.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_psi(-100.5), 4.615124601338064117 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi(-100.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_psi(-1.0e+5-0.5), 11.512935464924395337 ) > 1.0e-10 );
  gsl_test(s, "  gsl_sf_psi(-1.0e+5-0.5)");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_psi_1piy(0.8), -0.07088340212750589223 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi_1piy(0.8)");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_psi_1piy(1.0), 0.09465032062247697727 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi_1piy(1.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_psi_1piy(5.0), 1.6127848446157465854  ) > 1.0e-13 );
  gsl_test(s, "  gsl_sf_psi_1piy(5.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_psi_1piy(100.0),  4.605178519404762003 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi_1piy(100.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_psi_1piy(2000.0), 7.600902480375416216 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi_1piy(2000.0)");
  status += s;

  return status;
}


int check_synch(void)
{
  double y;
  int status = 0;
  int s;
  
  s = 0;
  gsl_sf_synchrotron_1_impl(0.01, &y);
  s += ( frac_diff( y, 0.444973 ) > 1.0e-05 );
  gsl_test(s, "  gsl_sf_synchrotron_1(0.01)");
  status += s;

  s = 0;
  gsl_sf_synchrotron_1_impl(1.0, &y);
  s += ( frac_diff( y, 0.651423 ) > 1.0e-05 );
  gsl_test(s, "  gsl_sf_synchrotron_1(1.0)");
  status += s;

  s = 0;
  gsl_sf_synchrotron_1_impl(10.0, &y);
  s += ( frac_diff( y, 0.000192238 ) > 1.0e-05 );
  gsl_test(s, "  gsl_sf_synchrotron_1(10.0)");
  status += s;

  s = 0;
  gsl_sf_synchrotron_1_impl(100.0, &y);
  s += ( frac_diff( y, 4.69759e-43 ) > 1.0e-05 );
  gsl_test(s, "  gsl_sf_synchrotron_1(100.0)");
  status += s;

  s = 0;
  gsl_sf_synchrotron_2_impl(0.01, &y);
  s += ( frac_diff( y, 0.23098077342226277732 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_synchrotron_2(0.01)");
  status += s;

  s = 0;
  gsl_sf_synchrotron_2_impl(1.0, &y);
  s += ( frac_diff( y, 0.4944750621042082670 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_synchrotron_2(1.0)");
  status += s;

  s = 0;
  gsl_sf_synchrotron_2_impl(10.0, &y);
  s += ( frac_diff( y, 0.00018161187569530204281 ) > 1.0e-05 );
  gsl_test(s, "  gsl_sf_synchrotron_2(10.0)");
  status += s;
  
  s = 0;
  gsl_sf_synchrotron_2_impl(100.0, &y);
  s += ( frac_diff( y, 4.666936458728046656e-43 ) > 1.0e-05 );
  gsl_test(s, "  gsl_sf_synchrotron_2(100.0)");
  status += s;

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

  yr = 1.0;
  yi = 5.0;
  gsl_sf_complex_sin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr, 62.44551846769653403 ) > 1.0e-14 );
  s += ( frac_diff( zi, 40.09216577799840254 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_sin_impl(1 + 5 I)");
  status += s;
  
  yr = 1.0;
  yi = 5.0;
  gsl_sf_complex_cos_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr,  40.09580630629882573 ) > 1.0e-14 );
  s += ( frac_diff( zi, -62.43984868079963017 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_cos_impl(1 + 5 I)");
  status += s;

  yr =   1.0;
  yi = 100.0;
  gsl_sf_complex_logsin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr, 99.3068528194400546900 ) > 1.0e-14 );
  s += ( frac_diff( zi,  0.5707963267948966192 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_logsin_impl(1 + 100 I)");
  status += s;
  
  yr =    1.0;
  yi = -100.0;
  gsl_sf_complex_logsin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr,  99.3068528194400546900 ) > 1.0e-14 );
  s += ( frac_diff( zi,  -0.5707963267948966192 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_logsin_impl(1 - 100 I)");
  status += s;

  yr = 5.0;
  yi = 5.0;
  gsl_sf_complex_logsin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr, 4.3068909128079757420 ) > 1.0e-14 );
  s += ( frac_diff( zi, 2.8540063315538773952 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_logsin_impl(5 + 5 I)");
  status += s;

  gsl_sf_polar_to_rect_impl(10.0, M_PI/6.0, &x, &y);
  s = 0;
  s += ( frac_diff( x, 10.0 * sqrt(3) / 2.0 ) > 1.0e-14 );
  s += ( frac_diff( y, 10.0 * 0.5           ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_polar_to_rect_impl(10, Pi/6)");
  status += s;
  
  gsl_sf_polar_to_rect_impl(10.0, -2.0/3.0*M_PI, &x, &y);
  s = 0;
  s += ( frac_diff( x, 10.0 * (-0.5)           ) > 1.0e-14 );
  s += ( frac_diff( y, 10.0 * (-sqrt(3) / 2.0) ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_polar_to_rect_impl(10, -2/3 Pi)");
  status += s;
  
  theta = 5.0*M_PI + M_PI/2.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, 3.0/2.0*M_PI ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_pos_impl: theta =  11/2 Pi");
  status += s;

  theta = -5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_pos_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_pos_impl: theta = -11/2 Pi");
  status += s;

  theta = 5.0*M_PI + M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, -M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta =  11/2 Pi");
  status += s;

  theta = -5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta = -11/2 Pi");
  status += s;

  theta =  5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta = -9/2 Pi");
  status += s;

  theta =  3.0/2.0*M_PI;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, -M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta =  3/2 Pi");
  status += s;

  theta = -3.0/2.0*M_PI;
  gsl_sf_angle_restrict_symm_impl(&theta);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta = -3/2 Pi");
  status += s;

  return status;
}


int check_zeta(void)
{
  int status = 0;
  int s;

  s = 0;
  s += ( frac_diff(gsl_sf_zeta_int(-61), -3.30660898765775767257e+34) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta_int(-5),  -0.003968253968253968253968) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta_int(5),   1.0369277551433699263313655) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta_int(31),  1.0000000004656629065033784) > 1.e-14 );
  gsl_test(s, "  gsl_sf_zeta_int");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_zeta(-51),    9.68995788746359406565e+24) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta(-5),    -0.003968253968253968253968) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta(-0.5),  -0.207886224977354566017307) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta(0.5),   -1.460354508809586812889499) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta(0.999), -999.4228571557887900      ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_zeta(1.0+1.0e-6), 1.0000005772157377174e+06      ) > 1.e-10 );
  s += ( frac_diff(gsl_sf_zeta(5),      1.036927755143369926331365) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta(25.5),   1.000000021074106110269959) > 1.e-14 );
  gsl_test(s, "  gsl_sf_zeta");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hzeta(2,  1.0),  1.6449340668482264365   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hzeta(5,  1.0),  1.0369277551433699263   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hzeta(9,  0.1),  1.0000000004253980e+9   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hzeta(2, 10.0),  0.1051663356816857461   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hzeta(5, 10.0),  0.0000304137986764702764) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hzeta");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_eta_int(-51), -4.363969073121683116e+40 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta_int(-5),  0.25 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta_int( 5),  0.9721197704469093059 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta_int( 6),  0.9855510912974351041 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta_int( 20), 0.9999990466115815221 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_eta_int");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_eta(-51.5), -1.2524184036924703656e+41 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta(-5),    0.25 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta(0.5),   0.6048986434216303702 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta(0.999), 0.6929872789683383574 ) > 1.0e-10 );
  s += ( frac_diff(gsl_sf_eta(1.0),   0.6931471805599453094 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta(1.0+1.0e-10), 0.6931471805759321998 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta( 5),    0.9721197704469093059 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta( 5.2),  0.9755278712546684682 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta( 6),    0.9855510912974351041 ) > 1.0e-14 );
  s += ( frac_diff(gsl_sf_eta( 20),   0.9999990466115815221 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_eta");
  status += s;

  return status;
}


int main(int argc, char * argv[])
{
/* test_coulomb(); */
/* test_recurse(); */

  gsl_test(check_airy(),       "Airy Functions");
  gsl_test(check_bessel(),     "Bessel Functions");
  gsl_test(check_cheb(),       "Chebyshev Evaluation");
  gsl_test(check_clausen(),    "Clausen Integral");
  gsl_test(check_coulomb(),    "Coulomb Wave Functions");
  gsl_test(check_coupling(),   "Coupling Coefficients");

  gsl_test(check_dawson(),     "Dawson Integral");
  gsl_test(check_debye(),      "Debye Functions");
  gsl_test(check_dilog(),      "Dilogarithm");

  gsl_test(check_ellint(),     "Elliptic Integrals");
  gsl_test(check_jac(),        "Jacobi Elliptic Functions");
  gsl_test(check_erf(),        "Error Functions");
  gsl_test(check_exp(),        "Exponential Functions");

  gsl_test(check_gamma(),      "Gamma Functions");
  gsl_test(check_gegen(),      "Gegenbauer Polynomials");
  gsl_test(check_hyperg(),     "Hypergeometric Functions");

  gsl_test(check_log(),        "Logarithm");
  gsl_test(check_poly(),       "Polynomial Evaluation");
  gsl_test(check_pow_int(),    "Integer Powers");
  gsl_test(check_psi(),        "Psi Functions");
  gsl_test(check_synch(),      "Synchrotron Functions");

  gsl_test(check_trig(),       "Trigonometric and Related Functions");
  gsl_test(check_zeta(),       "Zeta Functions");

  return gsl_test_summary();
}
