/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Calculate series for Y_nu and K_nu for small x and nu.
 * This is applicable for x < 2 and |nu|<=1/2.
 * These functions assume x > 0.
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "bessel_temme.h"


/* nu = (x+1)/4, -1<x<1, 1/(2nu)(1/Gamma[1-nu]-1/Gamma[1+nu]) */
static double g1_dat[14] = {
  -1.14516408366268311786898152867,
   0.0063608531134708423812295549463,
   0.00186245193007206848934643657413,
   0.000152833085873453507081227823964,
   0.0000170174640118020387953247317024,
  -6.4597502923347254354668326451e-7,
  -5.1819848432519380894104312968e-8,
   4.5189092894858183051123180797e-10,
   3.2433227371020873043666259180e-11,
   6.8309434024947522875432400828e-13,
   2.83535027551721015131196281296e-14,
  -7.9883905769323592875638087541e-16,
  -3.3726677300771949833341213457e-17,
  -3.6586334809210520744054437104e-20
};
static struct gsl_sf_cheb_series g1_cs = {
  g1_dat,
  13,
  -1, 1,
  (double *)0,
  (double *)0
};

/* nu = (x+1)/4, -1<x<1,  1/2 (1/Gamma[1-nu]+1/Gamma[1+nu]) */
static double g2_dat[15] = 
{
  1.88264552494967183501961697535,
 -0.077490658396167518329547945212,  
 -0.0182567148473249294195793409497,
  0.00063380302090748957959239717308,
  0.000076229054350872902119446117502,
 -9.5501647561720443519853993526e-7,
 -8.8927268107886351912431512955e-8,
 -1.95213347723196137405118801319e-9,
 -9.4003052735885162111769579771e-11,
  4.6875133849532393179290879101e-12,
  2.26585357469257595824475451450e-13,
 -1.17255096984880151118787352508e-15,
 -7.0441338200245222530843155877e-17,
 -2.43778783101076936506597402276e-18,
 -7.5225243218253901727164675011e-20
};
static struct gsl_sf_cheb_series g2_cs = {
  g2_dat,
  14,
  -1, 1,
  (double *)0,
  (double *)0
};


static
int
gsl_sf_temme_gamma(const double nu, double * g_1pnu, double * g_1mnu, double * g1, double * g2)
{
  const double anu = fabs(nu);    /* functions are even */
  const double x = 4.0*anu - 1.0;
  *g1 = gsl_sf_cheb_eval(&g1_cs, x);
  *g2 = gsl_sf_cheb_eval(&g2_cs, x);
  *g_1mnu = 1.0/(*g2 + nu * *g1);
  *g_1pnu = 1.0/(*g2 - nu * *g1);
  return GSL_SUCCESS;
}


int
gsl_sf_bessel_Y_temme(const double nu, const double x,
                      double * Y_nu, double * Y_nup1, double * Yp_nu)
{
  const int max_iter = 15000;
  
  const double half_x = 0.5 * x;
  const double ln_half_x = log(half_x);
  const double half_x_nu = exp(nu*ln_half_x);
  const double pi_nu   = M_PI * nu;
  const double alpha   = pi_nu / 2.0;
  const double sigma   = -nu * ln_half_x;
  const double sinrat  = (fabs(pi_nu) < GSL_MACH_EPS ? 1.0 : pi_nu/sin(pi_nu));
  const double sinhrat = (fabs(sigma) < GSL_MACH_EPS ? 1.0 : sinh(sigma)/sigma);
  const double sinhalf = (fabs(alpha) < GSL_MACH_EPS ? 1.0 : sin(alpha)/alpha);
  const double sin_sqr = nu*M_PI*M_PI*0.5 * sinhalf*sinhalf;
  
  double sum0, sum1;
  double fk, pk, qk, hk, ck;
  int k = 0;

  double g_1pnu, g_1mnu, g1, g2;
  gsl_sf_temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2);
  
  fk = 2.0/M_PI * sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
  pk = 1.0/M_PI /half_x_nu * g_1pnu;
  qk = 1.0/M_PI *half_x_nu * g_1mnu;
  hk = pk;
  ck = 1.0;

  sum0 = fk + sin_sqr * qk;
  sum1 = pk;

  while(k < max_iter) {
    double del0;
    double del1;
    double gk;
    k++;
    fk  = (k*fk + pk + qk)/(k*k-nu*nu);
    ck *= -half_x*half_x/k;
    pk /= (k - nu);
    qk /= (k + nu);
    gk  = fk + sin_sqr * qk;
    hk  = -k*gk + pk; 
    del0 = ck * gk;
    del1 = ck * hk;
    sum0 += del0;
    sum1 += del1;
    if (fabs(del0) < (1.0 + fabs(sum0))*GSL_MACH_EPS) break;
  }

  *Y_nu   = -sum0;
  *Y_nup1 = -sum1 * 2.0/x;
  *Yp_nu  = nu/x * *Y_nu - *Y_nup1;

  if(k == max_iter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}


int
gsl_sf_bessel_K_scaled_temme(const double nu, const double x,
                             double * K_nu, double * K_nup1, double * Kp_nu)
{
  const int max_iter = 15000;

  const double half_x    = 0.5 * x;
  const double ln_half_x = log(half_x);
  const double half_x_nu = exp(nu*ln_half_x);
  const double pi_nu   = M_PI * nu;
  const double sigma   = -nu * ln_half_x;
  const double sinrat  = (fabs(pi_nu) < GSL_MACH_EPS ? 1.0 : pi_nu/sin(pi_nu));
  const double sinhrat = (fabs(sigma) < GSL_MACH_EPS ? 1.0 : sinh(sigma)/sigma);
  const double ex = exp(x);

  double sum0, sum1;
  double fk, pk, qk, hk, ck;
  int k = 0;

  double g_1pnu, g_1mnu, g1, g2;
  gsl_sf_temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2);

  fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
  pk = 0.5/half_x_nu * g_1pnu;
  qk = 0.5*half_x_nu * g_1mnu;
  hk = pk;
  ck = 1.0;
  sum0 = fk;
  sum1 = hk;
  while(k < max_iter) {
    double del0;
    double del1;
    k++;
    fk  = (k*fk + pk + qk)/(k*k-nu*nu);
    ck *= half_x*half_x/k;
    pk /= (k - nu);
    qk /= (k + nu);
    hk  = -k*fk + pk;
    del0 = ck * fk;
    del1 = ck * hk;
    sum0 += del0;
    sum1 += del1;
    if(fabs(del0) < fabs(sum0)*GSL_MACH_EPS) break;
  }
  
  *K_nu   = sum0 * ex;
  *K_nup1 = sum1 * 2.0/x * ex;
  *Kp_nu  = - *K_nup1 + nu/x * *K_nu;
  if(k == max_iter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}
