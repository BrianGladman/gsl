/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#define locMin(a,b)  ((a) < (b) ? (a) : (b))
#define locMax(a,b)  ((a) > (b) ? (a) : (b))

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
temme_gamma(const double nu, double * g_1pnu, double * g_1mnu, double * g1, double * g2)
{
  const double anu = fabs(nu);    /* functions are even */
  const double x = 4.0*anu - 1.0;
  *g1 = gsl_sf_cheb_eval(&g1_cs, x);
  *g2 = gsl_sf_cheb_eval(&g2_cs, x);
  *g_1mnu = 1.0/(*g2 + nu * *g1);
  *g_1pnu = 1.0/(*g2 - nu * *g1);
  return GSL_SUCCESS;
}


/* Temme's series for K_nu and K_{nu+1}.
 * Assumes 0 < x < 2  and |nu| < 1/2.
 */
static
int
bessel_K_temme(const double nu, const double x, double * K_nu, double * K_nup1)
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
  temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2);

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
  if(k == max_iter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}


/* Evaluate continued fraction CF1 to obtain I'/I.
 */
static
int
bessel_I_CF1(const double nu, const double x, double * result)
{
  const int max_iter = 5000;
  int i = 0;
  double x_inv = 1.0/x;
  double r = nu * x_inv;
  double b = 2.0*nu*x_inv;
  double d = 0.0;
  double c;
  r = locMax(r, 1.0e-100);
  c = r;
  while(i < max_iter) {
    double del;
    b  += 2.0*x_inv;;
    d   = 1.0/(b+d);
    c   = b + 1.0/c;
    del = c*d;
    r  *= del;
    if (fabs(del-1.0) < GSL_MACH_EPS) break;
    ++i;
  }
  
  *result = r;
  if(i == max_iter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}


/* Downward recurrence to obtain I,I' values:
 * I_{nu_min + kmax}, I_{nu_min + kmax - 1}, ..., I_{nu_min}
 * and similarly for I'
 */
static
int
bessel_I_recur(const double nu_min, const double x, const int kmax,
               const double I_start, const double Ip_start,
	       double * I_end, double * Ip_end,
	       double * I_array, double * Ip_array
	       )
{
  double x_inv  = 1.0/x;
  double nu_max = nu_min + kmax;
  double I_nu  = I_start;
  double Ip_nu = Ip_start;
  double nu = nu_max;
  int k;
  if(I_array  != (double *)0)  I_array[kmax] = I_start;
  if(Ip_array != (double *)0) Ip_array[kmax] = Ip_start;
  for(k=kmax-1; k>=0; k--) {
    double nu_x_inv = nu*x_inv;
    double I_nu_save = I_nu;
    I_nu  = nu_x_inv*I_nu + Ip_nu;
    Ip_nu = (nu_x_inv - x_inv)*I_nu + I_nu_save;
    if(I_array  != (double *)0)  I_array[k] = I_nu;
    if(Ip_array != (double *)0) Ip_array[k] = Ip_nu;
    nu -= 1.0;
  }
  *I_end  = I_nu;
  *Ip_end = Ip_nu;
  return GSL_SUCCESS;
}


/* Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
 * to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
 */
static
int
bessel_K_CF2(const double nu, const double x, double * K_nu, double * K_nup1)
{
  const int maxiter = 10000;

  int i = 1;
  double bi = 2.0*(1.0 + x);
  double di = 1.0/bi;
  double delhi = di;
  double hi    = di;

  double qi   = 0.0;
  double qip1 = 1.0;

  double ai = -(0.25 - nu*nu);
  double a1 = ai;
  double ci = -ai;
  double Qi = -ai;

  double s = 1.0 + Qi*delhi;

  for(i=2; i<=maxiter; i++) {
    double dels;
    double tmp;
    ai -= 2.0*(i-1);
    ci  = -ai*ci/i;
    tmp  = (qi - bi*qip1)/ai;
    qi   = qip1;
    qip1 = tmp;
    Qi += ci*qip1;
    bi += 2.0;
    di  = 1.0/(bi + ai*di);
    delhi = (bi*di - 1.0) * delhi;
    hi += delhi;
    dels = Qi*delhi;
    s += dels;
    if(fabs(dels/s) < GSL_MACH_EPS) break;
  }
  
  hi *= -a1;
  
  *K_nu   = sqrt(M_PI/(2.0*x)) / s;
  *K_nup1 = *K_nu * (nu + x + 0.5 - hi)/x;
  if(i == maxiter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}


/* Forward recurrence for K,K' values.
 */
static
int
bessel_K_recur(const double nu_min, const double x, const int kmax,
               const double K_min, const double K_minp1,
               double * result_K_array, double * result_Kp_array)
{
  int k;

  if(result_K_array != (double *)0) {
    result_K_array[0] = K_min;
    if(kmax > 0) result_K_array[1] = K_minp1;
  }
  if(result_Kp_array != (double *)0) {
    result_Kp_array[0] = -K_minp1 + nu_min/x * K_min;
    if(kmax > 0) result_Kp_array[1] = -(nu_min + 1.0)/x * K_minp1 - K_min;
  }
  
  if(kmax > 0) {
    double K_k;
    double K_km1 = K_minp1;
    double K_km2 = K_min;
    for(k=2; k<=kmax; k++) {
      double nu   = nu_min + k - 1.0;
      double nuox = nu/x;
      K_k = K_km2 + 2.0*nuox * K_km1;
      if(result_K_array  != (double *)0)  result_K_array[k] = K_k;
      if(result_Kp_array != (double *)0) result_Kp_array[k] = -nuox*K_k-K_km1;
      K_km2 = K_km1;
      K_km1 = K_k;
    }
  }
}


/*
 * Assumes x > 0  and  nu_min >= 0.
 */
static
int
bessel_InuKnu_steed_scaled(const double nu_min, const double x, const int kmax,
                           double * result_I_array,  double * result_K_array,
                           double * result_Ip_array, double * result_Kp_array)
{
  int n;
  int k;
  int N     = (int)(nu_min + 0.5);
  double mu = nu_min - N;
  double I_small = 100.0*DBL_MIN;
  double I_ratio_numax;
  double I_ratio_numin;
  double I_scale_numin;
  double I_min, Ip_min;
  double I_min_save;
  double K_mu, K_mup1, K_min, K_minp1;
  double Kp_min;

  /* Evaluate I'/I at nu_max.
   */
  bessel_I_CF1(nu_min + kmax, x, &I_ratio_numax);

  /* Do backward recurrence for I and I' from nu_max to nu_min.
   */
  bessel_I_recur(nu_min, x, kmax,
                 I_small, I_small*I_ratio_numax,
		 &I_min, &Ip_min,
		 result_I_array, result_Ip_array
		 );
  I_ratio_numin = Ip_min/I_min;

  /* Evaluate K_mu and K_{mu+1}.
   */
  if(x < 2.0) {
    bessel_K_temme(mu, x, &K_mu, &K_mup1);
  }
  else {
    bessel_K_CF2(mu, x, &K_mu, &K_mup1);
  }

  /* Do forward recurrence to obtain K and K' at nu_min.
   */
  K_min = K_mu;
  K_minp1 = K_mup1;
  for(n=1; n<=N; n++) {
    double nu = mu + n;
    double K_minp1_save = K_minp1;
    K_minp1 = K_min + 2.0*nu/x*K_minp1;
    K_min   = K_minp1_save;
  }

  /* Use Wronskian relation to obtain correctly
   * normalized values for K,K' I,I' at nu_min,
   * and determine scaling factor that must be
   * applied to the I,I' values.
   */
  Kp_min = nu_min/x * K_min - K_minp1;
  I_min_save = I_min;
  I_min  = 1.0/(x*(I_ratio_numin*K_min - Kp_min));
  Ip_min = I_ratio_numin*I_min;
  I_scale_numin = I_min / I_min_save;

  /* Apply scaling factor to I,I'.
   */
  if(result_I_array != (double *)0) {
    for(k=0; k<=kmax; k++) {
      result_I_array[k] *= I_scale_numin;
    }
  }
  if(result_Ip_array != (double *)0) {
    for(k=0; k<=kmax; k++) {
      result_Ip_array[k] *= I_scale_numin;
    }
  }

  /* Do forward recurrence to obtain K,K' values.
   */
  bessel_K_recur(nu_min, x, kmax, K_min, K_minp1, result_K_array, result_Kp_array);
 
}

/* Assumes we can apply the uniform asymptotic result
 * to determine starting values for recursion.
 */
static
int
bessel_InuKnu_scaled_asymp(const double nu_min, const double x, const int kmax,
                           double * result_I_array,  double * result_K_array,
                           double * result_Ip_array, double * result_Kp_array)
{
  double nu_max = nu_min + kmax;
  double K_min, K_minp1;
  double I_max, I_maxp1;
  int stat_I   = gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nu_max,     x, &I_max);
  int stat_Ip1 = gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nu_max+1.0, x, &I_maxp1);
  int stat_K   = gsl_sf_bessel_Knu_scaled_asymp_unif_impl(nu_min,     x, &K_min);
  int stat_Kp1 = gsl_sf_bessel_Knu_scaled_asymp_unif_impl(nu_min+1.0, x, &K_minp1);
  double Ip_max = I_maxp1 + nu_max/x * I_max;
  double I_min, Ip_min;

  bessel_I_recur(nu_min, x, kmax,
                 I_max, Ip_max,
	         &I_min, &Ip_min,
	         result_I_array, result_Ip_array
	         );

  bessel_K_recur(nu_min, x, kmax,
                 K_min, K_minp1,
                 result_K_array, result_Kp_array);
}



int
gsl_sf_bessel_InuKnu_scaled_impl(double nu, double x, int kmax,
                                 double * result_I,  double * result_K,
                                 double * result_Ip, double * result_Kp)
{
  if(x < 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    int status = 0;
    int k;
    if(result_I != (double *)0) {
      for(k=0; k<=kmax; k++) result_I[k] = 0.0;
      if(nu == 0.0) result_I[k] = 1.0;
    }
    if(result_Ip != (double *)0) {
      for(k=0; k<=kmax; k++) result_I[k] = 0.0;
      if(nu < 1.0) status += 1;
      else if(nu == 1.0) result_I[0] = 0.5;
    }
    if(result_K != (double *)0) {
      for(k=0; k<=kmax; k++) result_K[k] = 0.0;
      status += 1;
    }
    if(result_Kp != (double *)0) {
      for(k=0; k<=kmax; k++) result_Kp[k] = 0.0;
      status += 1;
    }
    if(status) {
      return GSL_EDOM;
    }
    else {
      return GSL_SUCCESS;
    }
  }
  else if(locMin( 0.29/(nu*nu), 0.5/(nu*nu + x*x) ) < GSL_ROOT3_MACH_EPS) {
    return bessel_InuKnu_scaled_asymp(nu, x, kmax,
                                      result_I,  result_K,
                                      result_Ip, result_Kp
			              );
  }
  else {
    return bessel_InuKnu_steed_scaled(nu, x, kmax,
                                      result_I,  result_K,
                                      result_Ip, result_Kp
			              );
  }
}

int
gsl_sf_bessel_Inu_impl(double nu, double x, double * result)
{
  if(x < 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else if(x*x < 10.0*(nu+1.0)*GSL_ROOT5_MACH_EPS) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 4, &b);
    *result = ex * b;
    return stat;
  }
  else if(x*x < 10.0*(nu+1.0)) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 100, result);
    *result = ex * b;
    return stat;
  }
  else if(locMin( 0.29/(nu*nu), 0.5/(nu*nu + x*x) ) < GSL_ROOT3_MACH_EPS) {
    return gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nu, x, result);
  }
  else {
    /* FIXME */
  }
}


int
gsl_sf_bessel_Inu_e(double nu, double x, double * result)
{
  int status = gsl_sf_bessel_Inu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Inu_e", status);
  }
  return status;
}


double
gsl_sf_bessel_Inu(double nu, double x)
{
  double y;
  int status = gsl_sf_bessel_Inu_impl(nu, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Inu", status);
  }
  return y;
}
