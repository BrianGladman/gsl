/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_poly.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_ellint.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_bessel.h"
#include "gsl_sf_hyperg.h"
#include "gsl_sf_legendre.h"

#define Root_2OverPi_  0.797884560802865355879892
#define locEPS         (1000.0*GSL_MACH_EPS)

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

#define locMAX(a,b)    ((a) > (b) ? (a) : (b))

#define RECURSE_LARGE  (1.0e-5*DBL_MAX)
#define RECURSE_SMALL  (1.0e+5*DBL_MIN)


/* Backward recursion for
 * y_n := 1/Gamma[1+mu+n] P^{mu+n}_{-1/2 + I tau}(x)
 * x > 1
 * mu_min > -1
 * n = N, N-1, ..., 0
 *
 * checked OK: [GJ] Tue Sep 15 19:25:35 MDT 1998 
 */
static
int
backward_recurse_pos_mu_xgt1(const double mu_min, const double tau, const double x,
                             const int N,
                             const double y_N, const double y_Np1,
			     double * result_y_0, double * result_y_1)
{
  const double c1 = 1.0/sqrt(x*x-1.0);
  const double t2 = tau*tau;
  double y_np1 = y_Np1;
  double y_n   = y_N;
  double y_nm1;
  int n;

  for(n=N; n>=1; n--) {
    const double mu = mu_min + n;
    const double d  = (mu - 0.5)*(mu-0.5) + t2;
    y_nm1 = -mu*(mu+1.0)/d * (y_np1 + 2.0*mu*x*c1*y_n/(mu+1.0));
    y_np1 = y_n;
    y_n   = y_nm1;
    if(fabs(y_nm1) > RECURSE_LARGE) {
      *result_y_0 = 0.0;
      *result_y_1 = 0.0;
      return GSL_EOVRFLW;
    }
  }
  
  *result_y_0 = y_n;
  *result_y_1 = y_np1;
  return GSL_SUCCESS;
}
			 

/* Backward recursion for
 * y_n := 1/Gamma[1-mu-n] P^{-mu-n}_{-1/2 + I tau}(x)
 * x > 1
 * mu_max < 0
 * n = N, N-1, ..., 0
 * mu != integer
 *
 * checked OK: Tue Sep 15 19:53:20 MDT 1998 
 */
static
int
backward_recurse_neg_mu_xgt1(const double mu, const double tau, const double x,
                             const int N,
                             const double y_N, const double y_Np1,
			     double * result_y_0, double * result_y_1)
{ 
  const double c1 = 1.0/sqrt(x*x-1.0);
  const double t2 = tau*tau;
  double y_np1 = y_Np1;
  double y_n   = y_N;
  double y_nm1;
  int n;

  for(n=N; n>=1; n--) {
    const double mupn = mu + n;
    const double d = (mupn+0.5)*(mupn+0.5) + t2;
    y_nm1 = -2.0*mupn*x*c1/(mupn-1.0)*y_n - d/(mupn*(mupn-1.0))*y_np1;
    y_np1 = y_n;
    y_n   = y_nm1;
    if(fabs(y_nm1) > RECURSE_LARGE) {
      *result_y_0 = 0.0;
      *result_y_1 = 0.0;
      return GSL_EOVRFLW;
    }
  }
  
  *result_y_0 = y_n;
  *result_y_1 = y_np1;
  return GSL_SUCCESS;
}


/* Backward recursion for
 * y_n := 1/Gamma[1+mu+n] P^{mu+n}_{-1/2 + I tau}(x)
 * |x| < 1
 * mu_min > -1
 * n = N, N-1, ..., 0
 *
 * Note that this is not stable for x < 0, because it is dominated
 * by the Q part of the solution there. Use this for 0 < x < 1.
 *
 * checked OK: Wed Sep 16 19:14:27 MDT 1998 
 */
static
int
backward_recurse_pos_mu_xlt1(const double mu_min, const double tau, const double x,
                             const int N,
                             const double y_N, const double y_Np1,
			     double * result_y_0, double * result_y_1)
{
  const double c1 = 1.0/sqrt(1.0-x*x);
  const double t2 = tau*tau;
  double y_np1 = y_Np1;
  double y_n   = y_N;
  double y_nm1;
  int n;

  for(n=N; n>=1; n--) {
    const double mu = mu_min + n;
    const double d  = (mu - 0.5)*(mu-0.5) + t2;
    y_nm1 = mu*(mu+1.0)/d * (y_np1 + 2.0*mu*x*c1*y_n/(mu+1.0));
    y_np1 = y_n;
    y_n   = y_nm1;
    if(fabs(y_nm1) > RECURSE_LARGE) {
      *result_y_0 = 0.0;
      *result_y_1 = 0.0;
      return GSL_EOVRFLW;
    }
  }
  
  *result_y_0 = y_n;
  *result_y_1 = y_np1;
  return GSL_SUCCESS;
}


/* Forward recursion for
 * y_n := 1/Gamma[1+mu+n] P^{mu+n}_{-1/2 + I tau}(x)
 * |x| < 1
 * mu_min > -1
 * n = N, N-1, ..., 0
 *
 * Note that this is not stable for x > 0. Use this for -1 < x < 0.
 *
 * checked OK: Wed Sep 16 19:40:29 MDT 1998 
 */
static
int
forward_recurse_pos_mu_xlt1(const double mu_min, const double tau, const double x,
                            const int N,
                            const double y_0, const double y_1,
			    double * result_yN, double * result_yNp1)
{
  const double c1 = 1.0/sqrt(1.0-x*x);
  const double t2 = tau*tau;
  double y_nm1 = y_0;
  double y_n   = y_1;
  double y_np1;
  int n;

  for(n=1; n<=N; n++) {
    const double mu = mu_min + n;
    const double d  = (mu-0.5)*(mu-0.5) + t2;
    y_np1 = (d/mu*y_nm1 - 2.0*mu*x*c1*y_n)/(mu+1.0);
    y_nm1 = y_n;
    y_n   = y_np1;
    if(fabs(y_np1) > RECURSE_LARGE) {
      *result_yN   = 0.0;
      *result_yNp1 = 0.0;
      return GSL_EOVRFLW;
    }
  }

  *result_yN   = y_nm1;
  *result_yNp1 = y_n;
  return GSL_SUCCESS;
}


/* Backward recursion for
 * y_n := 1/Gamma[1-mu-n] P^{-mu-n}_{-1/2 + I tau}(x)
 * |x| < 1
 * mu_max < 0
 * n = N, N-1, ..., 0
 * mu != integer
 *
 * Note that this is not stable for x < 0, because it is dominated
 * by the Q part of the solution there. Use this for 0 < x < 1.
 *
 * checked OK: Wed Sep 16 20:11:21 MDT 1998
 */
static
int
backward_recurse_neg_mu_xlt1(const double mu, const double tau, const double x,
                             const int N,
                             const double y_N, const double y_Np1,
			     double * result_y_0, double * result_y_1)
{ 
  const double c1 = 1.0/sqrt(1.0-x*x);
  const double t2 = tau*tau;
  double y_np1 = y_Np1;
  double y_n   = y_N;
  double y_nm1;
  int n;

  for(n=N; n>=1; n--) {
    const double mupn = mu + n;
    const double d = (mupn+0.5)*(mupn+0.5) + t2;
    y_nm1 = -2.0*mupn*x*c1/(mupn-1.0)*y_n + d/(mupn*(mupn-1.0))*y_np1;
    y_np1 = y_n;
    y_n   = y_nm1;
    if(fabs(y_nm1) > RECURSE_LARGE) {
      *result_y_0 = 0.0;
      *result_y_1 = 0.0;
      return GSL_EOVRFLW;
    }
  }
  
  *result_y_0 = y_n;
  *result_y_1 = y_np1;
  return GSL_SUCCESS;
}


/* Forward recursion for
 * y_n := 1/Gamma[1-mu-n] P^{-mu-n}_{-1/2 + I tau}(x)
 * |x| < 1
 * mu_min > -1
 * n = N, N-1, ..., 0
 *
 * Note that this is not stable for x > 0. Use this for -1 < x < 0.
 *
 * checked OK: Wed Sep 16 20:32:55 MDT 1998
 */
static
int
forward_recurse_neg_mu_xlt1(const double mu, const double tau, const double x,
                            const int N,
                            const double y_0, const double y_1,
			    double * result_yN, double * result_yNp1)
{
  const double c1 = 1.0/sqrt(1.0-x*x);
  const double t2 = tau*tau;
  double y_nm1 = y_0;
  double y_n   = y_1;
  double y_np1;
  int n;

  for(n=1; n<=N; n++) {
    const double mupn = mu + n;
    const double d = (mupn+0.5)*(mupn+0.5) + t2;
    y_np1 = (y_nm1 + 2.0*mupn*x*c1/(mupn-1.0)*y_n) * mupn*(mupn-1.0)/d;
    y_nm1 = y_n;
    y_n   = y_np1;
    if(fabs(y_np1) > RECURSE_LARGE) {
      *result_yN   = 0.0;
      *result_yNp1 = 0.0;
      return GSL_EOVRFLW;
    }
  }

  *result_yN   = y_nm1;
  *result_yNp1 = y_n;
  return GSL_SUCCESS;
}

#if 0
void test_recurse(void)
{
  /*
  double mu_min = 0.0;
  double tau    = 5.0;
  double x      = -0.5;
  int N = 1000;
  double y_n   =  3.815801082323749413e+241;
  double y_np1 =  6.602600486662392861e+241;
  double y_0, y_1;
  backward_recurse_pos_mu_xlt1(mu_min, tau, x,
                               N,
                               y_n, y_np1,
			       &y_0, &y_1);
  printf("%24.18g   %24.18g    %24.18g  %24.18g\n", y_n, y_np1, y_0, y_1);
  exit(0);
  */

  /*
  double mu_min = 0.0;
  double tau    = 5.0;
  double x      = -0.5;
  int N = 1000;
  double y_n;
  double y_np1;
  double y_0 = 6693.103933518211835;
  double y_1 = 35192.90653374874763;
  forward_recurse_pos_mu_xlt1(mu_min, tau, x,
                               N,
                               y_0, y_1,
			       &y_n, &y_np1);
  printf("%24.18g   %24.18g    %24.18g  %24.18g\n", y_0, y_1, y_n, y_np1);
  exit(0);
  */

  /*
  double mu  = 1.0/3.0;
  double tau = 5.0;
  double x   = 3.0;
  int N = 500;
  double y_n   =   2.5799943898780086603e-79 ;
  double y_np1 =  -1.8208744513179016911e-79 ;
  double y_0, y_1;
  backward_recurse_neg_mu_xgt1(mu, tau, x,
                               N,
                               y_n, y_np1,
			       &y_0, &y_1);
  printf("%24.18g   %24.18g    %24.18g  %24.18g\n", y_n, y_np1, y_0, y_1);
  exit(0);
  */

  /*
  double mu  = 1.0/3.0;
  double tau = 5.0;
  double x   = 0.5;
  int N = 20;
  double y_n   = ;
  double y_np1 = ;
  double y_0, y_1;
  backward_recurse_neg_mu_xlt1(mu, tau, x,
                               N,
                               y_n, y_np1,
			       &y_0, &y_1);
  printf("%24.18g   %24.18g    %24.18g  %24.18g\n", y_n, y_np1, y_0, y_1);
  exit(0);
  */
  
  /*
  double mu_min = 1.0/3.0;
  double tau    = 5.0;
  double x      = -0.5;
  int N = 1000;
  double y_n;
  double y_np1;
  double y_0 =  2906.1619891014585047;
  double y_1 = -206.00701262658096387;
  forward_recurse_neg_mu_xlt1(mu_min, tau, x,
                              N,
                              y_0, y_1,
			      &y_n, &y_np1);
  printf("%24.18g   %24.18g   %24.18g  %24.18g\n", y_0, y_1, y_n, y_np1);
  exit(0);
  */
}
#endif /* 0 */


/* Implementation of large negative mu asymptotic
 * [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326]
 */

static inline double olver_U1(double beta2, double p)
{
  return (p-1.0)/(24.0*(1.0+beta2)) * (3.0 + beta2*(2.0 + 5.0*p*(1.0+p)));
}

static inline double olver_U2(double beta2, double p)
{
  double beta4 = beta2*beta2;
  double p2    = p*p;
  double poly1 =  4.0*beta4 + 84.0*beta2 - 63.0;
  double poly2 = 16.0*beta4 + 90.0*beta2 - 81.0;
  double poly3 = beta2*p2*(97.0*beta2 - 432.0 + 77.0*p*(beta2-6.0) - 385.0*beta2*p2*(1.0 + p));
  return (1.0-p)/(1152.0*(1.0+beta2)) * (poly1 + poly2 + poly3);
}

static const double U3c1[] = {   -1307.0,   -1647.0,    3375.0,    3675.0 };
static const double U3c2[] = {   29366.0,   35835.0, -252360.0, -272630.0,
                                276810.0,  290499.0 };
static const double U3c3[] = {  -29748.0,   -8840.0, 1725295.0, 1767025.0,
                              -7313470.0, -754778.0, 6309875.0, 6480045.0 };
static const double U3c4[] = {    2696.0,    -16740.0,   -524250.0,  -183975.0,
                              14670540.0,  14172939.0, -48206730.0, -48461985.0,
			      36756720.0,  37182145.0 };
static const double U3c5[] = {       9136.0,      22480.0,     12760.0,
                                  -252480.0,   -4662165.0,   -1705341.0,
				 92370135.0,   86244015.0, -263678415.0,
			       -260275015.0, 185910725.0,  185910725.0 };

static double olver_U3(double beta2, double p)
{
  double beta4 = beta2*beta2;
  double beta6 = beta4*beta2;
  double opb2s = (1.0+beta2)*(1.0+beta2);
  double den   = 39813120.0 * opb2s*opb2s;
  double poly1 = gsl_sf_poly_eval(U3c1, 4, p);
  double poly2 = gsl_sf_poly_eval(U3c2, 6, p);
  double poly3 = gsl_sf_poly_eval(U3c3, 8, p);
  double poly4 = gsl_sf_poly_eval(U3c4, 10, p);
  double poly5 = gsl_sf_poly_eval(U3c5, 12, p);
  
  return (p-1.0)*(     1215.0*poly1 + 324.0*beta2*poly2
                 + 54.0*beta4*poly3 +  12.0*beta6*poly4
		 + beta4*beta4*poly5
		 ) / den;
}

/* Large negative mu asymptotic
 * P^{-mu}_{-1/2 + I tau}, mu -> Inf
 * |x| < 1
 *
 * [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326]
 */
static
int
conicalP_xlt1_large_neg_mu(double mu, double tau, double x, double * result)
{
  double beta  = tau/mu;
  double beta2 = beta*beta;
  double S     = beta * acos((1.0-beta2)/(1.0+beta2));
  double p     = x/sqrt(beta2*(1.0-x*x) + 1.0);
  double lg_mup1;
  int lg_stat = gsl_sf_lngamma_impl(mu+1.0, &lg_mup1);
  double ln_pre_1 =  0.5*mu*(S - log(1.0+beta2) + log((1.0-p)/(1.0+p))) - lg_mup1;
  double ln_pre_2 = -0.25 * log(1.0 + beta2*(1.0-x));
  double ln_pre_3 = -tau * atan(p*beta);
  double ln_pre = ln_pre_1 + ln_pre_2 + ln_pre_3;
  
  if(ln_pre > GSL_LOG_DBL_MAX) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(ln_pre < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    double sum = 1.0 - olver_U1(beta2, p)/mu + olver_U2(beta2, p)/(mu*mu);
    *result = exp(ln_pre) * sum;
    return lg_stat;
  }
}



/* Implementation of large tau asymptotic
 *
 * A_n^{-mu}, B_n^{-mu}  [Olver, p.465, 469]
 */

#ifdef HAVE_INLINE
inline
#endif
static double olver_B0_xi(double mu, double xi)
{
  return (1.0 - 4.0*mu*mu)/(8.0*xi) * (1.0/tanh(xi) - 1.0/xi);
}

static double olver_A1_xi(double mu, double xi, double x)
{
  double B = olver_B0_xi(mu, xi);
  double psi = (4.0*mu*mu - 1.0)/16.0 * (1.0/(x*x-1.0) - 1.0/(xi*xi));
  return 0.5*xi*xi*B*B + (mu+0.5)*B - psi + mu/6.0*(0.25 - mu*mu);
}

#ifdef HAVE_INLINE
inline
#endif
static double olver_B0_th(double mu, double theta)
{
  return -(1.0 - 4.0*mu*mu)/(8.0*theta) * (1.0/tan(theta) - 1.0/theta);
}

static double olver_A1_th(double mu, double theta, double x)
{
  double B = olver_B0_th(mu, theta);
  double psi = (4.0*mu*mu - 1.0)/16.0 * (1.0/(x*x-1.0) + 1.0/(theta*theta));
  return -0.5*theta*theta*B*B + (mu+0.5)*B - psi + mu/6.0*(0.25 - mu*mu);
}


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}
 * 1 < x
 * tau -> Inf 
 * [Olver, p. 469]
 */
static
int
conicalP_xgt1_neg_mu_largetau_impl(const double mu, const double tau,
                                   const double x, double * result)
{
  double xi = acosh(x);
  double xi_pre;
  double pre;
  double sumA, sumB;
  double arg;
  double J_mup1, J_mu, J_mum1;

  if(xi < GSL_ROOT4_MACH_EPS) {
    xi_pre = 1.0 - xi*xi/6.0;
  }
  else {
    xi_pre = xi/sinh(xi);
  }

  pre = sqrt(xi_pre) * pow(1.0/tau, mu);
  if(pre == 0.0) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }

  arg = tau*xi;

  gsl_sf_bessel_Jnu_impl(mu + 1.0,   arg, &J_mup1);
  gsl_sf_bessel_Jnu_impl(mu,         arg, &J_mu);
  J_mum1 = -J_mup1 + 2.0*mu/arg; /* careful of mu < 1 */

  sumA = 1.0 - olver_A1_xi(-mu, xi, x)/(tau*tau);
  sumB = olver_B0_xi(-mu, xi);

  *result = pre * (J_mu * sumA - xi/tau * J_mum1 * sumB);
  return GSL_SUCCESS; /* FIXME: hmmm, success??? */
}


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}
 * -1 < x < 1
 * tau -> Inf 
 * [Olver, p. 473]
 */
static
int
conicalP_xlt1_neg_mu_largetau_impl(const double mu, const double tau,
                                   const double x, double * result)
{
  double theta = acos(x);
  double th_pre;
  double pre;
  double sumA, sumB;
  double arg;
  double I_mup1, I_mu, I_mum1;

  if(theta < GSL_ROOT4_MACH_EPS) {
    th_pre = 1.0 + theta*theta/6.0;
  }
  else {
    th_pre = theta/sin(theta);
  }

  pre = sqrt(th_pre) * pow(1.0/tau, mu);
  if(pre == 0.0) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }

  arg = tau*theta;

  gsl_sf_bessel_Inu_impl(mu + 1.0,   arg, &I_mup1);
  gsl_sf_bessel_Inu_impl(mu,         arg, &I_mu);
  I_mum1 = I_mup1 + 2.0*mu/arg; /* careful of mu < 1 */

  sumA = 1.0 - olver_A1_th(-mu, theta, x)/(tau*tau);
  sumB = olver_B0_th(-mu, theta);

  *result = pre * (I_mu * sumA - theta/tau * I_mum1 * sumB);
  return GSL_SUCCESS;
}


/* P^{mu}_{-1/2 + I tau}  first hypergeometric representation
 * -1 < x < 1
 * more effective for |x| small
 *
 * [Kolbig,   (3)] (note typo in args of gamma functions)
 * [Bateman, (22)] (correct form)
 */
static
int
conicalP_xlt1_hyperg_A(double mu, double tau, double x, double * result)
{
  double x2 = x*x;
  double pre  = M_SQRTPI / pow(0.5*sqrt(1-x2), mu);
  double ln_g1, ln_g2, arg_g1, arg_g2;
  double pre1, pre2, F1, F2;
  int status;

  int stat_F1 = gsl_sf_hyperg_2F1_conj_impl(0.25 - 0.5*mu, 0.5*tau, 0.5, x2, &F1);
  int stat_F2 = gsl_sf_hyperg_2F1_conj_impl(0.75 - 0.5*mu, 0.5*tau, 1.5, x2, &F2);

  if(stat_F1 == GSL_SUCCESS && stat_F2 == GSL_SUCCESS) {
    status = GSL_SUCCESS;
  }
  else if(stat_F1 == GSL_ELOSS || stat_F2 == GSL_ELOSS) {
    status = GSL_ELOSS;
  }
  else {
    status = GSL_EFAILED;
  }

  gsl_sf_lngamma_complex_impl(0.75 - 0.5*mu, -0.5*tau, &ln_g1, &arg_g1);
  gsl_sf_lngamma_complex_impl(0.25 - 0.5*mu, -0.5*tau, &ln_g2, &arg_g2);
  
  pre1 =        exp(-2.0*ln_g1);
  pre2 = -2.0*x*exp(-2.0*ln_g2);
  
  *result = pre * (pre1 * F1 + pre2 * F2);
  return status;
}


/* P^{mu}_{-1/2 + I tau}
 * defining hypergeometric representation
 * [Abramowitz+Stegun, 8.1.2]
 * 1 < x < 3
 * effective for x near 1
 *
 */
static
int
conicalP_def_hyperg(double mu, double tau, double x, double * result)
{
  double F;
  int stat_F = gsl_sf_hyperg_2F1_conj_renorm_impl(0.5, tau, 1.0-mu, 0.5*(1.0-x), &F);
  *result = pow((x+1.0)/(x-1.0), 0.5*mu) * F;
  return stat_F;
}


/* P^{mu}_{-1/2 + I tau}  second hypergeometric representation
 * [Zhurina+Karmazina, (3.1)] 
 * -1 < x < 3
 * effective for x near 1
 *
 */
static
int
conicalP_xnear1_hyperg_C(double mu, double tau, double x, double * result)
{
  double ln_pre, arg_pre;
  double ln_g1, arg_g1;
  double ln_g2, arg_g2;
  double F;

  int stat_F = gsl_sf_hyperg_2F1_conj_renorm_impl(0.5+mu, tau, 1.0+mu, 0.5*(1.0-x), &F);

  gsl_sf_lngamma_complex_impl(0.5+mu, tau, &ln_g1, &arg_g1);
  gsl_sf_lngamma_complex_impl(0.5-mu, tau, &ln_g2, &arg_g2);

  ln_pre  = mu*M_LN2 - 0.5*mu*log(fabs(x*x-1.0)) + ln_g1 - ln_g2;
  arg_pre = arg_g1 - arg_g2;

  *result = exp(ln_pre) * F;
  return stat_F;
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/


/* P^0_{-1/2 + I lambda}
 */
int
gsl_sf_conical_0_impl(const double lambda, const double x, double * result)
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(lambda == 0.0) {
    double K;
    int stat_K;
    if(x < 1.0) {
      const double th = acos(x);
      const double s  = sin(0.5*th);
      stat_K = gsl_sf_ellint_Kcomp_impl(s, locEPS, &K);
      *result = 2.0/M_PI * K;
      return stat_K;
    }
    else {
      const double xi = acosh(x);
      const double c  = cosh(0.5*xi);
      const double t  = tanh(0.5*xi);
      stat_K = gsl_sf_ellint_Kcomp_impl(t, locEPS, &K);
      *result = 2.0/M_PI / c * K;
      return stat_K;
    }
  }
  else if(    x <= 0.0
          || (x < 0.1 && lambda < 17.0)
	  || (x < 0.2 && lambda < 5.0 )
	  ) {
    return conicalP_xlt1_hyperg_A(0.0, lambda, x, result);
  }
  else if(   ((0.1 < x && x <= 0.2) && (5.0 < lambda && lambda < 17.0))
          || ((0.2 < x && x < 1.5)  && (lambda < 20.0))
          ) {
  }
  else if(1.5 <= x && lambda < locMax(x,20.0)) {
  }
  else {
  }
}


/* P^1_{-1/2 + I lambda}
 */
int
gsl_sf_conical_1_impl(const double lambda, const double x, double * result)
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(lambda == 0.0) {
    double K, E;
    int stat_K, stat_E;
    if(fabs(x-1.0) < 10.0*GSL_MACH_EPS) {
      /* FIXME */
    }
    else if(x < 1.0) {
      const double th = acos(x);
      const double s  = sin(0.5*th);
      const double c2 = 1.0 - s*s;
      const double sth = sin(th);
      stat_K = gsl_sf_ellint_Kcomp_impl(s, locEPS, &K);
      stat_E = gsl_sf_ellint_Ecomp_impl(s, locEPS, &E);
      *result = 2.0/(M_PI*sth) * (E - c2 * K);
      return stat_K;
    }
    else {
      const double xi = acosh(x);
      const double c  = cosh(0.5*xi);
      const double t  = tanh(0.5*xi);
      const double sxi = sinh(xi);
      stat_K = gsl_sf_ellint_Kcomp_impl(t, locEPS, &K);
      stat_E = gsl_sf_ellint_Ecomp_impl(t, locEPS, &E);
      *result = 2.0/(M_PI*sxi) * c * (E - K);
      return stat_K;
    }
  }
  else if(    x <= 0.0
          || (x < 0.1 && lambda < 17.0)
	  || (x < 0.2 && lambda < 5.0 )
	  ) {
    return conicalP_xlt1_hyperg_A(1.0, lambda, x, result);
  }
  else if(   ((0.1 < x && x <= 0.2) && (5.0 < lambda && lambda < 17.0))
          || ((0.2 < x && x < 1.5)  && (lambda < 20.0))
          ) {
  }
  else if(1.5 <= x && lambda < locMax(x,20.0)) {
  }
  else {
  }
}


/* P^{1/2}_{-1/2 + I lambda} (x)
 * [Abramowitz+Stegun 8.6.8, 8.6.12]
 * checked OK [GJ] Fri May  8 12:24:36 MDT 1998 
 */
int gsl_sf_conical_half_impl(const double lambda, const double x,
                             double * result
                             )
{
  /* FIXME: do something about this nonsense */
  double one_plus_x  = 1.0 + x;
  double one_minus_x = 1.0 - x;
  /* double x = 0.5 * (one_plus_x - one_minus_x); */
  if(fabs(x) < 1.0) {
    double ac  = acos(x);
    double den = sqrt(sqrt(one_minus_x*one_plus_x));
    *result = Root_2OverPi_ / den * cosh(ac * lambda);
    return GSL_SUCCESS;
  }
  else if(one_minus_x < 0.0) { /* x > 1 */
    double sq_term = sqrt(-one_minus_x*one_plus_x);
    double ln_term = log(x + sq_term);
    double den = sqrt(sq_term);
    *result = Root_2OverPi_ / den * cos(lambda * ln_term);
    return GSL_SUCCESS;
  }
  else {
    /* |x| == 1 or x < -1 */
    *result = 0.0;
    return GSL_EDOM;
  }
}


/* P^{-1/2}_{-1/2 + I lambda} (x)
 * [Abramowitz+Stegun 8.6.9, 8.6.14]
 * checked OK [GJ] Fri May  8 12:24:43 MDT 1998 
 */
int gsl_sf_conical_mhalf_impl(const double lambda, const double x,
                              double * result
                              )
{
  /* FIXME: do something about this nonsense */
  double one_plus_x  = 1.0 + x;
  double one_minus_x = 1.0 - x;
  /* double x = 0.5 * (one_plus_x - one_minus_x); */
  if(fabs(x) < 1.0) {
    double ac  = acos(x);
    double den = sqrt(sqrt(one_minus_x*one_plus_x));
    double arg = ac * lambda;
    if(fabs(arg) < GSL_SQRT_MACH_EPS) {
      *result = Root_2OverPi_ / den * ac;
    }
    else {
      *result = Root_2OverPi_ / (den*lambda) * sinh(arg);
    }
    return GSL_SUCCESS;
  }
  else if(one_minus_x < 0.0) { /* x > 1 */
    double sq_term = sqrt(-one_minus_x*one_plus_x);
    double ln_term = log(x + sq_term);
    double den = sqrt(sq_term);
    double arg = lambda * ln_term;
    if(arg < GSL_SQRT_MACH_EPS) {
      *result = Root_2OverPi_ / den * ln_term;
    }
    else {
      *result = Root_2OverPi_ / (den*lambda) * sin(arg);
    }
    return GSL_SUCCESS;
  }
  else if(one_minus_x == 0.0) { /* x == 1 */
    *result = 0.0;
    return GSL_SUCCESS;
  }  
  else {
    /* x <= -1 */
    *result = 0.0;
    return GSL_EDOM;
  }
}

int gsl_sf_conical_sph_reg_impl(const int lmax, const double lambda,
                                const double x,
                                double * result
				)
{
  /* FIXME: do something about this nonsense */
  double one_p_x  = 1.0 + x;
  double one_m_x = 1.0 - x;
  /* double x = 0.5 * (one_plus_x - one_minus_x); */

  if(fabs(x) < 1.0) {
    double f0;
    double p[2];
    int l_start = 20 + (int) ceil(lmax * (1. + (0.14 + 0.026*lambda)));
    p[0] = lambda*lambda;
    p[1] = x/sqrt(one_m_x*one_p_x);
    gsl_sf_conical_mhalf_impl(lambda, x, &f0);  /* l =  0  */ 
    /*recurse_backward_minimal_simple_conical_sph_reg_xlt1(l_start, lmax, 0, p, f0, harvest, result);
    */
    
  }
  else if(x > 1.0) {
    double f0;
    double p[2];
    int l_start = 10 + (int) ceil(lmax * (1. + (0.14 + 0.026*lambda)*(x-1.)));
    p[0] = lambda*lambda;
    p[1] = x/sqrt(-one_m_x*one_p_x);
    gsl_sf_conical_mhalf_impl(lambda, x, &f0);
    /*recurse_backward_minimal_simple_conical_sph_reg_xgt1(l_start, lmax, 0, p, f0, harvest, result);*/
  }
  else {
    return GSL_EDOM;
  }
}

int gsl_sf_hyper_0_impl(const double lambda, const double x, double * result)
{
  *result = sin(lambda*x)/(lambda*sinh(x));
}

int gsl_sf_hyper_1_impl(const double lambda, const double x, double * result)
{
  *result = sin(lambda*x)/(lambda*sinh(x)) 
	  /sqrt(lambda*lambda+1.) * (1./tanh(x) - lambda/tan(lambda*x));
}

int gsl_sf_hyper_array_impl(int lmax, double lambda, double x, double * result, double * harvest)
{
  double X = 1./tanh(x);
  double y2, y1, y0;
  int ell;

  gsl_sf_hyper_0_impl(lambda, x, &y2);
  gsl_sf_hyper_1_impl(lambda, x, &y1);

  harvest[0] = y2;
  harvest[1] = y1;

  for(ell=2; ell<=lmax; ell++) {
    double a = sqrt(lambda*lambda + ell*ell);
    double b = sqrt(lambda*lambda + (ell-1)*(ell-1));
    y0 = ((2*ell-1)*X*y1 - b*y2) / a;
    y2 = y1;
    y1 = y0;
    harvest[ell] = y0;
  }
  
  *result = y0;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_conical_half_e(const double lambda, const double x, double * result)
{
  int status = gsl_sf_conical_half_impl(lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_half_e", status);
  }
  return status;
}

int gsl_sf_conical_mhalf_e(const double lambda, const double x, double * result)
{
  int status = gsl_sf_conical_mhalf_impl(lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_mhalf_e", status);
  }
  return status;
}

int gsl_sf_conical_sph_reg_e(const int l, const double lambda, const double x, double * result)
{
  int status = gsl_sf_conical_sph_reg_impl(l, lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_sph_reg_e", status);
  }
  return status;
}

int gsl_sf_conical_sph_reg_array_e(const int l, const double lambda, const double x, double * result_array)
{
  int status = gsl_sf_conical_sph_reg_array_impl(l, lambda, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_sph_reg_array_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_conical_half(const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conical_half_impl(lambda, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conical_half", status);
  }
  return y;
}

double gsl_sf_conical_mhalf(const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conical_mhalf_impl(lambda, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conical_mhalf", status);
  }
  return y;
}

double gsl_sf_conical_sph_reg(const int l, const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conical_sph_reg_impl(l, lambda, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conical_sph_reg", status);
  }
  return y;
}
