/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "legendre.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_poly.h"
#include "gsl_sf_trig.h"
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


/* Continued fraction for f_{ell+1}/f_ell
 * f_ell := P^{-mu-ell}_{-1/2 + I tau}(x),  x < 1.0
 *
 * Uses standard CF method from Temme's book.
 */
static
int
conicalP_negmu_xlt1_CF1(const double mu, const int ell, const double tau,
                        const double x, double * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int n = 1;
  double xi = x/(sqrt(1.0-x)*sqrt(1.0+x));
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = 1.0;
  double b1 = 2.0*(mu + ell + 1.0) * xi;
  double An = b1*Anm1 + a1*Anm2;
  double Bn = b1*Bnm1 + a1*Bnm2;
  double an, bn;
  double fn = An/Bn;

  while(n < maxiter) {
    double old_fn;
    double del;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    an = tau*tau + (mu - 0.5 + ell + n)*(mu - 0.5 + ell + n);
    bn = 2.0*(ell + mu + n) * xi;
    An = bn*Anm1 + an*Anm2;
    Bn = bn*Bnm1 + an*Bnm2;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An /= RECUR_BIG;
      Bn /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
    }

    old_fn = fn;
    fn = An/Bn;
    del = old_fn/fn;
    
    if(fabs(del - 1.0) < 3.0*GSL_MACH_EPS) break;
  }

  *result = fn;

  if(n == maxiter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* Continued fraction for f_{ell+1}/f_ell
 * f_ell := P^{-mu-ell}_{-1/2 + I tau}(x),  x >= 1.0
 *
 * Uses Gautschi (Euler) equivalent series.
 */
static
int
conicalP_negmu_xgt1_CF1(const double mu, const int ell, const double tau,
                        const double x, double * result)
{ 
  const int maxk = 20000;
  const double gamma = 1.0-1.0/(x*x);
  double tk   = 1.0;
  double sum  = 1.0;
  double rhok = 0.0;
  int k;
 
  for(k=1; k<maxk; k++) {
    double tlk = 2.0*(ell + mu + k);
    double l1k = (ell + mu - 0.5 + 1.0 + k);
    double ak = -(tau*tau + l1k*l1k)/(tlk*(tlk+2.0)) * gamma;
    rhok = -ak*(1.0 + rhok)/(1.0 + ak*(1.0 + rhok));
    tk  *= rhok;
    sum += tk;
    if(fabs(tk/sum) < 3.0*GSL_MACH_EPS) break;
  }

  *result = sqrt(x-1.0)*sqrt(x+1.0) / (x*(2.0*(ell+mu+1.0))) * sum;

  if(k == maxk)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* Implementation of large negative mu asymptotic
 * [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326]
 */

#ifdef HAVE_INLINE
inline
#endif
static double olver_U1(double beta2, double p)
{
  return (p-1.0)/(24.0*(1.0+beta2)) * (3.0 + beta2*(2.0 + 5.0*p*(1.0+p)));
}

#ifdef HAVE_INLINE
inline
#endif
static double olver_U2(double beta2, double p)
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

#if 0
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
#endif /* 0 */


/* Large negative mu asymptotic
 * P^{-mu}_{-1/2 + I tau}, mu -> Inf
 * |x| < 1
 *
 * [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326]
 */
int
gsl_sf_conicalP_xlt1_large_neg_mu_impl(double mu, double tau, double x,
                                       double * result, double * ln_multiplier)
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
  double sum   = 1.0 - olver_U1(beta2, p)/mu + olver_U2(beta2, p)/(mu*mu);

  if(sum == 0.0) {
    *result = 0.0;
    *ln_multiplier = 0.0;
    return GSL_SUCCESS;
  }
  else {
    int stat_e = gsl_sf_exp_mult_impl(ln_pre, sum, result);
    if(stat_e != GSL_SUCCESS) {
      *result = sum;
      *ln_multiplier = ln_pre;
    }
    else {
      *ln_multiplier = 0.0;
    }
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
  double psi;
  if(fabs(x - 1.0) < GSL_ROOT4_MACH_EPS) {
    double y = x - 1.0;
    double s = -1.0/3.0 + y*(2.0/15.0 - y *(61.0/945.0 - 452.0/14175.0*y));
    psi = (4.0*mu*mu - 1.0)/16.0 * s;
  }
  else {
    psi = (4.0*mu*mu - 1.0)/16.0 * (1.0/(x*x-1.0) - 1.0/(xi*xi));
  }
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
  double psi;
  if(fabs(x - 1.0) < GSL_ROOT4_MACH_EPS) {
    double y = 1.0 - x;
    double s = -1.0/3.0 + y*(2.0/15.0 - y *(61.0/945.0 - 452.0/14175.0*y));
    psi = (4.0*mu*mu - 1.0)/16.0 * s;
  }
  else {
    psi = (4.0*mu*mu - 1.0)/16.0 * (1.0/(x*x-1.0) + 1.0/(theta*theta));
  }
  return -0.5*theta*theta*B*B + (mu+0.5)*B - psi + mu/6.0*(0.25 - mu*mu);
}


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}
 * 1 < x
 * tau -> Inf 
 * [Olver, p. 469]
 */
int
gsl_sf_conicalP_xgt1_neg_mu_largetau_impl(const double mu, const double tau,
                                          const double x, double acosh_x,
                                          double * result, double * ln_multiplier)
{
  double xi = acosh_x;
  double ln_xi_pre;
  double ln_pre;
  double sumA, sumB, sum;
  double arg;
  double J_mup1, J_mu, J_mum1;

  if(xi < GSL_ROOT4_MACH_EPS) {
    ln_xi_pre = -xi*xi/6.0;           /* log(1.0 - xi*xi/6.0) */
  }
  else {
    double lnshxi;
    gsl_sf_lnsinh_impl(xi, &lnshxi);
    ln_xi_pre = log(xi) - lnshxi;     /* log(xi/sinh(xi) */
  }

  ln_pre = 0.5*ln_xi_pre - mu*log(tau);

  arg = tau*xi;
  gsl_sf_bessel_Jnu_impl(mu + 1.0,   arg, &J_mup1);
  gsl_sf_bessel_Jnu_impl(mu,         arg, &J_mu);
  J_mum1 = -J_mup1 + 2.0*mu/arg*J_mu;      /* careful of mu < 1 */

  sumA = 1.0 - olver_A1_xi(-mu, xi, x)/(tau*tau);
  sumB = olver_B0_xi(-mu, xi);
  sum  = J_mu * sumA - xi/tau * J_mum1 * sumB;

  if(sum == 0.0) {
    *result = 0.0;
    *ln_multiplier = 0.0;
    return GSL_SUCCESS;
  }
  else {
    int stat_e = gsl_sf_exp_mult_impl(ln_pre, sum, result);
    if(stat_e != GSL_SUCCESS) {
      *result = sum;
      *ln_multiplier = ln_pre;
    }
    else {
      *ln_multiplier = 0.0;
    }
    return GSL_SUCCESS;
  }
}


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}
 * -1 < x < 1
 * tau -> Inf 
 * [Olver, p. 473]
 */
int
gsl_sf_conicalP_xlt1_neg_mu_largetau_impl(const double mu, const double tau,
                                          const double x, const double acos_x,
                                          double * result, double * ln_multiplier)
{
  double theta = acos_x;
  double ln_th_pre;
  double ln_pre;
  double sumA, sumB, sum;
  double arg;
  double I_mup1, I_mu, I_mum1;

  if(theta < GSL_ROOT4_MACH_EPS) {
    ln_th_pre = theta*theta/6.0;   /* log(1.0 + theta*theta/6.0) */
  }
  else {
    
    ln_th_pre = log(theta/sin(theta));
  }

  ln_pre = 0.5 * ln_th_pre - mu * log(tau);

  arg = tau*theta;
  gsl_sf_bessel_Inu_impl(mu + 1.0,   arg, &I_mup1);
  gsl_sf_bessel_Inu_impl(mu,         arg, &I_mu);
  I_mum1 = I_mup1 + 2.0*mu/arg; /* careful of mu < 1 */

  sumA = 1.0 - olver_A1_th(-mu, theta, x)/(tau*tau);
  sumB = olver_B0_th(-mu, theta);
  sum  = I_mu * sumA - theta/tau * I_mum1 * sumB;

  if(sum == 0.0) {
    *result = 0.0;
    *ln_multiplier = 0.0;
    return GSL_SUCCESS;
  }
  else {
    int stat_e = gsl_sf_exp_mult_impl(ln_pre, sum, result);
    if(stat_e != GSL_SUCCESS) {
      *result = sum;
      *ln_multiplier = ln_pre;
    }
    else {
      *ln_multiplier = 0.0;
    }
    return GSL_SUCCESS;
  }
}


/* Hypergeometric function which appears in the
 * large x expansion below:
 *
 *   2F1(1/4 - mu/2 - I tau/2, 3/4 - mu/2 - I tau/2, 1 - I tau, y)
 *
 * Note that for the usage below y = 1/x^2;
 */
static
int
conicalP_hyperg_large_x(const double mu, const double tau, const double y,
                        double * reF, double * imF)
{
  const int kmax = 1000;
  const double re_a = 0.25 - 0.5*mu;
  const double re_b = 0.75 - 0.5*mu;
  const double re_c = 1.0;
  const double im_a = -0.5*tau;
  const double im_b = -0.5*tau;
  const double im_c = -tau;

  double re_sum = 1.0;
  double im_sum = 0.0;
  double re_term = 1.0;
  double im_term = 0.0;
  int k;

  for(k=1; k<=kmax; k++) {
    double re_ak = re_a + k - 1.0;
    double re_bk = re_b + k - 1.0;
    double re_ck = re_c + k - 1.0;
    double im_ak = im_a;
    double im_bk = im_b;
    double im_ck = im_c;
    double den = re_ck*re_ck + im_ck*im_ck;
    double re_multiplier = ((re_ak*re_bk - im_ak*im_bk)*re_ck + im_ck*(im_ak*re_bk + re_ak*im_bk)) / den;
    double im_multiplier = ((im_ak*re_bk + re_ak*im_bk)*re_ck - im_ck*(re_ak*re_bk - im_ak*im_bk)) / den;
    double re_tmp = re_multiplier*re_term - im_multiplier*im_term;
    double im_tmp = im_multiplier*re_term + re_multiplier*im_term;
    double asum = fabs(re_sum) + fabs(im_sum);
    re_term = y/k * re_tmp;
    im_term = y/k * im_tmp;
    if(fabs(re_term/asum) < GSL_MACH_EPS && fabs(im_term/asum) < GSL_MACH_EPS) break;
    re_sum += re_term;
    im_sum += im_term;
  }

  *reF = re_sum;
  *imF = im_sum;

  if(k == kmax)
    return GSL_EMAXITER;
  else  
    return GSL_SUCCESS;
}


/* P^{mu}_{-1/2 + I tau}
 * x->Inf
 */
int
gsl_sf_conicalP_large_x_impl(const double mu, const double tau, const double x,
                             double * result, double * ln_multiplier)
{
  /* 2F1 term
   */
  double y = ( x < 0.5*GSL_SQRT_DBL_MAX ? 1.0/(x*x) : 0.0 );
  double reF, imF;
  int stat_F = conicalP_hyperg_large_x(mu, tau, y, &reF, &imF);

  /* f = Gamma(+i tau)/Gamma(1/2 - mu + i tau)
   * FIXME: shift so it's better for tau-> 0
   */
  double lgr_num, lgth_num;
  double lgr_den, lgth_den;
  int stat_gn = gsl_sf_lngamma_complex_impl(0.0,tau,&lgr_num,&lgth_num);
  int stat_gd = gsl_sf_lngamma_complex_impl(0.5-mu,tau,&lgr_den,&lgth_den);

  double angle = lgth_num - lgth_den + atan2(imF,reF);

  double lnpre_const = 0.5*M_LN2 - 0.5*M_LNPI;
  double lnpre_comm = (mu-0.5)*log(x) - 0.5*mu*(log(x+1.0) + log(x-1.0));

  /*  result = pre*|F|*|f| * cos(angle - tau * (log(x)+M_LN2))
   */
  double c = cos(angle + tau*(log(x) + M_LN2));
  int status = GSL_ERROR_SELECT_3(stat_gd, stat_gn, stat_F);
  if(c == 0.0) {
    *result = 0.0;
    return status;
  }
  else {
    double lnFf     = 0.5*log(reF*reF+imF*imF) + lgr_num - lgr_den;
    double lnnoc    = lnpre_const + lnpre_comm + lnFf;
    int stat_e = gsl_sf_exp_mult_impl(lnnoc, c, result);
    if(stat_e == GSL_SUCCESS) {
      *ln_multiplier = 0.0;
    }
    else {
      *result = c;
      *ln_multiplier = lnnoc;
    }
    return status;
  }
}


/* P^{mu}_{-1/2 + I tau}  first hypergeometric representation
 * -1 < x < 1
 * This is more effective for |x| small, however it will work w/o
 * reservation for any x < 0 because everything is positive
 * definite in that case.
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

  int stat_F1 = gsl_sf_hyperg_2F1_conj_impl(0.25 - 0.5*mu, 0.5*tau, 0.5, x2, &F1);
  int stat_F2 = gsl_sf_hyperg_2F1_conj_impl(0.75 - 0.5*mu, 0.5*tau, 1.5, x2, &F2);
  int status = GSL_ERROR_SELECT_2(stat_F1, stat_F2);

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
#if 0
static
int
conicalP_def_hyperg(double mu, double tau, double x, double * result)
{
  double F;
  int stat_F = gsl_sf_hyperg_2F1_conj_renorm_impl(0.5, tau, 1.0-mu, 0.5*(1.0-x), &F);
  *result = pow((x+1.0)/(x-1.0), 0.5*mu) * F;
  return stat_F;
}
#endif /* 0 */


/* P^{mu}_{-1/2 + I tau}  second hypergeometric representation
 * [Zhurina+Karmazina, (3.1)] 
 * -1 < x < 3
 * effective for x near 1
 *
 */
#if 0
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
#endif /* 0 */


/* V0, V1 from Kolbig, m = 0
 */
static
int
conicalP_0_V(const double t, const double f, const double tau, const double sgn,
             double * V0, double * V1)
{
  double C[8];
  double T[8];
  double H[8];
  double V[12];
  int i;
  T[0] = 1.0;
  H[0] = 1.0;
  V[0] = 1.0;
  for(i=1; i<=7; i++) {
    T[i] = T[i-1] * t;
    H[i] = H[i-1] * (t*f);
  }
  for(i=1; i<=11; i++) {
    V[i] = V[i-1] * tau;
  }

  C[0] = 1.0;
  C[1] = (H[1]-1.0)/(8.0*T[1]);
  C[2] = (9.0*H[2] + 6.0*H[1] - 15.0 - sgn*8.0*T[2])/(128.0*T[2]);
  C[3] = 5.0*(15.0*H[3] + 27.0*H[2] + 21.0*H[1] - 63.0 - sgn*T[2]*(16.0*H[1]+24.0))/(1024.0*T[3]);
  C[4] = 7.0*(525.0*H[4] + 1500.0*H[3] + 2430.0*H[2] + 1980.0*H[1] - 6435.0
              + 192.0*T[4] - sgn*T[2]*(720.0*H[2]+1600.0*H[1]+2160.0)
              ) / (32768.0*T[4]);
  C[5] = 21.0*(2835.0*H[5] + 11025.0*H[4] + 24750.0*H[3] + 38610.0*H[2]
               + 32175.0*H[1] - 109395.0 + T[4]*(1984.0*H[1]+4032.0)
               - sgn*T[2]*(4800.0*H[3]+15120.0*H[2]+26400.0*H[1]+34320.0)
	       ) / (262144.0*T[5]);
  C[6] = 11.0*(218295.0*H[6] + 1071630.0*H[5] + 3009825.0*H[4] + 6142500.0*H[3]
               + 9398025.0*H[2] + 7936110.0*H[1] - 27776385.0
	       + T[4]*(254016.0*H[2]+749952.0*H[1]+1100736.0)
	       - sgn*T[2]*(441000.0*H[4] + 1814400.0*H[3] + 4127760.0*H[2]
	                 + 6552000.0*H[1] + 8353800.0 + 31232.0*T[4]
			 )
               ) / (4194304.0*T[6]);

  *V0 = C[0] + (-4.0*C[3]/T[1]+C[4])/V[4]
             + (-192.0*C[5]/T[3]+144.0*C[6]/T[2])/V[8]
             + sgn * (-C[2]/V[2]
	              + (-24.0*C[4]/T[2]+12.0*C[5]/T[1]-C[6])/V[6] 
	              + (-1920.0*C[6]/T[4])/V[10]
		      );
  *V1 = C[1]/V[1] + (8.0*(C[3]/T[2]-C[4]/T[1])+C[5])/V[5]
                  + (384.0*C[5]/T[4] - 768.0*C[6]/T[3])/V[9]
                  + sgn * ((2.0*C[2]/T[1]-C[3])/V[3]
		           + (48.0*C[4]/T[3]-72.0*C[5]/T[2] + 18.0*C[6]/T[1])/V[7]
		           + (3840.0*C[6]/T[5])/V[11]
		           );

  return GSL_SUCCESS;
}


/* V0, V1 from Kolbig, m = 1
 */
static
int
conicalP_1_V(const double t, const double f, const double tau, const double sgn,
             double * V0, double * V1)
{
  double Cm1;
  double C[8];
  double T[8];
  double H[8];
  double V[12];
  int i;
  T[0] = 1.0;
  H[0] = 1.0;
  V[0] = 1.0;
  for(i=1; i<=7; i++) {
    T[i] = T[i-1] * t;
    H[i] = H[i-1] * (t*f);
  }
  for(i=1; i<=11; i++) {
    V[i] = V[i-1] * tau;
  }

  Cm1  = -1.0;
  C[0] = 3.0*(1.0-H[1])/(8.0*T[1]);
  C[1] = (-15.0*H[2]+6.0*H[1]+9.0+sgn*8.0*T[2])/(128.0*T[2]);
  C[2] = 3.0*(-35.0*H[3] - 15.0*H[2] + 15.0*H[1] + 35.0 + sgn*T[2]*(32.0*H[1]+8.0))/(1024.0*T[3]);
  C[3] = (-4725.0*H[4] - 6300.0*H[3] - 3150.0*H[2] + 3780.0*H[1] + 10395.0
          -1216.0*T[4] + sgn*T[2]*(6000.0*H[2]+5760.0*H[1]+1680.0)) / (32768.0*T[4]);
  C[4] = 7.0*(-10395.0*H[5] - 23625.0*H[4] - 28350.0*H[3] - 14850.0*H[2]
              +19305.0*H[1] + 57915.0 - T[4]*(6336.0*H[1]+6080.0)
	      + sgn*T[2]*(16800.0*H[3] + 30000.0*H[2] + 25920.0*H[1] + 7920.0)
	      ) / (262144.0*T[5]);
  C[5] = (-2837835.0*H[6] - 9168390.0*H[5] - 16372125.0*H[4] - 18918900*H[3]
          -10135125.0*H[2] + 13783770.0*H[1] + 43648605.0
	  -T[4]*(3044160.0*H[2] + 5588352.0*H[1] + 4213440.0)
	  +sgn*T[2]*(5556600.0*H[4] + 14817600.0*H[3] + 20790000.0*H[2]
	             + 17297280.0*H[1] + 5405400.0 + 323072.0*T[4]
		     )
          ) / (4194304.0*T[6]);
  C[6] = 0.0;

  *V0 = C[0] + (-4.0*C[3]/T[1]+C[4])/V[4]
             + (-192.0*C[5]/T[3]+144.0*C[6]/T[2])/V[8]
             + sgn * (-C[2]/V[2]
	              + (-24.0*C[4]/T[2]+12.0*C[5]/T[1]-C[6])/V[6] 
		      );
  *V1 = C[1]/V[1] + (8.0*(C[3]/T[2]-C[4]/T[1])+C[5])/V[5]
                  + (384.0*C[5]/T[4] - 768.0*C[6]/T[3])/V[9]
                  + sgn * (Cm1*V[1] + (2.0*C[2]/T[1]-C[3])/V[3]
		           + (48.0*C[4]/T[3]-72.0*C[5]/T[2] + 18.0*C[6]/T[1])/V[7]
		           );

  return GSL_SUCCESS;
}



/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* P^0_{-1/2 + I lambda}
 */
int
gsl_sf_conicalP_0_impl(const double lambda, const double x, double * result)
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x == 1.0) {
    *result = 1.0;
    return GSL_SUCCESS;
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
  else if(   (x <= 0.0 && lambda < 1000.0)
          || (x <  0.1 && lambda < 17.0)
	  || (x <  0.2 && lambda < 5.0 )
    ) {
    return conicalP_xlt1_hyperg_A(0.0, lambda, x, result);
  }
  else if(   (x <= 0.2 && lambda < 17.0)
          || (x <= 1.5 && lambda < 20.0)
    ) {
    return gsl_sf_hyperg_2F1_conj_impl(0.5, lambda, 1.0, (1.0-x)/2, result);
  }
  else if(1.5 < x && lambda < locMAX(x,20.0)) {
    double P;
    double lm;
    int stat_P = gsl_sf_conicalP_large_x_impl(0.0, lambda, x,
                                              &P, &lm
                                              );
    if(P != 0.0 && (stat_P == GSL_SUCCESS || stat_P == GSL_ELOSS)) {
      int stat_e = gsl_sf_exp_mult_impl(lm, P, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_P);
    }
    else {
      *result = 0.0;
      return stat_P;
    }
  }
  else {
    double V0, V1;
    if(x < 1.0) {
      double th  = acos(x);
      double sth = sqrt(1.0-x*x);  /* sin(th) */
      double I0, I1;
      int stat_I0 = gsl_sf_bessel_I0_scaled_impl(th * lambda, &I0);
      int stat_I1 = gsl_sf_bessel_I1_scaled_impl(th * lambda, &I1);
      int stat_I  = GSL_ERROR_SELECT_2(stat_I0, stat_I1);
      int stat_V  = conicalP_0_V(th, x/sth, lambda, -1.0, &V0, &V1);
      double bessterm = V0 * I0 + V1 * I1;
      int stat_e = gsl_sf_exp_mult_impl(th * lambda, sqrt(th/sth) * bessterm, result);
      return GSL_ERROR_SELECT_3(stat_e, stat_V, stat_I);
    }
    else {
      double sh = sqrt(x-1.0)*sqrt(x+1.0);  /* sinh(xi)      */
      double xi = log(x + sh);              /* xi = acosh(x) */
      double J0, J1;
      int stat_J0 = gsl_sf_bessel_J0_impl(xi * lambda, &J0);
      int stat_J1 = gsl_sf_bessel_J1_impl(xi * lambda, &J1);
      int stat_J  = GSL_ERROR_SELECT_2(stat_J0, stat_J1);
      int stat_V  = conicalP_0_V(xi, x/sh, lambda, 1.0, &V0, &V1);
      double bessterm = V0 * J0 + V1 * J1;
      *result = sqrt(xi/sh) * bessterm;
      return GSL_ERROR_SELECT_2(stat_V, stat_J);
    }
  }
}


/* P^1_{-1/2 + I lambda}
 */
int
gsl_sf_conicalP_1_impl(const double lambda, const double x, double * result)
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(lambda == 0.0) {
    double K, E;
    int stat_K, stat_E;
    if(x == 1.0) {
      *result = 0.0;
      return GSL_SUCCESS;
    }
    else if(x < 1.0) {
      if(1.0-x < GSL_SQRT_MACH_EPS) {
        *result = 0.25/M_SQRT2 * sqrt(1.0-x) * (1.0 + 5.0/16.0 * (1.0-x));
	return GSL_SUCCESS;
      }
      else {
        const double th = acos(x);
        const double s  = sin(0.5*th);
        const double c2 = 1.0 - s*s;
        const double sth = sin(th);
        stat_K = gsl_sf_ellint_Kcomp_impl(s, locEPS, &K);
        stat_E = gsl_sf_ellint_Ecomp_impl(s, locEPS, &E);
        *result = 2.0/(M_PI*sth) * (E - c2 * K);
        return stat_K;
      }
    }
    else {
      if(x-1.0 < GSL_SQRT_MACH_EPS) {
        *result = -0.25/M_SQRT2 * sqrt(x-1.0) * (1.0 - 5.0/16.0 * (x-1.0));
	return GSL_SUCCESS;
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
  }
  else if(   (x <= 0.0 && lambda < 1000.0)
          || (x <  0.1 && lambda < 17.0)
	  || (x <  0.2 && lambda < 5.0 )
    ) {
    return conicalP_xlt1_hyperg_A(1.0, lambda, x, result);
  }
  else if(   (x <= 0.2 && lambda < 17.0)
          || (x <  1.5 && lambda < 20.0)
    ) {
    const double arg = fabs(x*x - 1.0);
    const double sgn = GSL_SIGN(1.0 - x);
    const double pre = 0.5*(lambda*lambda + 0.25) * sgn * sqrt(arg);
    double F;
    int stat_F = gsl_sf_hyperg_2F1_conj_impl(1.5, lambda, 2.0, (1.0-x)/2, &F);
    *result = pre * F;
    return stat_F;
  }
  else if(1.5 <= x && lambda < locMAX(x,20.0)) {
    double P;
    double lm;
    int stat_P = gsl_sf_conicalP_large_x_impl(1.0, lambda, x,
                                              &P, &lm
                                              );
    if(P != 0.0 && (stat_P == GSL_SUCCESS || stat_P == GSL_ELOSS)) {
      int stat_e = gsl_sf_exp_mult_impl(lm, P, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_P);
    }
    else {
      *result = 0.0;
      return stat_P;
    }
  }
  else {
    double V0, V1;
    if(x < 1.0) {
      double th  = acos(x);
      double sth = sqrt(1.0-x*x);  /* sin(th) */
      double I0, I1;
      int stat_I0 = gsl_sf_bessel_I0_scaled_impl(th * lambda, &I0);
      int stat_I1 = gsl_sf_bessel_I1_scaled_impl(th * lambda, &I1);
      int stat_I  = GSL_ERROR_SELECT_2(stat_I0, stat_I1);
      int stat_V  = conicalP_1_V(th, x/sth, lambda, -1.0, &V0, &V1);
      double bessterm = V0 * I0 + V1 * I1;
      int stat_e = gsl_sf_exp_mult_impl(th * lambda, sqrt(th/sth) * bessterm, result);
      return GSL_ERROR_SELECT_3(stat_e, stat_V, stat_I);
    }
    else {
      double sh = sqrt(x-1.0)*sqrt(x+1.0);  /* sinh(xi)      */
      double xi = log(x + sh);              /* xi = acosh(x) */
      double J0, J1;
      int stat_J0 = gsl_sf_bessel_J0_impl(xi * lambda, &J0);
      int stat_J1 = gsl_sf_bessel_J1_impl(xi * lambda, &J1);
      int stat_J  = GSL_ERROR_SELECT_2(stat_J0, stat_J1);
      int stat_V  = conicalP_1_V(xi, x/sh, lambda, 1.0, &V0, &V1);
      double bessterm = V0 * J0 + V1 * J1;
      *result = sqrt(xi/sh) * bessterm;
      return GSL_ERROR_SELECT_2(stat_V, stat_J);
    }
  }
}


/* P^{1/2}_{-1/2 + I lambda} (x)
 * [Abramowitz+Stegun 8.6.8, 8.6.12]
 * checked OK [GJ] Fri May  8 12:24:36 MDT 1998 
 */
int gsl_sf_conicalP_half_impl(const double lambda, const double x,
                              double * result
                              )
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0) {
    double ac  = acos(x);
    double den = sqrt(sqrt(1.0-x)*sqrt(1.0+x));
    *result = Root_2OverPi_ / den * cosh(ac * lambda);
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* x > 1 */
    double sq_term = sqrt(x-1.0)*sqrt(x+1.0);
    double ln_term = log(x + sq_term);
    double den = sqrt(sq_term);
    *result = Root_2OverPi_ / den * cos(lambda * ln_term);
    return GSL_SUCCESS;
  }
}


/* P^{-1/2}_{-1/2 + I lambda} (x)
 * [Abramowitz+Stegun 8.6.9, 8.6.14]
 * checked OK [GJ] Fri May  8 12:24:43 MDT 1998 
 */
int gsl_sf_conicalP_mhalf_impl(const double lambda, const double x, double * result)
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0) {
    double ac  = acos(x);
    double den = sqrt(sqrt(1.0-x)*sqrt(1.0+x));
    double arg = ac * lambda;
    if(fabs(arg) < GSL_SQRT_MACH_EPS) {
      *result = Root_2OverPi_ / den * ac;
    }
    else {
      *result = Root_2OverPi_ / (den*lambda) * sinh(arg);
    }
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* x > 1 */
    double sq_term = sqrt(x-1.0)*sqrt(x+1.0);
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
}


int gsl_sf_conicalP_sph_reg_impl(const int l, const double lambda,
                                 const double x,
                                 double * result
                                 )
{
  if(x <= -1.0 || l < -1) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(l == -1) {
    return gsl_sf_conicalP_half_impl(lambda, x, result);
  }
  else if(l == 0) {
    return gsl_sf_conicalP_mhalf_impl(lambda, x, result);
  }
  else if(x == 1.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 0.0) {
    double c = 1.0/sqrt(1.0-x*x);
    double Pellm1;
    double Pell;
    double Pellp1;
    int ell;
    int stat_0 = gsl_sf_conicalP_half_impl(lambda, x, &Pellm1);  /* P^( 1/2) */
    int stat_1 = gsl_sf_conicalP_mhalf_impl(lambda, x, &Pell);   /* P^(-1/2) */
    int stat_P = GSL_ERROR_SELECT_2(stat_0, stat_1);

    for(ell=0; ell<l; ell++) {
      double d = (ell+1.0)*(ell+1.0) + lambda*lambda;
      Pellp1 = (Pellm1 - (2.0*ell+1.0)*c*x * Pell) / d;
      Pellm1 = Pell;
      Pell   = Pellp1;
    }

    *result = Pell;
    return stat_P;
  }
  else if(x < 1.0) {
    const double xi = x/(sqrt(1.0-x)*sqrt(1.0+x));
    double rat;
    double Phf;
    int stat_CF1 = conicalP_negmu_xlt1_CF1(0.5, l, lambda, x, &rat);
    int stat_Phf = gsl_sf_conicalP_half_impl(lambda, x, &Phf);
    double Pellp1 = rat * GSL_SQRT_DBL_MIN;
    double Pell   = GSL_SQRT_DBL_MIN;
    double Pellm1;
    int ell;

    for(ell=l; ell>=0; ell--) {
      double d = (ell+1.0)*(ell+1.0) + lambda*lambda;
      Pellm1 = (2.0*ell+1.0)*xi * Pell + d * Pellp1;
      Pellp1 = Pell;
      Pell   = Pellm1;
    }

    *result = GSL_SQRT_DBL_MIN * Phf / Pell;
    return GSL_ERROR_SELECT_2(stat_Phf, stat_CF1);
  }
  else if(x == 1.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* x > 1.0 */

    const double xi = x/sqrt((x-1.0)*(x+1.0));
    double rat;
    int stat_CF1 = conicalP_negmu_xgt1_CF1(0.5, l, lambda, x, &rat);
    int stat_P;
    double Pellp1 = rat * GSL_SQRT_DBL_MIN;
    double Pell   = GSL_SQRT_DBL_MIN;
    double Pellm1;
    int ell;

    for(ell=l; ell>=0; ell--) {
      double d = (ell+1.0)*(ell+1.0) + lambda*lambda;
      Pellm1 = (2.0*ell+1.0)*xi * Pell - d * Pellp1;
      Pellp1 = Pell;
      Pell   = Pellm1;
    }

    if(fabs(Pell) > fabs(Pellp1)){
      double Phf;
      stat_P = gsl_sf_conicalP_half_impl(lambda, x, &Phf);
      *result = GSL_SQRT_DBL_MIN * Phf / Pell;
    }
    else {
      double Pmhf;
      stat_P = gsl_sf_conicalP_mhalf_impl(lambda, x, &Pmhf);
      *result = GSL_SQRT_DBL_MIN * Pmhf / Pellp1;
    }

    return GSL_ERROR_SELECT_2(stat_P, stat_CF1);
  }
}


int gsl_sf_conicalP_cyl_reg_impl(const int m, const double lambda,
                                 const double x,
                                 double * result
                                 )
{
  if(x <= -1.0 || m < -1) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(m == -1) {
    return gsl_sf_conicalP_1_impl(lambda, x, result);
  }
  else if(m == 0) {
    return gsl_sf_conicalP_0_impl(lambda, x, result);
  }
  else if(x == 1.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 0.0) {
    double c = 1.0/sqrt(1.0-x*x);
    double Pkm1;
    double Pk;
    double Pkp1;
    int k;
    int stat_0 = gsl_sf_conicalP_1_impl(lambda, x, &Pkm1);  /* P^1 */
    int stat_1 = gsl_sf_conicalP_0_impl(lambda, x, &Pk);    /* P^0 */
    int stat_P = GSL_ERROR_SELECT_2(stat_0, stat_1);

    for(k=0; k<m; k++) {
      double d = (k+0.5)*(k+0.5) + lambda*lambda;
      Pkp1 = (Pkm1 - 2.0*k*c*x * Pk) / d;
      Pkm1 = Pk;
      Pk   = Pkp1;
    }

    *result = Pk;
    return stat_P;
  }
  else if(x < 1.0) {
    const double xi = x/(sqrt(1.0-x)*sqrt(1.0+x));
    double rat;
    double P0;
    int stat_CF1 = conicalP_negmu_xlt1_CF1(0.0, m, lambda, x, &rat);
    int stat_P0  = gsl_sf_conicalP_0_impl(lambda, x, &P0);
    double Pkp1 = rat * GSL_SQRT_DBL_MIN;
    double Pk   = GSL_SQRT_DBL_MIN;
    double Pkm1;
    int k;

    for(k=m; k>0; k--) {
      double d = (k+0.5)*(k+0.5) + lambda*lambda;
      Pkm1 = 2.0*k*xi * Pk + d * Pkp1;
      Pkp1 = Pk;
      Pk   = Pkm1;
    }

    *result = GSL_SQRT_DBL_MIN * P0 / Pk;
    return GSL_ERROR_SELECT_2(stat_P0, stat_CF1);
  }
  else if(x == 1.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* x > 1.0 */

    const double xi = x/sqrt((x-1.0)*(x+1.0));
    double rat;
    int stat_CF1 = conicalP_negmu_xgt1_CF1(0.0, m, lambda, x, &rat);
    int stat_P;
    double Pkp1 = rat * GSL_SQRT_DBL_MIN;
    double Pk   = GSL_SQRT_DBL_MIN;
    double Pkm1;
    int k;

    for(k=m; k>-1; k--) {
      double d = (k+0.5)*(k+0.5) + lambda*lambda;
      Pkm1 = 2.0*k*xi * Pk - d * Pkp1;
      Pkp1 = Pk;
      Pk   = Pkm1;
    }

    if(fabs(Pk) > fabs(Pkp1)){
      double P1;
      stat_P = gsl_sf_conicalP_1_impl(lambda, x, &P1);
      *result = GSL_SQRT_DBL_MIN * P1 / Pk;
    }
    else {
      double P0;
      stat_P = gsl_sf_conicalP_0_impl(lambda, x, &P0);
      *result = GSL_SQRT_DBL_MIN * P0 / Pkp1;
    }

    return GSL_ERROR_SELECT_2(stat_P, stat_CF1);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_conicalP_0_e(const double lambda, const double x, double * result)
{
  int status = gsl_sf_conicalP_0_impl(lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conicalP_0_e", status);
  }
  return status;
}

int gsl_sf_conicalP_1_e(const double lambda, const double x, double * result)
{
  int status = gsl_sf_conicalP_1_impl(lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conicalP_1_e", status);
  }
  return status;
}

int gsl_sf_conicalP_half_e(const double lambda, const double x, double * result)
{
  int status = gsl_sf_conicalP_half_impl(lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conicalP_half_e", status);
  }
  return status;
}

int gsl_sf_conicalP_mhalf_e(const double lambda, const double x, double * result)
{
  int status = gsl_sf_conicalP_mhalf_impl(lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conicalP_mhalf_e", status);
  }
  return status;
}

int gsl_sf_conicalP_sph_reg_e(const int l, const double lambda, const double x, double * result)
{
  int status = gsl_sf_conicalP_sph_reg_impl(l, lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conicalP_sph_reg_e", status);
  }
  return status;
}

int gsl_sf_conicalP_cyl_reg_e(const int m, const double lambda, const double x, double * result)
{
  int status = gsl_sf_conicalP_cyl_reg_impl(m, lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conicalP_cyl_reg_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/


double gsl_sf_conicalP_0(const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conicalP_0_impl(lambda, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conicalP_0", status);
  }
  return y;
}

double gsl_sf_conicalP_1(const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conicalP_1_impl(lambda, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conicalP_1", status);
  }
  return y;
}

double gsl_sf_conicalP_half(const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conicalP_half_impl(lambda, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conicalP_half", status);
  }
  return y;
}

double gsl_sf_conicalP_mhalf(const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conicalP_mhalf_impl(lambda, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conicalP_mhalf", status);
  }
  return y;
}

double gsl_sf_conicalP_sph_reg(const int l, const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conicalP_sph_reg_impl(l, lambda, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conicalP_sph_reg", status);
  }
  return y;
}


double gsl_sf_conicalP_cyl_reg(const int m, const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conicalP_cyl_reg_impl(m, lambda, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conicalP_cyl_reg", status);
  }
  return y;
}
