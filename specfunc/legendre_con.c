/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_poly.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_bessel.h"
#include "gsl_sf_legendre.h"

#define Root_2OverPi_  0.797884560802865355879892


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326 */

static inline double olver_U1(double beta2, double p)
{
  return (p-1.)/(24.*(1.+beta2)) * (3. + beta2*(2. + 5.*p*(1.+p)));
}
static inline double olver_U2(double beta2, double p)
{
  double beta4 = beta2*beta2;
  double p2    = p*p;
  double poly1 =  4.*beta4 + 84.*beta2 - 63.;
  double poly2 = 16.*beta4 + 90.*beta2 - 81.;
  double poly3 = beta2*p2*(97.*beta2 - 432. + 77.*p*(beta2-6.) - 385.*beta2*p2*(1. + p));
  return (1.-p)/(1152.*(1.+beta2)) * (poly1 + poly2 + poly3);
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
  double den   = 39813120. * opb2s*opb2s;
  double poly1 = gsl_sf_poly_eval(U3c1, 4, p);
  double poly2 = gsl_sf_poly_eval(U3c2, 6, p);
  double poly3 = gsl_sf_poly_eval(U3c3, 8, p);
  double poly4 = gsl_sf_poly_eval(U3c4, 10, p);
  double poly5 = gsl_sf_poly_eval(U3c5, 12, p);
  
  return (p-1.)*(      1215*poly1 + 324*beta2*poly2
                 + 54*beta4*poly3 +  12*beta6*poly4
		 + beta4*beta4*poly5
		 ) / den;
}

/* P^{-mu}_{-1/2 + I tau}, mu -> Inf */
static int conicalP_xlt1_large_neg_mu(double mu, double tau, double x, double * result)
{
  double beta  = tau/mu;
  double beta2 = beta*beta;
  double S     = beta * acos((1.-beta2)/(1.+beta2));
  double p     = x/sqrt(beta2*(1.-x*x) + 1.);
  double ln_pre_1 =  0.5*mu*(S - log(1.+beta2) + log((1.-p)/(1.+p))) - gsl_sf_lngamma(mu+1.);
  double ln_pre_2 = -0.25 * log(1. + beta2*(1.-x));
  double ln_pre_3 = -tau * atan(p*beta);
  double ln_pre = ln_pre_1 + ln_pre_2 + ln_pre_3;
  
  if(ln_pre > GSL_LOG_DBL_MAX) {
    *result = 0.; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(ln_pre < GSL_LOG_DBL_MIN) {
    *result = 0.;
    return GSL_EUNDRFLW;
  }
  else {
    double sum = 1. - olver_U1(beta2, p)/mu + olver_U2(beta2, p)/(mu*mu);
    *result = exp(ln_pre) * sum;
    return GSL_SUCCESS;
  }
}


/* A_n^{-mu}, B_n^{-mu}   [Olver, p.465, 469] */

static double olver_B0_xi(double mu, double xi)
{
  return (1. - 4.*mu*mu)/(8.*xi) * (coth(xi) - 1./xi);
}
static double olver_A1_xi(double mu, double xi, double x)
{
  double B = olver_B0_xi(mu, xi);
  double psi = (4.*mu*mu - 1.)/16. * (1./(x*x-1.) - 1./(xi*xi));
  return 0.5*xi*xi*B*B + (mu+0.5)*B - psi + mu/6.*(0.25 - mu*mu);
}
static double olver_B0_th(double mu, double theta)
{
  return -(1. - 4.*mu*mu)/(8.*theta) * (cot(theta) - 1./theta);
}
static double olver_A1_th(double mu, double theta, double x)
{
  double B = olver_B0_th(mu, theta);
  double psi = (4.*mu*mu - 1.)/16. * (1./(x*x-1.) + 1./(theta*theta));
  return -0.5*theta*theta*B*B + (mu+0.5)*B - psi + mu/6.*(0.25 - mu*mu);
}


/* P^{-m}_{-1/2 + I tau}, tau -> Inf   [Olver, p. 473] */
static int cylconicalP_xlt1_large_tau_impl(int m, double tau, double x, double * result)
{
  double theta = acos(x);
  double th_pre;
  double pre;
  double sumA, sumB;
  double arg;
  double Im, Imm1;
  
  if(theta < GSL_ROOT4_MACH_EPS) {
    th_pre = 1. + theta*theta/6.;
  }
  else {
    th_pre = theta/sin(theta);
  }
  
  pre = sqrt(th_pre) * gsl_sf_pow_int(1./tau, m);
  if(pre == 0.0) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  
  arg = tau*theta;
  
  gsl_sf_bessel_In_impl(m,   arg, &Im);
  gsl_sf_bessel_In_impl(m-1, arg, &Imm1);
  
  sumA = 1. - olver_A1_th(m, theta, x)/(tau*tau);
  sumB = olver_B0_th(m, theta)/tau;
  
  *result = pre * (Im * sumA - theta/tau * Imm1 * sumB);
  return GSL_SUCCESS; /* FIXME: hmmm, success??? */
}

/* P^{-1/2 - ell}_{-1/2 + I tau}, tau -> Inf   [Olver, p. 473]
 * 
 */
static int sphconicalP_xlt1_large_tau_impl(int ell, double tau, double x, double * result)
{
  double mu = ell + 0.5;
  double theta = acos(x);
  double th_pre;
  double pre;
  double sumA, sumB;
  double arg;
  double rt_term;
  double Im, Imm1;
  
  if(theta < GSL_ROOT4_MACH_EPS) {
    th_pre = 1. + theta*theta/6.;
  }
  else {
    th_pre = theta/sin(theta);
  }
  
  pre = sqrt(th_pre) * pow(1./tau, mu);
  if(pre == 0.0) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  
  arg = tau*theta;

  gsl_sf_bessel_il_impl(ell,   arg, &Im);
  gsl_sf_bessel_il_impl(ell-1, arg, &Imm1);
  rt_term = sqrt(2.0*arg/M_PI);
  Im   *= rt_term;
  Imm1 *= rt_term;

  sumA = 1. - olver_A1_th(mu, theta, x)/(tau*tau);
  sumB = olver_B0_th(mu, theta)/tau;
  
  *result = pre * (Im * sumA - theta/tau * Imm1 * sumB);
  return GSL_SUCCESS; /* FIXME: hmmm, success??? */
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* P^{1/2}_{-1/2 + I lambda} (x)
 * [Abramowitz+Stegun 8.6.8, 8.6.12]
 * checked OK [GJ] Fri May  8 12:24:36 MDT 1998 
 */
int gsl_sf_conical_sph_irr_1_impl(const double lambda,
                                  const double one_minus_x,
                                  const double one_plus_x,
                                  double * result
                                  )
{
  double x = 0.5 * (one_plus_x - one_minus_x);
  if(fabs(x) < 1.) {
    double ac  = acos(x);
    double den = sqrt(sqrt(one_minus_x*one_plus_x));
    *result = Root_2OverPi_ / den * cosh(ac * lambda);
    return GSL_SUCCESS;
  }
  else if(one_minus_x < 0.) { /* x > 1. */
    double sq_term = sqrt(-one_minus_x*one_plus_x);
    double ln_term = log(x + sq_term);
    double den = sqrt(sq_term);
    *result = Root_2OverPi_ / den * cos(lambda * ln_term);
    return GSL_SUCCESS;
  }
  else {
    /* |x| == 1 or x < -1 */
    return GSL_EDOM;
  }
}


/* P^{-1/2}_{-1/2 + I lambda} (x)
 * [Abramowitz+Stegun 8.6.9, 8.6.14]
 * checked OK [GJ] Fri May  8 12:24:43 MDT 1998 
 */
int gsl_sf_conical_sph_reg_0_impl(const double lambda,
                                  const double one_minus_x,
                                  const double one_plus_x,
                                  double * result
                                  )
{
  double x = 0.5 * (one_plus_x - one_minus_x);
  if(fabs(x) < 1.) {
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
  else if(one_minus_x < 0.) { /* x > 1. */
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
  else if(one_minus_x == 0.) { /* x == 1. */
    *result = 0.;
    return GSL_SUCCESS;
  }  
  else {
    /* x <= -1 */
    return GSL_EDOM;
  }
}

int gsl_sf_conical_sph_reg_impl(const int lmax, const double lambda,
                                const double one_m_x, const double one_p_x,
                                double * result, double * harvest
				)
{
  double x = 0.5 * (one_p_x - one_m_x);

  if(fabs(x) < 1.) {
    double f0;
    double p[2];
    int l_start = 20 + (int) ceil(lmax * (1. + (0.14 + 0.026*lambda)));
    p[0] = lambda*lambda;
    p[1] = x/sqrt(one_m_x*one_p_x);
    gsl_sf_conical_sph_reg_0_impl(lambda, one_m_x, one_p_x, &f0);  /* l =  0  */ 
    recurse_backward_minimal_simple_conical_sph_reg_xlt1(l_start, lmax, 0, p, f0, harvest, result);
    
  }
  else if(x > 1.) {
    double f0;
    double p[2];
    int l_start = 10 + (int) ceil(lmax * (1. + (0.14 + 0.026*lambda)*(x-1.)));
    p[0] = lambda*lambda;
    p[1] = x/sqrt(-one_m_x*one_p_x);
    gsl_sf_conical_sph_reg_0_impl(lambda, one_m_x, one_p_x, &f0);
    recurse_backward_minimal_simple_conical_sph_reg_xgt1(l_start, lmax, 0, p, f0, harvest, result);
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

int gsl_sf_conical_sph_irr_1_e(const double lambda, const double x, double * result)
{
  int status = gsl_sf_conical_sph_irr_1_impl(lambda, 1.-x, 1.+x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_sph_irr_1_e", status);
  }
  return status;
}

int gsl_sf_conical_sph_reg_0_e(const double lambda, const double x, double * result)
{
  int status = gsl_sf_conical_sph_reg_0_impl(lambda, 1.-x, 1.+x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_sph_reg_0_e", status);
  }
  return status;
}

int gsl_sf_conical_sph_reg_e(const int l, const double lambda, const double x, double * result)
{
  int status = gsl_sf_conical_sph_reg_impl(l, lambda, 1.-x, 1.+x, result, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_sph_reg_e", status);
  }
  return status;
}

int gsl_sf_conical_sph_reg_array_e(const int l, const double lambda, const double x, double * result_array)
{
  double y;
  int status = gsl_sf_conical_sph_reg_impl(l, lambda, 1.-x, 1.+x, &y, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_conical_sph_reg_array_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_conical_sph_irr_1(const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conical_sph_irr_1_impl(lambda, 1.-x, 1.+x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conical_sph_irr_1_e", status);
  }
  return y;
}

double gsl_sf_conical_sph_reg_0(const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conical_sph_reg_0_impl(lambda, 1.-x, 1.+x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conical_sph_reg_0_e", status);
  }
  return y;
}

double gsl_sf_conical_sph_reg(const int l, const double lambda, const double x)
{
  double y;
  int status = gsl_sf_conical_sph_reg_impl(l, lambda, 1.-x, 1.+x, &y, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_conical_sph_reg", status);
  }
  return y;
}
