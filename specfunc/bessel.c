/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Miscellaneous support functions for Bessel function evaluations.
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_airy.h"
#include "gsl_sf_elementary.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_trig.h"
#include "bessel_amp_phase.h"
#include "bessel.h"

#define CubeRoot2_  1.25992104989487316476721060728


/* sum which occurs in Taylor series for J_nu(x) * Gamma(nu+1)/((x/2)^nu)
 * or for I_nu(x) * Gamma(nu+1)/((x/2)^nu)
 *  sign = -1  ==> Jnu
 *  sign = +1  ==> Inu
 * [Abramowitz+Stegun, 9.1.10]
 * [Abramowitz+Stegun, 9.6.7]
 *
 * Assumes: nu >= 0
 *
 * checked OK [GJ] Mon May  4 00:11:54 EDT 1998 
 */
static int Inu_Jnu_taylorsum(const double nu, const double x,
                             const int sign,
                             const int kmax,
			     const double threshold,
                             double * result
                             )
{
  int k;
  double y = sign * 0.25 * x*x;

  double sum  = 1.0;
  double term = 1.0;

  for(k=1; k<=kmax; k++) {
    term *= y/(nu+k)/k;
    sum  += term;
    if(fabs(term/sum) < threshold) break;
  }

  *result = sum;

  if(k == kmax)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* Debye functions [Abramowitz+Stegun, 9.3.9-10] */

#ifdef HAVE_INLINE
inline
#endif
static double debye_u1(const double * tpow)
{
  return (3.0*tpow[1] - 5.0*tpow[3])/24.0;
}
#ifdef HAVE_INLINE
inline
#endif
static double debye_u2(const double * tpow)
{
  return (81.0*tpow[2] - 462.0*tpow[4] + 385.0*tpow[6])/1152.0;
}
#ifdef HAVE_INLINE
inline
#endif
static double debye_u3(const double * tpow)
{
  return (30375.0*tpow[3] - 369603.0*tpow[5] + 765765.0*tpow[7] - 425425.0*tpow[9])/414720.0;
}
#ifdef HAVE_INLINE
inline
#endif
static double debye_u4(const double * tpow)
{
  return (4465125.0*tpow[4] - 94121676.0*tpow[6] + 349922430.0*tpow[8] - 
          446185740.0*tpow[10] + 185910725.0*tpow[12])/39813120.0;
}
#ifdef HAVE_INLINE
inline
#endif
static double debye_u5(const double * tpow)
{
  return (1519035525.0*tpow[5]     - 49286948607.0*tpow[7] + 
          284499769554.0*tpow[9]   - 614135872350.0*tpow[11] + 
          566098157625.0*tpow[13]  - 188699385875.0*tpow[15])/6688604160.0;
}
#if 0
#ifdef HAVE_INLINE
inline
#endif
static double debye_u6(const double * tpow)
{
  return (2757049477875.0*tpow[6] - 127577298354750.0*tpow[8] + 
          1050760774457901.0*tpow[10] - 3369032068261860.0*tpow[12] + 
          5104696716244125.0*tpow[14] - 3685299006138750.0*tpow[16] + 
          1023694168371875.0*tpow[18])/4815794995200.0;
}
#endif

/* and Debye functions for imaginary argument */

#ifdef HAVE_INLINE
inline
#endif
static double debye_u1_im(const double * tpow)  /* u_1(i t)/i */
{
  return (3.0*tpow[1] + 5.0*tpow[3])/24.0;
}
#ifdef HAVE_INLINE
inline
#endif
static double debye_u2_im(const double * tpow)  /* u_2(i t)   */
{
  return (-81.0*tpow[2] - 462.0*tpow[4] - 385.0*tpow[6])/1152.0;
}
#ifdef HAVE_INLINE
inline
#endif
static double debye_u3_im(const double * tpow)  /* u_3(i t)/i */
{
  return (-30375.0*tpow[3] - 369603.0*tpow[5] - 765765.0*tpow[7] - 425425.0*tpow[9])/414720.0;
}
#ifdef HAVE_INLINE
inline
#endif
static double debye_u4_im(const double * tpow)  /* u_4(i t)   */
{
  return (4465125.0*tpow[4] + 94121676.0*tpow[6] + 349922430.0*tpow[8] + 
          446185740.0*tpow[10] + 185910725.0*tpow[12])/39813120.0;
}
#ifdef HAVE_INLINE
inline
#endif
static double debye_u5_im(const double * tpow)  /* u_5(i t)/i */
{
  return (1519035525.0*tpow[5]     + 49286948607.0*tpow[7] + 
          284499769554.0*tpow[9]   + 614135872350.0*tpow[11] + 
          566098157625.0*tpow[13]  + 188699385875.0*tpow[15])/6688604160.0;
}
#if 0
#ifdef HAVE_INLINE
inline
#endif
static double debye_u6_im(const double * tpow)
{
  return (-2757049477875.0*tpow[6]     - 127577298354750.0*tpow[8] -
           1050760774457901.0*tpow[10] - 3369032068261860.0*tpow[12] -
           5104696716244125.0*tpow[14] - 3685299006138750.0*tpow[16] - 
           1023694168371875.0*tpow[18])/4815794995200.0;
}
#endif


#if 0
/* [Abramowitz+Stegun, 9.3.17] */
static double debye_L(const double nu, const double * tpow)
{
  double nu2 = nu*nu;
  double nu4 = nu2*nu2;
  return 1. + debye_u2_im(tpow)/nu2 + debye_u4_im(tpow)/nu4;
}
#endif


#if 0
/* [Abramowitz+Stegun, 9.3.18] */
static double debye_M(const double nu, const double * tpow)
{
  double nu3 = nu*nu*nu;
  double nu5 = nu3*nu*nu;
  return debye_u1_im(tpow)/nu + debye_u3_im(tpow)/nu3 + debye_u5_im(tpow)/nu5;
}
#endif


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* Taylor expansion for J_nu(x) or I_nu(x)
 *   sign = -1  ==> Jnu
 *   sign = +1  ==> Inu
 *
 * error ~ o( (x/2)^(2N) / N! / (nu+1)^N )   N = kmax + 1
 *
 * empirical error analysis:
 *   for kmax=4, choose x*x < 10.*(n+1)*GSL_ROOT5_MACH_EPS
 *
 * Checks: nu >= 0; x >= 0
 *
 */
int gsl_sf_bessel_Inu_Jnu_taylor_impl(const double nu, const double x,
                                      const int sign,
                                      const int kmax,
				      const double threshold,
                                      double * result
                                      )
{
  if(nu < 0.0 || x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    if(nu == 0.0) {
      *result = 1.0;
    }
    else {
      *result = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if(nu == 0.0) {
    /* avoid unnecessary gamma computation */
    return Inu_Jnu_taylorsum(nu, x, sign, kmax, threshold, result);
  }
  else {
    /* nu > 0 and x > 0 */
    double pre;
    double sum;
    int stat_sum = Inu_Jnu_taylorsum(nu, x, sign, kmax, threshold, &sum);
    int stat_pre;
    int stat_m;

    /* prefactor = (x/2)^nu / Gamma(nu+1)
     */
    if(nu > INT_MAX-1) {
      double ln_pre;
      double lg;
      gsl_sf_lngamma_impl(nu+1.0, &lg);  /* ok by construction */
      ln_pre = nu*log(0.5*x) - lg;
      stat_pre = gsl_sf_exp_impl(ln_pre, &pre);
    }
    else {
      /* y^nu / Gamma(nu+1) = y^N /N! y^f / (N+1)_f */
      int    N = (int)floor(nu + 0.5);
      double f = nu - N;
      double poch_factor;
      double tc_factor;
      int stat_poch = gsl_sf_poch_impl(N+1.0, f, &poch_factor);
      int stat_tc   = gsl_sf_taylorcoeff_impl(N, 0.5*x, &tc_factor);
      pre = tc_factor * pow(0.5*x,f) / poch_factor;
      stat_pre = GSL_ERROR_SELECT_2(stat_tc, stat_poch);
    }

    stat_m = gsl_sf_multiply_impl(pre, sum, result);
    return GSL_ERROR_SELECT_3(stat_m, stat_pre, stat_sum);
  }
}


/* x >> nu*nu+1
 * error ~ O( ((nu*nu+1)/x)^3 )
 *
 * empirical error analysis:
 *   choose  GSL_ROOT3_MACH_EPS * x > (nu*nu + 1)
 *
 * This is not especially useful. When the argument gets
 * large enough for this to apply, the cos() and sin()
 * start loosing digits. However, this seems inevitable
 * for this particular method.
 *
 * checked OK [GJ] Sun May  3 22:35:00 EDT 1998 
 */
int gsl_sf_bessel_Jnu_asympx_impl(const double nu, const double x, double * result)
{
  double mu   = 4.0*nu*nu;
  double mum1 = mu-1.0;
  double mum9 = mu-9.0;
  double chi = x - (0.5*nu + 0.25)*M_PI;
  double P   = 1.0 - mum1*mum9/(128.0*x*x);
  double Q   = mum1/(8.0*x);
  *result = sqrt(2.0/(M_PI*x)) * (cos(chi)*P - sin(chi)*Q);
  return GSL_SUCCESS;
}

/* x >> nu*nu+1
 */
int gsl_sf_bessel_Ynu_asympx_impl(const double nu, const double x, double * result)
{
  double ampl;
  double theta;
  double alpha = x;
  double beta  = -0.5*nu*M_PI;
  int stat_a = gsl_sf_bessel_asymp_Mnu_impl(nu, x, &ampl);
  int stat_t = gsl_sf_bessel_asymp_thetanu_corr_impl(nu, x, &theta);
  int stat_red1 = gsl_sf_angle_restrict_pos_impl(&alpha);
  int stat_red2 = gsl_sf_angle_restrict_pos_impl(&beta);
  *result = ampl * sin(alpha+beta+theta);
  return GSL_ERROR_SELECT_4(stat_red1, stat_red2, stat_t, stat_a);
}

/* x >> nu*nu+1
 */
int gsl_sf_bessel_Inu_scaled_asympx_impl(const double nu, const double x, double * result)
{
  double mu   = 4.0*nu*nu;
  double mum1 = mu-1.0;
  double mum9 = mu-9.0;
  *result = 1.0/sqrt(2.0*M_PI*x) * (1.0 - mum1/(8.0*x) + mum1*mum9/(128.0*x*x));
  return GSL_SUCCESS;
}

/* x >> nu*nu+1
 */
int gsl_sf_bessel_Knu_scaled_asympx_impl(const double nu, const double x, double * result)
{
  double mu   = 4.0*nu*nu;
  double mum1 = mu-1.0;
  double mum9 = mu-9.0;
  *result = sqrt(M_PI/(2.0*x)) * (1.0 + mum1/(8.0*x) + mum1*mum9/(128.0*x*x));
  return GSL_SUCCESS;
}


/* nu -> Inf; x < nu   [Abramowitz+Stegun, 9.3.7]
 *
 * error:
 *   In the discussion below of the error in the Inu uniform expansion,
 *   we see how the error of the form u_N(t)/nu^N behaves. The difference
 *   here is that t > 1, so the polynomials are not well controlled. We
 *   must restrict ourselves to regions with x < nu since
 *   t = coth_a = 1/sqrt(1-(x/nu)^2). We can make some rough observations
 *   as follows.
 *
 *   We are going up to coth_a^{15}, ie debye_u5(), so we want coth_a to be
 *   near enough to 1 that this polynomial does not start to dominate.
 *   Roughly then, we want 1 + 15/2 (x/nu)^2 ~ 1, so then x << C nu,
 *   with C = sqrt(2/15) = 0.365.
 *
 *   Let's say that we only allow x < 0.5 nu. Then |u_6(t)| < 40, empirically.
 *   So the error is less than 40/nu^6. However, we are in a nonlinear
 *   regime already for x = 0.5nu, so we need to check things numerically.
 *   We will choose a condition valid for x < 0.5 nu, with an educated guess.
 *
 * empirical error analysis, assuming 14 digit requirement:
 *   choose  x < 0.500 nu   ==>  nu > 340
 *   choose  x < 0.365 nu   ==>  nu > 200
 *   choose  x < 0.250 nu   ==>  nu > 125
 *   choose  x < 0.100 nu   ==>  nu >  85
 *
 *   Our educated guess for the condition is, with x < 0.5 nu,
 *      0.4 / (nu^2 (1-(x/nu)^2)^6) < MACH_EPS^{1/3}
 *   or, equivalently,
 *      0.63 / (nu (1-(x/nu)^2)^3)  < MACH_EPS^{1/6}
 *
 *   The various powers here make sense since they correspond to the
 *   correct power of t appearing in the highest order term in u_6(t),
 *   which is t^{18}.
 *
 * This is superseded by the Olver uniform asymptotics.
 */
#if 0
int gsl_sf_bessel_Jnu_asymp_Debye_impl(const double nu, const double x, double * result)
{
  int i;
  double cosh_a = nu/x;        
  double sinh_a = sqrt(cosh_a*cosh_a-1.0);
  double coth_a = cosh_a/sinh_a;
  double tanh_a = 1.0/coth_a;
  double ea  = cosh_a + sinh_a;
  double lnt = nu*(tanh_a - log(ea));
  double pre = 1.0/sqrt(2.0*M_PI*nu*tanh_a) * exp(lnt);
  double tpow[16];
  double sum;
  tpow[0] = 1.;
  for(i=1; i<16; i++) tpow[i] = coth_a * tpow[i-1];
  sum = 1.0 + debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) + debye_u3(tpow)/(nu*nu*nu)
        + debye_u4(tpow)/(nu*nu*nu*nu) + debye_u5(tpow)/(nu*nu*nu*nu*nu);
  *result = pre * sum;
  return GSL_SUCCESS;
}
#endif /* 0 */


/* nu -> Inf; x < nu   [Abramowitz+Stegun, 9.3.7]
 *
 * error:
 *   same as discussion above for Jnu
 *
 * This is superseded by the Olver uniform asymptotics.
 */
#if 0
int gsl_sf_bessel_Ynu_asymp_Debye_impl(const double nu, const double x, double * result)
{
  int i;
  double cosh_a = nu/x;                         
  double sinh_a = sqrt(cosh_a*cosh_a-1.);       
  double coth_a = cosh_a/sinh_a;                
  double tanh_a = 1.0/coth_a;                    
  double ea  = cosh_a + sinh_a;                 
  double lnt = nu*(tanh_a - log(ea));
  double pre = sqrt(2.0/(M_PI*nu*tanh_a)) * exp(-lnt);
  double tpow[16];
  double sum;
  tpow[0] = 1.0;
  for(i=1; i<16; i++) tpow[i] = coth_a * tpow[i-1];
  sum = 1.0 - debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) - debye_u3(tpow)/(nu*nu*nu)
        + debye_u4(tpow)/(nu*nu*nu*nu) - debye_u5(tpow)/(nu*nu*nu*nu*nu);
  *result = -pre * sum;
  return GSL_SUCCESS;
}
#endif /* 0 */


/* nu -> Inf; x > nu   [Abramowitz+Stegun, 9.3.15]
 *
 * error:
 *
 *   Our educated guess for the condition is, with x > 2 nu,
 *      0.4 / (nu^2 (1-(nu/x)^2)^6) < MACH_EPS^{1/3}
 *
 * This is superseded by the Olver uniform asymptotics.
 */
#if 0
int gsl_sf_bessel_Jnu_asymp_Debye_osc_impl(const double nu, const double x, double * result)
{
  int i;
  double sec_b = x/nu;
  double tan_b = sqrt(sec_b*sec_b - 1.0);
  double t   = 1.0/tan_b;
  double Psi = -0.25*M_PI + nu*(tan_b - atan(tan_b));
  double pre = sqrt(2.0/(M_PI*nu*tan_b));
  double tpow[16];
  double L, M;
  tpow[0] = 1.0;
  for(i=0; i<16; i++) tpow[i] = t * tpow[i-1];
  L = debye_L(nu, tpow);
  M = debye_M(nu, tpow);
  *result = pre * (L*cos(Psi) + M*sin(Psi));
  return GSL_SUCCESS;
}
#endif /* 0 */


/* nu -> Inf; x > nu   [Abramowitz+Stegun, 9.3.16]
 *
 * error:
 *  same analysis as for Jnu above
 *
 * This is superseded by the Olver uniform asymptotics.
 */
#if 0
int gsl_sf_bessel_Ynu_asymp_Debye_osc_impl(const double nu, const double x, double * result)
{
  int i;
  double sec_b = x/nu;
  double tan_b = sqrt(sec_b*sec_b - 1.);
  double t   = 1./tan_b;
  double Psi = -0.25*M_PI + nu*(tan_b - atan(tan_b));
  double pre = sqrt(2./(M_PI*nu*tan_b));
  double tpow[16];
  double L, M;
  tpow[0] = 1.;
  for(i=0; i<16; i++) tpow[i] = t * tpow[i-1];
  L = debye_L(nu, tpow);
  M = debye_M(nu, tpow);
  *result = pre * (L*sin(Psi) - M*cos(Psi));
  return GSL_SUCCESS;
}
#endif /* 0 */


/* nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.7]
 *
 * error:
 *   The error has the form u_N(t)/nu^N  where  0 <= t <= 1.
 *   It is not hard to show that |u_N(t)| is small for such t.
 *   We have N=6 here, and |u_6(t)| < 0.025, so the error is clearly
 *   bounded by 0.025/nu^6. This gives the asymptotic bound on nu
 *   seen below as nu ~ 100. For general MACH_EPS it will be 
 *                     nu > 0.5 / MACH_EPS^(1/6)
 *   When t is small, the bound is even better because |u_N(t)| vanishes
 *   as t->0. In fact u_N(t) ~ C t^N as t->0, with C ~= 0.1.
 *   We write
 *                     err_N <= min(0.025, C(1/(1+(x/nu)^2))^3) / nu^6
 *   therefore
 *                     min(0.29/nu^2, 0.5/(nu^2+x^2)) < MACH_EPS^{1/3}
 *   and this is the general form.
 *
 * empirical error analysis, assuming 14 digit requirement:
 *   choose   x > 50.000 nu   ==>  nu >   3
 *   choose   x > 10.000 nu   ==>  nu >  15
 *   choose   x >  2.000 nu   ==>  nu >  50
 *   choose   x >  1.000 nu   ==>  nu >  75
 *   choose   x >  0.500 nu   ==>  nu >  80
 *   choose   x >  0.100 nu   ==>  nu >  83
 *
 * This makes sense. For x << nu, the error will be of the form u_N(1)/nu^N,
 * since the polynomial term will be evaluated near t=1, so the bound
 * on nu will become constant for small x. Furthermore, decreasing x with
 * nu fixed will decrease the error.
 *
 * checked OK [GJ] Sun May  3 21:31:53 EDT 1998 
 */
int gsl_sf_bessel_Inu_scaled_asymp_unif_impl(const double nu, const double x, double * result)
{
  int i;
  double z = x/nu;
  double root_term = sqrt(1.0 + z*z);
  double pre = 1.0/sqrt(2.0*M_PI*nu * root_term);
  double eta = root_term + log(z/(1.0+root_term));
  double ex  = ( z < 1.0/GSL_ROOT3_DBL_EPSILON ? exp(nu*(-z + eta)) : exp(-0.5*nu/z*(1.0 - 1.0/(12.0*z*z))) );
  double t = 1.0/root_term;
  double sum;
  double tpow[16];
  tpow[0] = 1.0;
  for(i=1; i<16; i++) tpow[i] = t * tpow[i-1];
  sum = 1.0 + debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) + debye_u3(tpow)/(nu*nu*nu)
        + debye_u4(tpow)/(nu*nu*nu*nu) + debye_u5(tpow)/(nu*nu*nu*nu*nu);
  *result = pre * ex * sum;
  return GSL_SUCCESS;
}

/* nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.8]
 *
 * error:
 *   identical to that above for Inu_scaled
 *
 * checked OK [GJ] Sun May  3 21:27:11 EDT 1998 
 */
int gsl_sf_bessel_Knu_scaled_asymp_unif_impl(const double nu, const double x, double * result)
{
  int i;
  double z = x/nu;
  double root_term = sqrt(1.0 + z*z);
  double pre = sqrt(M_PI/(2.0*nu*root_term));
  double eta = root_term + log(z/(1.0+root_term));
  double ex  = ( z < 1.0/GSL_ROOT3_DBL_EPSILON ? exp(nu*(z - eta)) : exp(0.5*nu/z*(1.0 + 1.0/(12.0*z*z))) );
  double t = 1.0/root_term;
  double sum;
  double tpow[16];
  tpow[0] = 1.0;
  for(i=1; i<16; i++) tpow[i] = t * tpow[i-1];
  sum = 1.0 - debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) - debye_u3(tpow)/(nu*nu*nu)
        + debye_u4(tpow)/(nu*nu*nu*nu) - debye_u5(tpow)/(nu*nu*nu*nu*nu);
  *result = pre * ex * sum;
  return GSL_SUCCESS;
}


/* Safely evaluate J_nu, Y_nu, J'_nu, Y'_nu at x = 0.
 * Assumes nu >= 0.
 */
int
gsl_sf_bessel_JnuYnu_zero(const double nu,
                          double * Jnu,  double * Ynu,
                          double * Jpnu, double * Ypnu
		          )
{
  int status = 0;

  if(Jnu != (double *)0) {
    if(nu == 0.0) *Jnu = 1.0;
    else *Jnu = 0.0;
  }
  if(Jpnu != (double *)0) {
    if(nu < 1.0) {
      *Jpnu = 0.0;
      status += 1;
    }
    else if(nu == 1.0) {
      *Jpnu = 0.5;
    }
    else {
      *Jpnu = 0.0;
    }
  }
  if(Ynu != (double *)0) {
    *Ynu = 0.0;
    status += 1;
  }
  if(Ypnu != (double *)0) {
    *Ypnu = 0.0;
    status += 1;
  }
  if(status)
    return GSL_EDOM;
  else
    return GSL_SUCCESS;
}


/* Evaluate the Steed method continued fraction CF2 for
 *
 * (J' + i Y')/(J + i Y) := P + i Q
 */
int
gsl_sf_bessel_JY_steed_CF2(const double nu, const double x,
                           double * P, double * Q)
{
  const int max_iter = 10000;
  const double SMALL = 1.0e-100;

  int i = 1;

  double x_inv = 1.0/x;
  double a = 0.25 - nu*nu;
  double p = -0.5*x_inv;
  double q = 1.0;
  double br = 2.0*x;
  double bi = 2.0;
  double fact = a*x_inv/(p*p + q*q);
  double cr = br + q*fact;
  double ci = bi + p*fact;
  double den = br*br + bi*bi;
  double dr = br/den;
  double di = -bi/den;
  double dlr = cr*dr - ci*di;
  double dli = cr*di + ci*dr;
  double temp = p*dlr - q*dli;
  q = p*dli + q*dlr;
  p = temp;
  for (i=2; i<=max_iter; i++) {
    a  += 2*(i-1);
    bi += 2.0;
    dr = a*dr + br;
    di = a*di + bi;
    if(fabs(dr)+fabs(di) < SMALL) dr = SMALL;
    fact = a/(cr*cr+ci*ci);
    cr = br + cr*fact;
    ci = bi - ci*fact;
    if(fabs(cr)+fabs(ci) < SMALL) cr = SMALL;
    den = dr*dr + di*di;
    dr /= den;
    di /= -den;
    dlr = cr*dr - ci*di;
    dli = cr*di + ci*dr;
    temp = p*dlr - q*dli;
    q = p*dli + q*dlr;
    p = temp;
    if(fabs(dlr-1.0)+fabs(dli) < GSL_DBL_EPSILON) break;
  }

  *P = p;
  *Q = q;

  if(i == max_iter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
 * to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
 *
 * This is unstable for small x; x > 2 is a good cutoff.
 * Also requires |nu| < 1/2.
 */
int
gsl_sf_bessel_K_scaled_steed_temme_CF2(const double nu, const double x,
                                       double * K_nu, double * K_nup1,
                                       double * Kp_nu)
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
    if(fabs(dels/s) < GSL_DBL_EPSILON) break;
  }
  
  hi *= -a1;
  
  *K_nu   = sqrt(M_PI/(2.0*x)) / s;
  *K_nup1 = *K_nu * (nu + x + 0.5 - hi)/x;
  *Kp_nu  = - *K_nup1 + nu/x * *K_nu;
  if(i == maxiter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/************************************************************************
 *                                                                      *
  Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from
  G.N.Watson, A Treatise on the Theory of Bessel Functions,
  2nd Edition (Cambridge University Press, 1944).
  Higher terms in expansion for x near l given by
  Airey in Phil. Mag. 31, 520 (1916).

  This approximation is accurate to near 0.1% at the boundaries
  between the asymptotic regions; well away from the boundaries
  the accuracy is better than 10^{-5}.
 *                                                                      *
 ************************************************************************/
#if 0
double besselJ_meissel(double nu, double x)
{
  double beta = pow(nu, 0.325);
  double result;

  /* Fitted matching points.   */
  double llimit = 1.1 * beta;
  double ulimit = 1.3 * beta;

  double nu2 = nu * nu;

  if (nu < 5. && x < 1.)
    {
      /* Small argument and order. Use a Taylor expansion. */
      int k;
      double xo2 = 0.5 * x;
      double gamfactor = pow(nu,nu) * exp(-nu) * sqrt(nu * 2. * M_PI)
	* (1. + 1./(12.*nu) + 1./(288.*nu*nu));
      double prefactor = pow(xo2, nu) / gamfactor;
      double C[5];

      C[0] = 1.;
      C[1] = -C[0] / (nu+1.);
      C[2] = -C[1] / (2.*(nu+2.));
      C[3] = -C[2] / (3.*(nu+3.));
      C[4] = -C[3] / (4.*(nu+4.));
      
      result = 0.;
      for(k=0; k<5; k++)
	result += C[k] * pow(xo2, 2.*k);

      result *= prefactor;
    }
  else if(x < nu - llimit)
    {
      /* Small x region: x << l.    */
      double z = x / nu;
      double z2 = z*z;
      double rtomz2 = sqrt(1.-z2);
      double omz2_2 = (1.-z2)*(1.-z2);

      /* Calculate Meissel exponent. */
      double term1 = 1./(24.*nu) * ((2.+3.*z2)/((1.-z2)*rtomz2) -2.);
      double term2 = - z2*(4. + z2)/(16.*nu2*(1.-z2)*omz2_2);
      double V_nu = term1 + term2;
      
      /* Calculate the harmless prefactor. */
      double sterlingsum = 1. + 1./(12.*nu) + 1./(288*nu2);
      double harmless = 1. / (sqrt(rtomz2*2.*M_PI*nu) * sterlingsum);

      /* Calculate the logarithm of the nu dependent prefactor. */
      double ln_nupre = rtomz2 + log(z) - log(1. + rtomz2);

      result = harmless * exp(nu*ln_nupre - V_nu);
    } 
  else if(x < nu + ulimit)
    {         
      /* Intermediate region 1: x near nu. */
      double eps = 1.-nu/x;
      double eps_x = eps * x;
      double eps_x_2 = eps_x * eps_x;
      double xo6 = x/6.;
      double B[6];
      static double gam[6] = {2.67894, 1.35412, 1., 0.89298, 0.902745, 1.};
      static double sf[6] = {0.866025, 0.866025, 0., -0.866025, -0.866025, 0.};
      
      /* Some terms are identically zero, because sf[] can be zero.
       * Some terms do not appear in the result.
       */
      B[0] = 1.;
      B[1] = eps_x;
      /* B[2] = 0.5 * eps_x_2 - 1./20.; */
      B[3] = eps_x * (eps_x_2/6. - 1./15.);
      B[4] = eps_x_2 * (eps_x_2 - 1.)/24. + 1./280.;
      /* B[5] = eps_x * (eps_x_2*(0.5*eps_x_2 - 1.)/60. + 43./8400.); */

      result  = B[0] * gam[0] * sf[0] / pow(xo6, 1./3.);
      result += B[1] * gam[1] * sf[1] / pow(xo6, 2./3.);
      result += B[3] * gam[3] * sf[3] / pow(xo6, 4./3.);
      result += B[4] * gam[4] * sf[4] / pow(xo6, 5./3.);

      result /= (3.*M_PI);
    }
  else 
    {
      /* Region of very large argument. Use expansion
       * for x>>l, and we need not be very exacting.
       */
      double secb = x/nu;
      double sec2b= secb*secb;
      
      double cotb = 1./sqrt(sec2b-1.);      /* cotb=cot(beta) */

      double beta = acos(nu/x);
      double trigarg = nu/cotb - nu*beta - 0.25 * M_PI;
      
      double cot3b = cotb * cotb * cotb;
      double cot6b = cot3b * cot3b;

      double sum1, sum2, expterm, prefactor, trigcos;

      sum1  = 2.0 + 3.0 * sec2b;
      trigarg -= sum1 * cot3b / (24.0 * nu);

      trigcos = cos(trigarg);

      sum2 = 4.0 + sec2b;
      expterm = sum2 * sec2b * cot6b / (16.0 * nu2);

      expterm = exp(-expterm);
      prefactor = sqrt(2. * cotb / (nu * M_PI));
      
      result = prefactor * expterm * trigcos;
    }

  return  result;
}
#endif
