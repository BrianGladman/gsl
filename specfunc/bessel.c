/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_airy.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_bessel.h"

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
                             double * result
                             )
{
  int k;
  double y = sign * 0.25 * x*x;
  
  double kfact  = 1.;
  double nuprod = 1.;
  double yprod  = 1.;
  double ans    = 1.;
  double delta;
 
  for(k=1; k<=kmax; k++) {
    nuprod *= (nu + k);
    kfact  *= k;
    yprod  *= y;
    delta = yprod / (kfact * nuprod);
    ans += delta;
  }
  
  *result = ans;
  
  if(fabs(delta) > 10. * GSL_MACH_EPS) {
    return GSL_ELOSS;
  }
  else {
    return GSL_SUCCESS;
  }
}


/* Debye functions [Abramowitz+Stegun, 9.3.9-10] */

static inline double debye_u1(const double * tpow)
{
  return (3*tpow[1] - 5*tpow[3])/24.;
}
static inline double debye_u2(const double * tpow)
{
  return (81*tpow[2] - 462*tpow[4] + 385*tpow[6])/1152.;
}
static inline double debye_u3(const double * tpow)
{
  return (30375.*tpow[3] - 369603.*tpow[5] + 765765.*tpow[7] - 425425.*tpow[9])/414720.;
}
static inline double debye_u4(const double * tpow)
{
  return (4465125.*tpow[4] - 94121676.*tpow[6] + 349922430.*tpow[8] - 
          446185740.*tpow[10] + 185910725.*tpow[12])/39813120.;
}
static inline double debye_u5(const double * tpow)
{
  return (1519035525.*tpow[5]     - 49286948607.*tpow[7] + 
          284499769554.*tpow[9]   - 614135872350.*tpow[11] + 
          566098157625.*tpow[13]  - 188699385875.*tpow[15])/6688604160.;
}
static inline double debye_u6(const double * tpow)
{
  return (2757049477875.*tpow[6] - 127577298354750.*tpow[8] + 
          1050760774457901.*tpow[10] - 3369032068261860.*tpow[12] + 
          5104696716244125.*tpow[14] - 3685299006138750.*tpow[16] + 
          1023694168371875.*tpow[18])/4815794995200.;
}

/* and Debye functions for imaginary argument */

static inline double debye_u1_im(const double * tpow)  /* u_1(i t)/i */
{
  return (3*tpow[1] + 5*tpow[3])/24.;
}
static inline double debye_u2_im(const double * tpow)  /* u_2(i t)   */
{
  return (-81*tpow[2] - 462*tpow[4] - 385*tpow[6])/1152.;
}
static inline double debye_u3_im(const double * tpow)  /* u_3(i t)/i */
{
  return (-30375.*tpow[3] - 369603.*tpow[5] - 765765.*tpow[7] - 425425.*tpow[9])/414720.;
}
static inline double debye_u4_im(const double * tpow)  /* u_4(i t)   */
{
  return (4465125.*tpow[4] + 94121676.*tpow[6] + 349922430.*tpow[8] + 
          446185740.*tpow[10] + 185910725.*tpow[12])/39813120.;
}
static inline double debye_u5_im(const double * tpow)  /* u_5(i t)/i */
{
  return (1519035525.*tpow[5]     + 49286948607.*tpow[7] + 
          284499769554.*tpow[9]   + 614135872350.*tpow[11] + 
          566098157625.*tpow[13]  + 188699385875.*tpow[15])/6688604160.;
}
static inline double debye_u6_im(const double * tpow)
{
  return (-2757049477875.*tpow[6] - 127577298354750.*tpow[8] -
           1050760774457901.*tpow[10] - 3369032068261860.*tpow[12] -
           5104696716244125.*tpow[14] - 3685299006138750.*tpow[16] - 
           1023694168371875.*tpow[18])/4815794995200.;
}

/* [Abramowitz+Stegun, 9.3.17] */
static double debye_L(const double nu, const double * tpow)
{
  double nu2 = nu*nu;
  double nu4 = nu2*nu2;
  return 1. + debye_u2_im(tpow)/nu2 + debye_u4_im(tpow)/nu4;
}

/* [Abramowitz+Stegun, 9.3.18] */
static double debye_M(const double nu, const double * tpow)
{
  double nu3 = nu*nu*nu;
  double nu5 = nu3*nu*nu;
  return debye_u1_im(tpow)/nu + debye_u3_im(tpow)/nu3 + debye_u5_im(tpow)/nu5;
}


/* functions for transition region [Abramowitz+Stegun, 9.3.25-26] */

static double trans_f1(const double * zpow)
{
  return -0.2*zpow[1];
}
static double trans_f2(const double * zpow)
{
  return -.09*zpow[5] + 3./35.*zpow[2];
}
static double trans_f3(const double * zpow)
{
  return 957./7000.*zpow[6] - 173./3150.*zpow[3] - 1./225.;
}
static double trans_f4(const double * zpow)
{
  return 27./20000.*zpow[10] - 23573./147000.*zpow[7] + 5903./138600.*zpow[4] + 947./346500.*zpow[1];
}
static double trans_g0(const double * zpow)
{
  return .3*zpow[2];
}
static double trans_g1(const double * zpow)
{
  return -17./70.*zpow[3] + 1./70.;
}
static double trans_g2(const double * zpow)
{
  return -9./1000.*zpow[7] + 611./3150.*zpow[4] - 37./3150.*zpow[1];
}
static double trans_g3(const double * zpow)
{
  return 549./28000.*zpow[8] - 110767./693000.*zpow[5] + 79./12375.*zpow[2];
}


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
 * checked OK [GJ] Sun May  3 22:34:46 EDT 1998 
 */
int gsl_sf_bessel_Inu_Jnu_taylor_impl(const double nu, const double x,
                                      const int sign,
                                      const int kmax,
                                      double * result
                                      )
{
  if(nu < 0.0 || x < 0.0) {
    return GSL_EDOM;
  }
  
  if(nu > 0. && x > 0.) {
    double g;
    int status = gsl_sf_lngamma_impl(nu+1., &g);  /* ok by construction */
    double ln_pre = nu*log(0.5*x) - g;
    if(ln_pre > GSL_LOG_DBL_MIN + 1.) {
      double pre = exp(ln_pre);
      double ts;
      status = Inu_Jnu_taylorsum(nu, x, sign, kmax, &ts);
      *result = pre * ts;
      return status;
    }
    else {
      *result = 0.;
      return GSL_EUNDRFLW;
    }
  }
  
  if(x == 0. && nu == 0.) {
    *result= 1.;
    return GSL_SUCCESS;
  }
  else if(x == 0.) {
    *result = 0.;
    return GSL_SUCCESS;
  }
  else if(nu == 0.) {
    int status = Inu_Jnu_taylorsum(nu, x, sign, kmax, result);
    return status;
  }
  else {
    return GSL_SUCCESS; /* NOT REACHED */
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
  double mu   = 4.*nu*nu;
  double mum1 = mu-1.;
  double mum9 = mu-9.;
  double chi = x - (0.5*nu + 0.25)*M_PI;
  double P   = 1. - mum1*mum9/(128.*x*x);
  double Q   = mum1/(8.*x);
  *result = sqrt(2./(M_PI*x)) * (cos(chi)*P - sin(chi)*Q);
  return GSL_SUCCESS;
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
 * checked OK [GJ] Mon May  4 00:06:52 EDT 1998 
 */
int gsl_sf_bessel_Ynu_asympx_impl(const double nu, const double x, double * result)
{
  double mu   = 4.*nu*nu;
  double mum1 = mu-1.;
  double mum9 = mu-9.;
  double chi = x - (0.5*nu + 0.25)*M_PI;
  double P   = 1. - mum1*mum9/(128.*x*x);
  double Q   = mum1/(8.*x);
  *result = sqrt(2./(M_PI*x)) * (sin(chi)*P + cos(chi)*Q);
  return GSL_SUCCESS;
}

/* x >> nu*nu+1
 * error ~ O( ((nu*nu+1)/x)^3 )
 *
 * empirical error analysis:
 *   choose  GSL_ROOT3_MACH_EPS * x > 0.25 * (nu*nu + 1)
 *
 * checked OK [GJ] Sun May  3 21:40:12 EDT 1998 
 */
int gsl_sf_bessel_Inu_scaled_asympx_impl(const double nu, const double x, double * result)
{
  double mu   = 4.*nu*nu;
  double mum1 = mu-1.;
  double mum9 = mu-9.;
  *result = 1./sqrt(2.*M_PI*x) * (1. - mum1/(8.*x) + mum1*mum9/(128.*x*x));
  return GSL_SUCCESS;
}

/* x >> nu*nu+1
 * error ~ O( ((nu*nu+1)/x)^3 )
 *
 * empirical error analysis:
 *   choose  GSL_ROOT3_MACH_EPS * x > 0.25 * (nu*nu + 1)
 *
 * checked OK [GJ] Sun May  3 21:37:16 EDT 1998 
 */
int gsl_sf_bessel_Knu_scaled_asympx_impl(const double nu, const double x, double * result)
{
  double mu   = 4.*nu*nu;
  double mum1 = mu-1.;
  double mum9 = mu-9.;
  *result = sqrt(M_PI/(2.*x)) * (1. + mum1/(8.*x) + mum1*mum9/(128.*x*x));
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
 */
int gsl_sf_bessel_Jnu_asymp_Debye_impl(const double nu, const double x, double * result)
{
  int i;
  double cosh_a = nu/x;                           
  double sinh_a = sqrt(cosh_a*cosh_a-1.);         
  double coth_a = cosh_a/sinh_a;                  
  double tanh_a = 1./coth_a;                      
  double ea  = cosh_a + sinh_a;                   
  double lnt = nu*(tanh_a - log(ea));
  double pre = 1./sqrt(2.*M_PI*nu*tanh_a) * exp(lnt);
  double tpow[16];
  double sum;
  tpow[0] = 1.;
  for(i=1; i<16; i++) tpow[i] = coth_a * tpow[i-1];
  sum = 1. + debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) + debye_u3(tpow)/(nu*nu*nu)
        + debye_u4(tpow)/(nu*nu*nu*nu) + debye_u5(tpow)/(nu*nu*nu*nu*nu);
  *result = pre * sum;
  return GSL_SUCCESS;
}

/* nu -> Inf; x < nu   [Abramowitz+Stegun, 9.3.7]
 *
 * error:
 *   same as discussion above for Jnu
 */
int gsl_sf_bessel_Ynu_asymp_Debye_impl(const double nu, const double x, double * result)
{
  int i;
  double cosh_a = nu/x;                         
  double sinh_a = sqrt(cosh_a*cosh_a-1.);       
  double coth_a = cosh_a/sinh_a;                
  double tanh_a = 1./coth_a;                    
  double ea  = cosh_a + sinh_a;                 
  double lnt = nu*(tanh_a - log(ea));
  double pre = sqrt(2./(M_PI*nu*tanh_a)) * exp(-lnt);
  double tpow[16];
  double sum;
  tpow[0] = 1.;
  for(i=1; i<16; i++) tpow[i] = coth_a * tpow[i-1];
  sum = 1. - debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) - debye_u3(tpow)/(nu*nu*nu)
        + debye_u4(tpow)/(nu*nu*nu*nu) - debye_u5(tpow)/(nu*nu*nu*nu*nu);
  *result = -pre * sum;
  return GSL_SUCCESS;
}


/* nu -> Inf; x > nu   [Abramowitz+Stegun, 9.3.15]
 *
 * error:
 *
 *   Our educated guess for the condition is, with x > 2 nu,
 *      0.4 / (nu^2 (1-(nu/x)^2)^6) < MACH_EPS^{1/3}
 */
int gsl_sf_bessel_Jnu_asymp_Debye_osc_impl(const double nu, const double x, double * result)
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
  printf("--  %g    %g  %g  %g  --  %g  %g\n",
         nu, L, M, Psi, tan_b, atan(tan_b)/((double)M_PI/2.)
	 );
  *result = pre * (L*cos(Psi) + M*sin(Psi));
  return GSL_SUCCESS;
}

/* nu -> Inf; x > nu   [Abramowitz+Stegun, 9.3.16]
 *
 * error:
 *  same analysis as for Jnu above
 */
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
  double root_term = sqrt(1. + z*z);
  double pre = 1./sqrt(2.*M_PI*nu * root_term);
  double eta = root_term + log(z/(1.+root_term));
  double ex  = ( z < 1./GSL_ROOT3_MACH_EPS ? exp(nu*(-z + eta)) : exp(-0.5*nu/z*(1. + 1./(12.*z*z))) );
  double t = 1./root_term;
  double sum;
  double tpow[16];
  tpow[0] = 1.;
  for(i=1; i<16; i++) tpow[i] = t * tpow[i-1];
  sum = 1. + debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) + debye_u3(tpow)/(nu*nu*nu)
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
  double root_term = sqrt(1. + z*z);
  double pre = sqrt(M_PI/(2.*nu*root_term));
  double eta = root_term + log(z/(1.+root_term));
  double ex  = ( z < 1./GSL_ROOT3_MACH_EPS ? exp(nu*(z - eta)) : exp(0.5*nu/z*(1. + 1./(12.*z*z))) );
  double t = 1./root_term;
  double sum;
  double tpow[16];
  tpow[0] = 1.;
  for(i=1; i<16; i++) tpow[i] = t * tpow[i-1];
  sum = 1. - debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) - debye_u3(tpow)/(nu*nu*nu)
        + debye_u4(tpow)/(nu*nu*nu*nu) - debye_u5(tpow)/(nu*nu*nu*nu*nu);
  *result = pre * ex * sum;
  return GSL_SUCCESS;
}

/* transition region
 *
 * error:
 *   The polynomials f and g are empirically bounded by 0.2 for 
 *   -1 <= z <= 1. So the maximum error is <= 0.2/nu^{10/3}. However,
 *   as one would expect, the polynomials die rapidly for |z| -> 0.
 *   For -1/2 <= z <= 1/2 we have |f|,|g| <= 0.01.
 *   For -1/4 <= z <= 1/4 we have |f|,|g| <= 0.001.
 *   
 */
int gsl_sf_bessel_Jnu_asymp_trans_impl(const double nu, const double x, double * result)
{
  int i;
  double nu_3 = pow(nu, 1./3.);
  double z  = (x - nu)/nu_3;
  double nu_3_2 = nu_3 * nu_3;
  double nu_3_4 = nu_3_2 * nu_3_2;
  double nu_3_8 = nu_3_4 * nu_3_4;
  double ai, aip;
  double sum1;
  double sum2;
  double zpow[11];
  zpow[0] = 1.;
  for(i=1; i<11; i++) zpow[i] = zpow[i-1] * z;
  sum1 = 1. + trans_f1(zpow)/(nu_3_2) 
            + trans_f2(zpow)/(nu_3_4)
            + trans_f3(zpow)/(nu*nu)
	    + trans_f4(zpow)/(nu_3_8);
  sum2 = trans_g0(zpow)
         + trans_g1(zpow)/(nu_3_2)
         + trans_g2(zpow)/(nu_3_4)
	 + trans_g3(zpow)/(nu*nu);
  gsl_sf_airy_Ai_impl(-CubeRoot2_*z, &ai);
  gsl_sf_airy_Ai_deriv_impl(-CubeRoot2_*z, &aip);
  *result = CubeRoot2_/nu_3 * ai * sum1
           + CubeRoot2_*CubeRoot2_/nu * aip * sum2;
  return GSL_SUCCESS; 
}

/* transition region
 */
int gsl_sf_bessel_Ynu_asymp_trans_impl(const double nu, const double x, double * result)
{
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


double besselJprime_meissel(double nu, double x, double J_nu)
{
  return  besselJ_meissel(nu-1., x) - nu * J_nu / x;
}


double sphbesselj_meissel(double el, double x)
{
  return sqrt(M_PI/(2. * x)) * besselJ_meissel(el+0.5, x);
}


double sphbesseljprime_meissel(double el, double x, double j_l)
{
  return sphbesselj_meissel(el-1., x) - j_l * (el+1.)/x;
}


void asymp_sphbesselj_meissel(double l, double x,
			      double *jl, double *jlp, int dflag)
{
  double arg = x - (l+0.5)*0.5*M_PI;
  double x2 = x*x;
  double c = cos(arg);
  double s = sin(arg);

  *jl = (c + s/(8.*x))/x;
  if(dflag){
    *jlp = -(s*(8.*x2 + 2.) + 7.*x*c)/ (8.*x2*x);
  }
}


/* Safely evaluate J_nu, Y_nu, J'_nu, Y'_nu at x = 0.
 * Assumes nu >= 0.
 */
int
gsl_sf_bessel_JnuYnu_zero(const double nu, const double x,
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


/* Perform backward recurrence for J_nu(x) and J'_nu(x)
 *
 *        J_{nu-1} = nu/x J_nu + J'_nu
 *       J'_{nu-1} = (nu-1)/x J_{nu-1} - J_nu
 *
 * J_array[kmax] = J_{nu + kmax}(x) = J_start
 *   ...
 * J_array[0]    = J_nu(x)          = J_end
 *
 * and similarly for Jp_array = J'
 */
int
gsl_sf_bessel_J_recur(const double nu_min, const double x, const int kmax,
                      const double J_start, const double Jp_start,
	              double * J_end, double * Jp_end,
	              double * J_array, double * Jp_array)
{
  double x_inv  = 1.0/x;
  double nu_max = nu_min + kmax;
  double nu     = nu_max;
  double J_nu  = J_start;
  double Jp_nu = Jp_start;
  int k;

  if(J_array  != (double *)0)  J_array[kmax] = J_start;
  if(Jp_array != (double *)0) Jp_array[kmax] = Jp_start;
  for(k=kmax-1; k>=0; k--) {
    double nu_x_inv = nu*x_inv;
    double J_nu_save = J_nu;
    J_nu  = nu_x_inv*J_nu + Jp_nu;
    Jp_nu = (nu_x_inv - x_inv)*J_nu - J_nu_save;
    if(J_array  != (double *)0)  J_array[k] = J_nu;
    if(Jp_array != (double *)0) Jp_array[k] = Jp_nu;
    nu -= 1.0;
  }
  *J_end  = J_nu;
  *Jp_end = Jp_nu;
  return GSL_SUCCESS;
}


/* Perform forward recurrence for Y_nu(x) and Y'_nu(x)
 *
 *        Y_{nu+1} =  nu/x Y_nu - Y'_nu
 *       Y'_{nu+1} = -nu/x Y_{nu+1} + Y_nu
 *
 * Y_array[0]    = Y_nu(x)          = Y_start
 *   ...
 * Y_array[kmax] = Y_{nu + kmax}(x) = Y_end
 *
 * and similarly for Yp_array = Y'
 */
int
gsl_sf_bessel_Y_recur(const double nu_min, const double x, const int kmax,
                      const double Y_start, const double Yp_start,
		      double * Y_end, double * Yp_end,
                      double * Y_array, double * Yp_array)
{
  double x_inv = 1.0/x;
  double nu = nu_min;
  double Y_nu  = Y_start;
  double Yp_nu = Yp_start;
  int k;

  if(Y_array  != (double *)0)  Y_array[0] = Y_start;
  if(Yp_array != (double *)0) Yp_array[0] = Yp_start;

  for(k=1; k<=kmax; k++) {
    double nuox = nu*x_inv;
    double Y_nu_save = Y_nu;
    Y_nu  = -Yp_nu + nuox * Y_nu;
    Yp_nu = Y_nu_save - nuox * Y_nu;
    if(Y_array  != (double *)0)  Y_array[k] = Y_nu;
    if(Yp_array != (double *)0) Yp_array[k] = Yp_nu;
    nu += 1.0;
  }
  *Y_end  = Y_nu;
  *Yp_end = Yp_nu;
  return GSL_SUCCESS;
}


/* Perform backward recurrence for J_nu(x) and J'_nu(x)
 *
 *        I_{nu-1} = nu/x I_nu + I'_nu
 *       I'_{nu-1} = (nu-1)/x I_{nu-1} + I_nu
 *
 * I_array[kmax] = I_{nu + kmax}(x) = I_start
 *   ...
 * I_array[0]    = I_nu(x)          = I_end
 *
 * and similarly for Ip_array = I'
 */
int
gsl_sf_bessel_I_recur(const double nu_min, const double x, const int kmax,
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
    double nuox = nu*x_inv;
    double I_nu_save = I_nu;
    I_nu  = nuox*I_nu + Ip_nu;
    Ip_nu = (nuox - x_inv)*I_nu + I_nu_save;
    if(I_array  != (double *)0)  I_array[k] = I_nu;
    if(Ip_array != (double *)0) Ip_array[k] = Ip_nu;
    nu -= 1.0;
  }
  *I_end  = I_nu;
  *Ip_end = Ip_nu;
  return GSL_SUCCESS;
}


/* Perform forward recurrence for K_nu(x) and K'_nu(x)
 *
 *        K_{nu+1} =  nu/x K_nu - K'_nu
 *       K'_{nu+1} = -nu/x K_{nu+1} - K_nu
 *
 * K_array[0]    = K_nu(x)          = K_start
 *   ...
 * K_array[kmax] = K_{nu + kmax}(x) = K_end
 *
 * and similarly for Kp_array = K'
 */
int
gsl_sf_bessel_K_recur(const double nu_min, const double x, const int kmax,
                      const double K_start, const double Kp_start,
		      double * K_end, double * Kp_end,
                      double * K_array, double * Kp_array)
{
  double x_inv = 1.0/x;
  double nu = nu_min;
  double K_nu  = K_start;
  double Kp_nu = Kp_start;
  int k;

  if(K_array  != (double *)0)  K_array[0] = K_start;
  if(Kp_array != (double *)0) Kp_array[0] = Kp_start;

  for(k=1; k<=kmax; k++) {
    double nuox = nu*x_inv;
    double K_nu_save = K_nu;
    K_nu  = -Kp_nu + nuox * K_nu;
    Kp_nu = -K_nu_save - nuox * K_nu;
    if(K_array  != (double *)0)  K_array[k] = K_nu;
    if(Kp_array != (double *)0) Kp_array[k] = Kp_nu;
    nu += 1.0;
  }
  *K_end  = K_nu;
  *Kp_end = Kp_nu;
  return GSL_SUCCESS;
}
