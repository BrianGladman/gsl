/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Evaluation of Coulomb wave functions F_L(eta, x), G_L(eta, x),
 * and their derivatives. Uses Steed's method, following the
 * analysis of [Barnett, Comp. Phys. Comm., 21, 297 (1981)].
 * Refinements for asymptotic regimes are also used.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_psi.h"
#include "gsl_sf_airy.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_coulomb.h"

#define locMax(a, b) ((a) > (b) ? (a) : (b))
#define locMin(a, b) ((a) < (b) ? (a) : (b))


/* the L=0 normalization constant */
static
double
C0sq(double eta)
{
  double twopieta = 2.*M_PI*eta;

  if(fabs(eta) < GSL_MACH_EPS) {
    return 1.0;
  }
  else if(twopieta > GSL_LOG_DBL_MAX) {
    return 0.0;
  }
  else {
    double scale = 0.0;
    gsl_sf_expm1_impl(twopieta, &scale);
    return twopieta/scale;
  }
}


/* the full definition of C_L(eta) for any valid L and eta
 * [Abramowitz and Stegun 14.1.7]
 * This depends on the complex gamma function. For large
 * arguments the phase of the complex gamma function is not
 * very accurately determined. However the modulus is, and that
 * is all that we need to calculate C_L.
 */
static
double
CLeta(double L, double eta)
{
  double ln1; /* log of numerator Gamma function */
  double ln2; /* log of denominator Gamma function */

  if(fabs(eta) < GSL_MACH_EPS) {
    gsl_sf_lngamma_impl(L+1.0, &ln1);
  }
  else {
    double p1;                 /* phase of numerator Gamma -- not used */
    gsl_sf_lngamma_complex_impl(L+1.0, eta, &ln1, &p1); /* should be ok */
  }
  gsl_sf_lngamma_impl(2*L+2.0, &ln2);
  
  return exp(L*M_LN2 - 0.5*eta*M_PI + ln1 - ln2);
}


int
gsl_sf_coulomb_CL_impl(double lam, double eta, double * result)
{
  if(lam <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(fabs(lam) < 10.0*GSL_MACH_EPS) {
    /* saves a calculation of complex_lngamma(), otherwise not necessary */
    *result = sqrt(C0sq(eta));
    return GSL_SUCCESS;
  }
  else {
    *result = CLeta(lam, eta);
    return GSL_SUCCESS;
  }
}


/* cl[0] .. cl[kmax] = C_{lam_min}(eta) .. C_{lam_min+kmax}(eta)
 */
int
gsl_sf_coulomb_CL_list(double lam_min, int kmax, double eta, double * cl)
{
  int k;
  gsl_sf_coulomb_CL_impl(lam_min, eta, cl);

  for(k=1; k<=kmax; k++) {
    double L = lam_min + k;
    cl[k] = cl[k-1] * sqrt(L*L + eta*eta)/(L*(2.0*L+1.0));
  }
  
  return GSL_SUCCESS;
}


/* Evaluate the series for Phi_L(eta,x) and Phi_L*(eta,x)
 * [Abramowitz+Stegun 14.1.5]
 * [Abramowitz+Stegun 14.1.13]
 *
 * The sequence of coefficients A_k^L is
 * manifestly well-controlled for L >= -1/2
 * and eta < 10.
 *
 * This makes sense since this is the region
 * away from threshold, and you expect
 * the evaluation to become easier as you
 * get farther from threshold.
 *
 * Empirically, this is quite well-behaved for
 *   L >= -1/2
 *   eta < 10
 *   x   < 10
 */
static
int
coulomb_Phi_series(const double lam, const double eta, const double x,
                   double * result, double * result_star)
{
  int kmin =   5;
  int kmax = 200;
  int k;
  double Akm2 = 1.0;
  double Akm1 = eta/(lam+ 1.0);
  double Ak;

  double xpow = x;
  double sum  = Akm2 + Akm1*x;
  double sump = (lam+1.0)*Akm2 + (lam+2.0)*Akm1*x;
  double prev_abs_del   = fabs(Akm1*x);
  double prev_abs_del_p = (lam+2.0) * prev_abs_del;

  for(k=2; k<kmax; k++) {
    double del;
    double del_p;
    double abs_del;
    double abs_del_p;

    Ak = (2.0*eta*Akm1 - Akm2)/(k*(2.0*lam + 1.0 + k));

    xpow *= x;
    del   = Ak*xpow;
    del_p = (k+lam+1.0)*del;
    sum  += del;
    sump += del_p;

    abs_del   = fabs(del);
    abs_del_p = fabs(del_p);

    if(          abs_del/(fabs(sum)+abs_del)          < GSL_MACH_EPS
       &&   prev_abs_del/(fabs(sum)+prev_abs_del)     < GSL_MACH_EPS
       &&      abs_del_p/(fabs(sump)+abs_del_p)       < GSL_MACH_EPS
       && prev_abs_del_p/(fabs(sump)+prev_abs_del_p)  < GSL_MACH_EPS
       && k > kmin
       ) break;

    /* We need to keep track of the previous delta because when
     * eta is near zero the odd terms of the sum are very small
     * and this could lead to premature termination.
     */
    prev_abs_del   = abs_del;
    prev_abs_del_p = abs_del_p;

    Akm2 = Akm1;
    Akm1 = Ak;
  }

  *result      = sum;
  *result_star = sump;

  if(k==kmax) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}


/* Evaluate the Frobenius series for G_lam(eta,x).
 *
 * G_lam = x^(-lam) Sum[ v_m, {m,0,Infinity}]
 *
 * v_m := b_m x^m
 * b_m m (m-1-2lam) - 2 eta b_{m-1} + b_{m-2} = 0
 *
 * Not appropriate for lam <= -1/2, lam = 0, or lam >= 1/2.
 */
static
int
coulomb_G_series(const double lam, const double eta, const double x, double *result)
{
  const int max_iter = 500;
  double Clam = CLeta(lam, eta);
  double vmm2 = 1.0/(2.0*lam+1.0) / Clam;  /* v_0 */
  double vmm1 = -x*eta/lam * vmm2;         /* v_1 */
  double vm;
  double sum  = vmm2 + vmm1;
  int m = 2;
  
  while(m < max_iter) {
    vm = x*(2.0*eta*vmm1 - x*vmm2)/(m*(m-1-2.0*lam));
    sum += vm;
    if(fabs(vm/sum) < 100.0*GSL_MACH_EPS) break;
    vmm2 = vmm1;
    vmm1 = vm;
    m++;
  }

  *result = pow(x, -lam) * sum / Clam;
  if(m == max_iter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* Evaluate the Frobenius series for G_0(eta,x)
 */
static
int
coulomb_G0_series(const double eta, const double x, double * result)
{
  const int max_iter = 500;
  double C0 = CLeta(0.0, eta);
  double r1pie;
  int psi_stat = gsl_sf_psi_1piy_impl(eta, &r1pie);
  double nu_mm2 = 0.0;                                 /* nu_0 */
  double nu_mm1 = 2.0*eta*x*(2.0*M_EULER-1.0+r1pie);   /* nu_1 */
  double nu_m;
  double nusum  = nu_mm2 + nu_mm1;
  double u_mm2 = 0.0;  /* u_0 */
  double u_mm1 = x;    /* u_1 */
  double u_m;
  double usum = u_mm2 + u_mm1;
  int m = 2;
  
  while(m < max_iter) {
    u_m  = (2.0*eta*x*u_mm1 - x*x*u_mm2)/(m*(m-1.0));
    nu_m = (2.0*eta*x*nu_mm1 - x*x*nu_mm2 - 2.0*eta*(2.0*m-1.0)*u_m)/(m*(m-1.0));
    nusum += nu_m;
    usum  += u_m;
    if(fabs(nu_m/nusum) + fabs(u_m/usum) < 100.0*GSL_MACH_EPS) break;
    nu_mm2 = nu_mm1;
    nu_mm1 = nu_m;
    u_mm2 = u_mm1;
    u_mm1 = u_m;
    m++;
  }
  
  *result = (nusum + 2.0*eta*usum * log(2.0*x)) / C0;
  return psi_stat;
}


/* Evolve the backwards recurrence for F,F'.
 *
 *    F_{lam-1}  = (S_lam F_lam + F_lam') / R_lam
 *    F_{lam-1}' = (S_lam F_{lam-1} - R_lam F_lam)
 * where
 *    R_lam = sqrt(1 + (eta/lam)^2)
 *    S_lam = lam/x + eta/lam
 *
 */
static
int
coulomb_F_recur(double lam_min, int kmax,
                double eta, double x,
                double F_lam_max, double Fp_lam_max,
		double * F_lam_min, double * Fp_lam_min
                )
{
  double x_inv = 1.0/x;
  double fcl = F_lam_max;
  double fpl = Fp_lam_max;
  double lam_max = lam_min + kmax;
  double lam = lam_max;
  int k;

  for(k=kmax-1; k>=0; k--) {
    double el = eta/lam;
    double rl = sqrt(1.0 + el*el);
    double sl = el  + lam*x_inv;
    double fc_lm1;
    fc_lm1 = (fcl*sl + fpl)/rl;
    fpl    =  fc_lm1*sl - fcl*rl;
    fcl    =  fc_lm1;
    lam -= 1.0;
  }

  *F_lam_min  = fcl;
  *Fp_lam_min = fpl;  
  return GSL_SUCCESS;
}


/* Evolve the forward recurrence for G,G'.
 *
 *   G_{lam+1}  = (S_lam G_lam - G_lam')/R_lam
 *   G_{lam+1}' = R_{lam+1} G_lam - S_lam G_{lam+1}
 *
 * where S_lam and R_lam are as above in the F recursion.
 */
static
int
coulomb_G_recur(const double lam_min, const int kmax,
                const double eta, const double x,
                const double G_lam_min, const double Gp_lam_min,
		double * G_lam_max, double * Gp_lam_max
                )
{
  double x_inv = 1.0/x;
  double gcl = G_lam_min;
  double gpl = Gp_lam_min;
  double lam = lam_min + 1.0;
  int k;

  for(k=1; k<=kmax; k++) {
    double el = eta/lam;
    double rl = sqrt(1.0 + el*el);
    double sl = el + lam*x_inv;
    double gcl1 = (sl*gcl - gpl)/rl;
    gpl   = rl*gcl - sl*gcl1;
    gcl   = gcl1;
    lam += 1.0;
  }
  
  *G_lam_max  = gcl;
  *Gp_lam_max = gpl;
  return GSL_SUCCESS;
}


/* Evaluate the first continued fraction, giving
 * the ratio F'/F at the upper lambda value.
 * We also determine the sign of F at that point,
 * since it is the sign of the last denominator
 * in the continued fraction.
 */
static
int
my_coulomb_CF1(double lambda,
               double eta, double x,
               double * fcl_sign,
	       double * result
               )
{
  const double CF1_small = 1.e-30;
  const double CF1_abort = 1.0e+05;
  const double CF1_acc   = 10.0*GSL_MACH_EPS;
  const double x_inv     = 1.0/x;
  const double px        = lam_max + 1.0 + CF1_abort;

  double pk = lambda + 1.0;
  double F  = eta/pk + pk*x_inv;
  double D, C;
  double df;

  *fcl_sign = 1.0;

  if(fabs(F) < CF1_small) F = CF1_small;
  D = 0.0;
  C = F;

  do {
    double pk1 = pk + 1.0;
    double ek  = eta / pk;
    double rk2 = 1.0 + ek*ek;
    double tk  = (pk + pk1)*(x_inv + ek/pk1);
    D   =  tk - rk2 * D;
    C   =  tk - rk2 / C;
    if(fabs(C) < CF1_small) C = CF1_small;
    if(fabs(D) < CF1_small) D = CF1_small;
    D = 1.0/D;
    df = D * C;
    F  = F * df;
    if(D < 0.0) {
      /* sign of result depends on sign of denominator */
      *fcl_sign = - *fcl_sign;
    }
    pk = pk1;
    if( pk > px ) {
      *result = F;
      return GSL_ERUNAWAY;
    }
  }
  while(fabs(df-1.0) > CF1_acc);
  
  *result = F;
  return GSL_SUCCESS;
}


static
int
coulomb_CF1(const double lambda,
            double eta, double x,
            double * fcl_sign,
	    double * result
            )
{
  const double CF1_abort = 1.e5;
  const double CF1_acc   = 10.0*GSL_MACH_EPS;
  const double x_inv     = 1.0/x;
  const double px        = lambda + 1.0 + CF1_abort;
  
  double pk = lambda + 1.0;
  
  double D;
  double df;

  double F;
  double p;
  double pk1;
  double ek;
  
  double fcl = 1.0;

  double tk;

  while(1) {
    ek = eta/pk;
    F = (ek + pk*x_inv)*fcl + (fcl - 1.0)*x_inv;
    pk1 = pk + 1.0;
    if(fabs(eta*x + pk*pk1) > CF1_acc) break;
    fcl = (1.0 + ek*ek)/(1.0 + eta*eta/(pk1*pk1));
    pk = 2.0 + pk;
  }

  D  = 1.0/((pk + pk1)*(x_inv + ek/pk1));
  df = -fcl*(1.0 + ek*ek)*D;
  
  if(fcl != 1.0) fcl = -1.0;
  if(D    < 0.0) fcl = -fcl;
  
  F = F + df;

  p = 1.0;
  do {
    pk = pk1;
    pk1 = pk + 1.0;
    ek  = eta / pk;
    tk  = (pk + pk1)*(x_inv + ek/pk1);
    D   =  tk - D*(1.0+ek*ek);
    if(fabs(D) < sqrt(CF1_acc)) {
      p += 1.0;
      if(p > 2.0) {
        printf("HELP............\n");
      }
    }
    D = 1.0/D;
    if(D < 0.0) {
      /* sign of result depends on sign of denominator */
      fcl = -fcl;
    }
    df = df*(D*tk - 1.0);
    F  = F + df;
    if( pk > px ) {
      printf("OH NO.....\n");
      return GSL_ERUNAWAY;
    }
  }
  while(fabs(df) > fabs(F)*CF1_acc);
  
  *fcl_sign = fcl;
  *result = F;
  return GSL_SUCCESS;
}


/* Evaluate the second continued fraction to 
 * obtain the ratio
 *    (G' + i F')/(G + i F) := P + i Q
 * at the specified lambda value.
 */
static
int
coulomb_CF2(const double lambda, const double eta, const double x,
            double * result_P, double * result_Q
            )
{
  int status = GSL_SUCCESS;

  const double CF2_acc   = 10.0*GSL_MACH_EPS;
  const double CF2_abort = 1.0e5;

  const double wi    = 2.0*eta;
  const double x_inv = 1.0/x;
  const double e2mm1 = eta*eta + lambda*(lambda + 1.0);
  
  double ar = -e2mm1;
  double ai =  eta;

  double br =  2.0*(x - eta);
  double bi =  2.0;

  double dr =  br/(br*br + bi*bi);
  double di = -bi/(br*br + bi*bi);

  double dp = -x_inv*(ar*di + ai*dr);
  double dq =  x_inv*(ar*dr - ai*di);

  double A, B, C, D;

  double pk =  0.0;
  double P  =  0.0;
  double Q  =  1.0 - eta*x_inv;
  
  do {
    P += dp;
    Q += dq;
    pk += 2.0;
    ar += pk;
    ai += wi;
    bi += 2.0;
    D  = ar*dr - ai*di + br;
    di = ai*dr + ar*di + bi;
    C  = 1.0/(D*D + di*di);
    dr =  C*D;
    di = -C*di;
    A  = br*dr - bi*di - 1.;
    B  = bi*dr + br*di;
    C  = dp*A  - dq*B;
    dq = dp*B  + dq*A;
    dp = C;
    if(pk > 2.0*CF2_abort) {
      status = GSL_ERUNAWAY;
      break;
    }
  }
  while(fabs(dp)+fabs(dq) > (fabs(P)+fabs(Q))*CF2_acc);

  if(Q < CF2_abort*GSL_MACH_EPS*fabs(P)) {
    status = GSL_ELOSS;
  }

  *result_P = P;
  *result_Q = Q;
  return status;
}


/* WKB evaluation of F, G. Assumes  0 < x < turning point.
 * Overflows are trapped, GSL_EOVRFLW is signalled,
 * and an exponent is returned such that:
 *
 *   result_F = fjwkb * exp(-exponent)
 *   result_G = gjwkb * exp( exponent)
 *
 * See [Biedenharn et al. Phys. Rev. 97, 542-554 (1955)]
 */
static
int
coulomb_jwkb(double lam, double eta, double x,
             double * fjwkb, double * gjwkb,
	     double * exponent)
{
  double gh2  = x*(2.0*eta - x);
  double llp1 = lam*(lam + 1.0);
  double xll1 = locMax(llp1, 0.0);
  double hll  = xll1 + 6.0/35.0;
  double hl   = sqrt(hll);
  double sl   = eta/hl + hl/x;
  double rl2  = 1.0 + eta*eta/hll;
  double gh   = sqrt(gh2 + hll)/x;
  double phi  = x*gh - 0.5*(hl*log((gh+sl)*(gh+sl)/rl2) - log(gh));

  if(eta != 0.0) phi -= eta*atan2(x*gh, x - eta);

  if(-phi < GSL_LOG_DBL_MAX) {
    *gjwkb = exp(-phi);
    *fjwkb = 0.5/(gh * *gjwkb);
    *exponent = 0.0;
    return GSL_SUCCESS;
  }
  else {
    *exponent = -phi;
    *gjwkb = 1.0;
    *fjwkb = 0.5/gh;
    return GSL_EOVRFLW;
  }
}


static
int
new_jwkb(const double lam, const double eta, const double x,
         double * fjwkb, double * gjwkb,
	 double * exponent)
{
  const double llp1      = lam*(lam+1.0) + 6.0/35.0;
  const double llp1_eff  = locMax(llp1, 0.0);
  const double rho_ghalf = sqrt(x*(2.0*eta - x) + llp1_eff);
  const double sinh_arg  = sqrt(llp1_eff/(eta*eta+llp1_eff)) * rho_ghalf / x;
  const double sinh_inv  = log(sinh_arg + sqrt(1.0 + sinh_arg*sinh_arg));

  const double phi = fabs(rho_ghalf - eta*atan2(rho_ghalf,x-eta) - sqrt(llp1_eff) * sinh_inv);

  const double zeta_half = pow(3.0*phi/2.0, 1.0/3.0);
  const double prefactor = sqrt(M_PI*phi*x/(6.0 * rho_ghalf));
  
  double F = prefactor *     3.0/zeta_half;
  double G = prefactor * M_SQRT3/zeta_half;
  double F_exp;
  double G_exp;
  
  const double airy_scale_exp = phi;
  double ai, bi;
  gsl_sf_airy_Ai_scaled_impl(zeta_half*zeta_half, &ai);
  gsl_sf_airy_Bi_scaled_impl(zeta_half*zeta_half, &bi);
  F *= ai;
  G *= bi;
  F_exp = log(F) - airy_scale_exp;
  G_exp = log(G) + airy_scale_exp;

  if(G_exp >= GSL_LOG_DBL_MAX) {
    *fjwkb = F;
    *gjwkb = G;
    *exponent = airy_scale_exp;
    return GSL_EOVRFLW;
  }
  else {
    *fjwkb = exp(F_exp);
    *gjwkb = exp(G_exp);
    *exponent = 0.0;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_coulomb_wave_FG_impl(const double eta, const double x,
                            const double lam_F,
			    const int  k_lam_G,      /* lam_G = lam_F - k_lam_G */
                            double * F, double * Fp,
			    double * G, double * Gp,
			    double * exp_F, double * exp_G)
{
  const double lam_G = lam_F - k_lam_G;

  if(x < 0.0 || lam_F <= -0.5 || lam_G <= -0.5) {
    *F  = 0.0;
    *Fp = 0.0;
    *G  = 0.0;
    *Gp = 0.0;
    *exp_F = 0.0;
    *exp_G = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    double C0 = CLeta(0.0, eta);
    *exp_F = 0.0;
    *exp_G = 0.0;
    *F  = 0.0;
    *Fp = 0.0;
    *G  = 0.0; /* FIXME: should be Inf */
    *Gp = 0.0; /* FIXME: should be Inf */
    if(lam_F == 0.0){
      *F  = 0.0;
      *Fp = C0;
    }
    if(lam_G == 0.0) {
      *G  = 1.0/C0;
    }
    return GSL_EDOM;
    /* After all, since we are asking for G, this is a domain error... */
  }
  else if(x < 1.0 && fabs(eta*x) < 5.0) {
    /* Reduce to a small lambda value and use the series
     * representations for F and G.
     */
    const double SMALL = 1.0e-100;
    const int N    = (int)(lam_F + 0.5);
    const int span = locMax(k_lam_G, N);
    const double lam_min = lam_F - N;    /* -1/2 <= lam_min <= 1/2 */
    const double C_lam_min = CLeta(lam_min, eta);
    const double pow_x = pow(x, lam_min);
    double F_lam_F, Fp_lam_F;
    double G_lam_G, Gp_lam_G;
    double Fp_over_F_lam_F;
    double F_sign_lam_F;
    double F_lam_min_unnorm, Fp_lam_min_unnorm;
    double Phi, Phi_star;
    double F_lam_min, Fp_lam_min;
    double G_lam_min, Gp_lam_min;
    double F_scale;
    
    /* Determine F'/F at lam_F. */
    coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F);

    /* Recurse down with unnormalized F,F' values. */
    F_lam_F  = SMALL;
    Fp_lam_F = Fp_over_F_lam_F * F_lam_F;
    coulomb_F_recur(lam_min, span, eta, x,
                    F_lam_F, Fp_lam_F,
		    &F_lam_min_unnorm, &Fp_lam_min_unnorm
		    );

    /* Determine F,F' at lam_min. */
    coulomb_Phi_series(lam_min, eta, x, &Phi, &Phi_star);
    F_lam_min  = C_lam_min * pow_x * x * Phi;
    Fp_lam_min = C_lam_min * pow_x * Phi_star;
    F_scale = F_lam_min / F_lam_min_unnorm;

    /* Determine G and G' at lam_min. */
    if(lam_min == 0.0) {
      coulomb_G0_series(eta, x, &G_lam_min);
    }
    else {
      coulomb_G_series(lam_min, eta, x, &G_lam_min);
    }
    Gp_lam_min = (Fp_lam_min * G_lam_min - 1.0)/F_lam_min;
 
    /* Apply scale to the original F,F' values. */
    F_lam_F  *= F_scale;
    Fp_lam_F *= F_scale;

    /* Recurse up to get the required G,G' values. */
    coulomb_G_recur(lam_min, locMax(N-k_lam_G,0), eta, x,
                    G_lam_min, Gp_lam_min,
		    &G_lam_G, &Gp_lam_G
		    );

    *F  = F_lam_F;
    *Fp = Fp_lam_F;
    *G  = G_lam_G;
    *Gp = Gp_lam_G;
    return GSL_SUCCESS;
  }
  else if(x < 2.0*eta) {
    /* Use WKB approximation to obtain F and G at the two
     * lambda values, and use the Wronskian and the
     * continued fractions for F'/F to obtain F' and G'.
     */
    double F_lam_F, G_lam_F;
    double F_lam_G, G_lam_G;
    double exp_lam_F, exp_lam_G;
    int stat_lam_F = coulomb_jwkb(lam_F, eta, x, &F_lam_F, &G_lam_F, &exp_lam_F);
    int stat_lam_G = coulomb_jwkb(lam_G, eta, x, &F_lam_G, &G_lam_G, &exp_lam_G);
    double Fp_over_F_lam_F;
    double Fp_over_F_lam_G;
    double F_sign_lam_F;
    double F_sign_lam_G;
    int stat_CF1_lam_F = coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F);
    int stat_CF1_lam_G = coulomb_CF1(lam_G, eta, x, &F_sign_lam_G, &Fp_over_F_lam_G);
    *F = F_lam_F;
    *G = G_lam_G;
    *Fp = Fp_over_F_lam_F * F_lam_F;
    *Gp = Fp_over_F_lam_G * G_lam_G - 1.0/F_lam_G;
    *exp_F = exp_lam_F;
    *exp_G = exp_lam_G;
    if(stat_CF1_lam_F != GSL_SUCCESS) {
      return stat_CF1_lam_F;
    }
    else if(stat_CF1_lam_G != GSL_SUCCESS) {
      return stat_CF1_lam_G;
    }
    else if(stat_lam_F == GSL_EOVRFLW || stat_lam_G == GSL_EOVRFLW) {
      return GSL_EOVRFLW;
    }
    else {
      return GSL_SUCCESS;
    }
  }
  else {
    /* x > 2 eta, so we know that we can find a lambda value such
     * that x is above the turning point. We do this, evaluate
     * using Steed's method at that oscillatory point, then
     * use recursion on F and G to obtain the required values.
     *
     * lam_0   = a value of lambda such that x is below the turning point
     * lam_min = minimum of lam_0 and the requested lam_G, since
     *           we must go at least as low as lam_G
     */
    const double SMALL = 1.0e-100;
    const double C = sqrt(1.0 + 4.0*x*(x-2.0*eta));
    const int N = ceil(lam_F - C + 0.5);
    const double lam_0   = lam_F - locMax(N, 0);
    const double lam_min = locMin(lam_0, lam_G);
    double F_lam_F, Fp_lam_F;
    double G_lam_G, Gp_lam_G;
    double F_lam_min_unnorm, Fp_lam_min_unnorm;
    double F_lam_min, Fp_lam_min;
    double G_lam_min, Gp_lam_min;
    double Fp_over_F_lam_F;
    double Fp_over_F_lam_min;
    double F_sign_lam_F;
    double P_lam_min, Q_lam_min;
    double alpha;
    double gamma;
    double F_scale;

    my_coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F);

    F_lam_F  = SMALL;
    Fp_lam_F = Fp_over_F_lam_F * F_lam_F;

    coulomb_F_recur(lam_min, locMax(k_lam_G, N), eta, x,
                    F_lam_F, Fp_lam_F,
		    &F_lam_min_unnorm, &Fp_lam_min_unnorm
		    );
    Fp_over_F_lam_min = Fp_lam_min_unnorm / F_lam_min_unnorm;

    coulomb_CF2(lam_min, eta, x, &P_lam_min, &Q_lam_min);

    alpha = Fp_over_F_lam_min - P_lam_min;
    gamma = alpha/Q_lam_min;
    F_lam_min  = F_sign_lam_F * sqrt(alpha*alpha/Q_lam_min + Q_lam_min);
    Fp_lam_min = Fp_over_F_lam_min * F_lam_min;
    G_lam_min  = gamma * F_lam_min;
    Gp_lam_min = P_lam_min * G_lam_min - Q_lam_min * F_lam_min;

    F_scale = F_lam_min / F_lam_min_unnorm;    
    F_lam_F  *= F_scale;
    Fp_lam_F *= F_scale;

    coulomb_G_recur(lam_min, locMax(N-k_lam_G,0), eta, x,
                    G_lam_min, Gp_lam_min,
		    &G_lam_G, &Gp_lam_G
		    );

    *F  = F_lam_F;
    *Fp = Fp_lam_F;
    *G  = G_lam_G;
    *Gp = Gp_lam_G;
    *exp_F = 0.0;
    *exp_G = 0.0;
    return GSL_SUCCESS;
  }
}



int
gsl_sf_coulomb_wave_F__array_impl(double lam_min, int kmax,
                                  double eta, double x, 
                                  double * fc,
			          double * F_exp)
{
  double G_exp;
}


int
gsl_sf_coulomb_wave_FG_array_impl(double lam_min, int kmax,
                                  double eta, double x,
                                  double * fc, double * gc,
			          double * F_exp, double * G_exp)
{
}

int
gsl_sf_coulomb_wave_FGp_array_impl(double lam_min, int kmax,
                                   double eta, double x,
		                   double * fc, double * fcp,
		                   double * gc, double * gcp,
			           double * F_exp, double * G_exp)

{
}


int
gsl_sf_coulomb_wave_sphF_array_impl(double lam_min, int kmax,
                                    double eta, double x,
		                    double * fc,
			            double * F_exp)
{
  int k;

  if(x < 0.0 || lam_min < -0.5) {
    return GSL_EDOM;
  }
  else if(x < 10.0/DBL_MAX) {
    for(k=0; k<=kmax; k++) {
      fc[k] = 0.0;
    }
    if(lam_min == 0.0) {
      fc[0] = sqrt(C0sq(eta));
    }
    *F_exp = 0.0;
    return GSL_SUCCESS;
  }
  else {
    double G_exp;
    /*
    for(k=0; k<=kmax; k++) {
      fc[k] = fc[k] / x;
    }
    */
    /* FIXME */
    return GSL_SUCCESS;
  }
}




int
gsl_sf_coulomb_wave_F_array_e(double lam_min, int kmax,
                        double eta, double x,
                        double * fc,
			double * F_exponent
                        )
{
/*
  int status = gsl_sf_coulomb_wave_F_array_impl(lam_min, kmax, x, eta,
                                          fc,
					  F_exponent
					  );
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_coulomb_wave_F_e", status);
  }
  return status;
  */
}


int
gsl_sf_coulomb_wave_FG_array_e(double lam_min, int kmax,
                         double eta, double x,
                         double * fc,
			 double * gc,
			 double * F_exponent,
			 double * G_exponent
                         )
{
/*
  int status = gsl_sf_coulomb_wave_FG_array_impl(lam_min, kmax,
                                           eta, x,
                                           fc, gc,
					   F_exponent, G_exponent
					   );
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_coulomb_wave_FG_e", status);
  }
  return status;
  */
}


int
gsl_sf_coulomb_wave_FGp_array_e(double lam_min, int kmax,
                          double eta, double x,
                          double * fc, double * fcp,
                          double * gc, double * gcp,
                          double * F_exponent,
                          double * G_exponent
                          )
{
/*
  int status = gsl_sf_coulomb_wave_FGp_array_impl(lam_min, kmax,
                                            eta, x,
                                            fc, fcp,
                                            gc, gcp,
                                            F_exponent,
                                            G_exponent
                                            );
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_coulomb_wave_FGp_e", status);
  }
  return status;
  */
}
