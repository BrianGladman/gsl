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
#include "gsl_sf_pow_int.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_coulomb.h"

#define locMax(a, b) ((a) > (b) ? (a) : (b))


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
    return twopieta/expm1(twopieta);
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
coulomb_Phi_series(double lam, double eta, double x,
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


/* Evolve the backwards recurrence for F,F'.
 *
 *    F_{lam-1}  = (S_lam F_lam + F_lam') / R_lam
 *    F_{lam-1}' = (S_lam F_{lam-1} - R_lam F_lam)
 * where
 *    R_lam = sqrt(1 + (eta/lam)^2)
 *    S_lam = lam/x + eta/lam
 *
 * The starting values F,Fp are given.
 * The array gc[] is preloaded with the
 * values of R_lam for use later, as in
 * the Barnett COULFG() code.
 */
static
int
coulomb_F_recur(double lam_min, int kmax,
                double eta, double x,
                double F_lam_max, double Fp_lam_max,
		double * F_lam_min, double * Fp_lam_min,
                double * fc, double * fcp,
                double * gc
                )
{
  double x_inv = 1.0/x;
  double fcl = F_lam_max;
  double fpl = Fp_lam_max;
  double lam_max = lam_min + kmax;
  double lam = lam_max;
  int k;

  fc[kmax] = fcl;
  if(fcp != (double *)0) fcp[kmax] = fpl;

  for(k=kmax-1; k>=0; k--) {
    double el = eta/lam;
    double rl = sqrt(1. + el*el);
    double sl = el  + lam*x_inv;
    double fc_lm1;
    fc_lm1 = (fcl*sl + fpl)/rl;
    fpl    =  fc_lm1*sl - fcl*rl;
    fcl    =  fc_lm1;
    fc[k]  =  fcl;
    if(fcp != (double *)0) fcp[k]   = fpl;
    if(gc  != (double *)0)  gc[k+1] = rl;
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
 * It is assumed that the values of R_lam are stored
 * in the gc[] array upon entry, as discussed above.
 */
static
int
coulomb_G_recur(double lam_min, int kmax,
                double eta, double x,
                double G, double Gp,
                double * gc, double * gcp
                )
{
  double x_inv = 1.0/x;
  double gcl = G;
  double gpl = Gp;
  double lam = lam_min + 1.0;
  int k;

  for(k=1; k<=kmax; k++) {
    double el = eta/lam;
    double rl = gc[k];
    double sl = el + lam*x_inv;
    double gcl1 = (sl*gcl - gpl)/rl;
    gpl   = rl*gcl - sl*gcl1;
    gcl   = gcl1;
    gc[k] = gcl1;
    if(gcp != (double *)0) gcp[k] = gpl;
    lam += 1.0;
  }
  
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
my_coulomb_CF1(double lam_min, int kmax,
            double eta, double x,
            double * fcl_sign,
	    double * result
            )
{
  const double CF1_small = 1.e-30;
  const double CF1_abort = 1.e5;
  const double CF1_acc   = 10.0*GSL_MACH_EPS;
  const double x_inv     = 1.0/x;
  const double lam_max   = lam_min + kmax;
  const double px        = lam_max + 1.0 + CF1_abort;
  
  double pk = lam_max + 1.0;
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
      return GSL_ERUNAWAY;
    }
  }
  while(fabs(df-1.0) > CF1_acc);
  
  *result = F;
  return GSL_SUCCESS;
}

static
int
coulomb_CF1(double lam_min, int kmax,
            double eta, double x,
            double * fcl_sign,
	    double * result
            )
{
  const double CF1_small = 1.e-30;
  const double CF1_abort = 1.e5;
  const double CF1_acc   = 10.0*GSL_MACH_EPS;
  const double x_inv     = 1.0/x;
  const double lam_max   = lam_min + kmax;
  const double px        = lam_max + 1.0 + CF1_abort;
  
  double pk = lam_max + 1.0;
  
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
 * at the minimum lambda value.
 */
static
int
coulomb_CF2(double lam_min, double eta, double x,
            double * result_P, double * result_Q
            )
{
  int status = GSL_SUCCESS;

  const double CF2_acc   = 10.0*GSL_MACH_EPS;
  const double CF2_abort = 1.0e5;

  const double wi    = 2.0*eta;
  const double x_inv = 1.0/x;
  const double e2mm1 = eta*eta + lam_min*(lam_min + 1.0);
  
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
  while(fabs(dp)+fabs(dq) >= (fabs(P)+fabs(Q))*CF2_acc);

  if(Q <= CF2_abort*GSL_MACH_EPS*fabs(P)) {
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
  double xll1 = locMax(lam*(lam + 1.0), 0.0);
  double hll  = xll1 + 6.0/35.0;
  double hl   = sqrt(hll);
  double sl   = eta/hl + hl/x;
  double rl2  = 1.0 + eta*eta/hll;
  double gh   = sqrt(gh2 + hll)/x;
  double phi  = x*gh - 0.5*(hl*log((gh+sl)*(gh+sl)/rl2) - log(gh));

  if(eta != 0.0) phi -= eta*atan2(x*gh,x - eta);

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

/* FIXME: Does jwkb() really return F and G, or
 * are they scaled in some sick way? The fortran
 * is totally confusing.
 * It would make a difference only here, where
 * we try to use the WKB thing.
 */
/* Determine the values of F and G at the
 * minimum lambda, given the non-normalized
 * values of F,F' and the sign of F there.
 *
 * If G is very large (and consequently F is very small),
 * we signal GSL_EOVRFLW and set an exponent such that:
 *
 *   result_F = F_lam_min * exp(F_exponent)
 *   result_G = G_lam_min * exp(G_exponent)
 */
static
int
coulomb_lam_min_values(double lam_min,
                       double eta, double x,
		       double F_sign,
                       double * F_lam_min,   double * G_lam_min,
                       double * Fp_lam_min,  double * Gp_lam_min,
		       double * F_exponent,  double * G_exponent
                       )
{
  int x_less_turn = (x*(x - 2.*eta) < lam_min*(lam_min+ 1.0) ? 1 : 0);
  int stat_WKB;
  int stat_CF2;

  double P, Q;
  double fjwkb, gjwkb;
  double exponent;

  *F_exponent = 0.0;
  *G_exponent = 0.0;

  if(x_less_turn) {
    double lam = locMax(lam_min, 0.0);
    stat_WKB = coulomb_jwkb(lam, eta, x, &fjwkb, &gjwkb, &exponent);
    if(stat_WKB == GSL_EOVRFLW) {
      /* WKB detects a large G value. Do not attempt
       * the continued fraction.
       */
      *F_lam_min  = fjwkb;
      *G_lam_min  = gjwkb;
      *F_exponent = -exponent;
      *G_exponent =  exponent;
      return GSL_EOVRFLW;
    }
  }

  /* WKB did not detect an overflow, so we are not very far
   * from the turning point, or we are above the
   * turning point. Attempt the continued fraction.
   */
  stat_CF2 = coulomb_CF2(lam_min, eta, x, &P, &Q);

  if(stat_CF2 == GSL_SUCCESS) {
    if(*F_lam_min != 0.0) {
      double F_ratio = *Fp_lam_min / *F_lam_min;
      double gamma = (F_ratio - P)/Q;
      double omega = 1.0/sqrt((F_ratio - P)*gamma + Q);
      *F_lam_min  =  F_sign * omega;
      *G_lam_min  = *F_lam_min * gamma;
      *Fp_lam_min = *F_lam_min * F_ratio;
      *Gp_lam_min = *G_lam_min * (P - Q/gamma);
      return GSL_SUCCESS;
    }
    else {
      double Fp_sign_lam_min = (*Fp_lam_min > 0.0 ? 1.0 : -1.0);
      *Fp_lam_min = Fp_sign_lam_min * sqrt(fabs(Q));
      *G_lam_min  = *Fp_lam_min / Q;
      *Gp_lam_min = *G_lam_min * P;
      return GSL_SUCCESS;
    }
  }
  else {
    if(x_less_turn) {
      /* Fallback on WKB values if CF2 claims it failed. */
      *F_lam_min  = fjwkb;
      *G_lam_min  = gjwkb;
      return GSL_SUCCESS;
    }
    else {
      /* This is bad news. We are int the oscillating
       * region, but CF2 claims it failed. This should
       * not happen.
       */
      return GSL_EFAILED;
    }
  }
}


/* Evaluate the F and G functions in the
 * small argument regime, where we can use
 * the series representation of F,F' to
 * avoid the continued fraction.
 */
static
int
coulomb_small_args(double lam_min, int kmax,
                   double eta, double x,
		   double * fc, double * fcp,
		   double * gc, double * gcp
		   )
{
  double CL, r;
  double P, Q;
  double F_lam_max, Fp_lam_max;
  double F_lam_min, Fp_lam_min;
  double G_lam_min, Gp_lam_min;
  double lam_max = lam_min + kmax; 
  
  int stat_CF2 = GSL_SUCCESS;
  int stat_Phi;
  int stat_Fre;
  
  /* explicitly evaluate F,F' at lam_max */
  stat_Phi = coulomb_Phi_series(lam_max, eta, x, &F_lam_max, &Fp_lam_max);
  gsl_sf_coulomb_CL_impl(lam_max, eta, &CL);
  r = CL * gsl_sf_pow_int(x, lam_max + 1.0);
  F_lam_max  *= r;
  Fp_lam_max *= r/x;
  
  /* backward recurrence to get values at lam_min and fill fc[] and fcp[] */
  stat_Fre = coulomb_F_recur(lam_min, kmax, eta, x,
                             F_lam_max, Fp_lam_max,
                             &F_lam_min, &Fp_lam_min,
                             fc, fcp, gc);
  
  if(gc != (double *)0) {
  
    /* continued fraction to get the ratio (G' + i F')/(G + i F) */
    stat_CF2 = coulomb_CF2(lam_min, eta, x, &P, &Q);
  
    /* solve for the values of G,G' at lam_min */
    G_lam_min  = (Fp_lam_min - P * F_lam_min)/Q;
    Gp_lam_min = P*G_lam_min - Q*F_lam_min;
    
    /* forward recurrence to get gc[] and gcp[] */
    coulomb_G_recur(lam_min,  kmax, eta, x, G_lam_min, Gp_lam_min, gc, gcp);
  }
  
  if(stat_Phi != GSL_SUCCESS) {
    return stat_Phi;
  }
  else if(stat_CF2 != GSL_SUCCESS) {
    return stat_CF2;
  }
  else {
    return GSL_SUCCESS;
  }
}


/* Evaluate sequences of Coulomb
 * wave functions for x = 0.
 */
static
int
coulomb_zero_x(double lam_min, int kmax,
               double eta,
               double * fc, double * fcp,
	       double * gc, double * gcp
               )
{
  const double C0 = sqrt(C0sq(eta));
  int status = GSL_SUCCESS;
  int k;

  if(fc != (double *)0) {
    for(k=0; k<=kmax; k++) fc[k] = 0.0;
  }

  if(fcp != (double *)0) {
    for(k=0; k<=kmax; k++) fcp[k] = 0.0;
    if(lam_min == 0.0) fcp[0] = C0;
  }

  if(gc != (double *)0) {
    if(lam_min == 0.0 && kmax == 0) {
      gc[0] = 1.0/C0;
    }
    else {
      status = GSL_EDOM;
    }
  }

  if(gcp != (double *)0) {
    if(lam_min == 0.0 && kmax == 0) {
      gcp[0] = 0.0;
    }
    else {
      status = GSL_EDOM;
    }
  }
  
  return status;
}


/* Coulomb wave functions F_L(eta,x), G_L(eta,x) and
 * their derivatives. Apply Steed's method in the
 * standard regime. Catch asymptotic regimes and
 * apply series methods when possible.
 *
 * Note the definition is such that:
 *     F_L(0, x) =  x j_L(x) =  sqrt(x) sqrt(pi/2) J_{L+1/2}(x)
 *     G_L(0, x) = -x y_L(x) = -sqrt(x) sqrt(pi/2) Y_{L+1/2}(x)
 *  
 *  ==>
 *  
 *     J_nu(x) =  sqrt(2/(pi x)) F_{nu-1/2}(0,x)
 *     Y_nu(x) = -sqrt(2/(pi x)) G_{nu-1/2}(0,x)
 *
 */
/* Original fortran code from COULFG() subroutine, heavily
 * modified by yours truly. The troglodyte fortran code makes
 * some sense if you read the paper:
 * [Barnett, Comp. Phys. Comm., 21, 297 (1981)]
 * But not much.
 */
static
int
gsl_sf_coulomb_wave_impl(double lam_min, int kmax,
                         double eta, double x,
                         double * fc, double * fcp,
	                 double * gc, double * gcp,
	                 double * F_exponent,
	                 double * G_exponent
                         )
{
  const double small = 100.0*GSL_SQRT_DBL_MIN;

  if(fc == (double *)0) {
    return GSL_EFAULT;
  }
  if(gc == (double *)0 && gcp != (double *)0) {
    return GSL_EFAULT;
  }

  if(x < 0.0 || lam_min <= -0.5 || kmax < 0) {
    return GSL_EDOM;
  }

  if(x == 0.0) {
    return coulomb_zero_x(lam_min, kmax, eta, fc, fcp, gc, gcp);
  }
  else if(x < 2.0 && fabs(eta) < 10.0) {
    return coulomb_small_args(lam_min, kmax, eta, x, fc, fcp, gc, gcp);
  }
  else {
    int k;
    double F_sign_lam_max = 1.0;
    double F_ratio_lam_max;
    double F_lam_max, Fp_lam_max;
    double F_lam_min, Fp_lam_min;
    double G_lam_min, Gp_lam_min;
    double Fp_unnorm_lam_min;
    double F_scale;
    int stat_CF1;

    /* Obtain F'/F and sign(F) at the upper lambda value.
     */
    stat_CF1 = coulomb_CF1(lam_min, kmax, eta, x,
                           &F_sign_lam_max, &F_ratio_lam_max
                           );

    /* Evolve the downward recurrence to the minimum lambda value.
     * Obtain "F_lam_min" and "Fp_lam_min", which do not have
     * the correct overall normalization. Note however that
     * the signs are correct because we have applied the sign
     * of F above, so the overall factor is positive.
     * Arrays fc[],fcp[] are filled similarly.
     */
    F_lam_max  = F_sign_lam_max  * small;
    Fp_lam_max = F_ratio_lam_max * F_lam_max;
    coulomb_F_recur(lam_min, kmax, eta, x,
                    F_lam_max, Fp_lam_max,
		    &F_lam_min, &Fp_lam_min, 
                    fc, fcp, gc);
    Fp_unnorm_lam_min = Fp_lam_min;

    /* Obtain the properly normalized values at the minimum lambda. */
    coulomb_lam_min_values(lam_min, eta, x,
                           F_sign_lam_max,
                           &F_lam_min,  &G_lam_min,
                           &Fp_lam_min, &Gp_lam_min,
                           F_exponent,  G_exponent);

    /* Determine the scaling factor for F,F'. */
    if(fc[0] != 0) {
      F_scale = F_lam_min / fc[0];
    }
    else {
      F_scale = Fp_lam_min / Fp_unnorm_lam_min;
    }

    /* Apply scaling to F,F'. */
    for(k=0; k<=kmax; k++) fc[k] *= F_scale;
    if(fcp != (double *)0) {
      for(k=0; k<=kmax; k++) fcp[k] *= F_scale;
    }

    /* Evolve upward recurrence for G and G'. */
    if(gc  != (double *)0)  gc[0] =  G_lam_min;
    if(gcp != (double *)0) gcp[0] = Gp_lam_min;
    if(gc  != (double *)0) {
      coulomb_G_recur(lam_min, kmax, eta, x, G_lam_min, Gp_lam_min, gc, gcp);
    }

    return GSL_SUCCESS;
  }
}


int
gsl_sf_coulomb_wave_F_impl(double lam_min, int kmax,
                           double eta, double x, 
                           double * fc,
			   double * F_exp)
{
  double G_exp;
  return gsl_sf_coulomb_wave_impl(lam_min, kmax,
                                  eta, x,
	                          fc, (double *) 0,
                                  (double *) 0, (double *) 0,
				  F_exp, &G_exp
	                          );
}


int
gsl_sf_coulomb_wave_FG_impl(double lam_min, int kmax,
                            double eta, double x,
                            double * fc, double * gc,
			    double * F_exp, double * G_exp)
{
  return gsl_sf_coulomb_wave_impl(lam_min, kmax,
                                  eta, x,
                                  fc, (double *) 0,
                                  gc, (double *) 0,
				  F_exp, G_exp
                                  );
}

int
gsl_sf_coulomb_wave_FGp_impl(double lam_min, int kmax,
                             double eta, double x,
		             double * fc, double * fcp,
		             double * gc, double * gcp,
			     double * F_exp, double * G_exp)

{
  return gsl_sf_coulomb_wave_impl(lam_min, kmax,
                                  eta, x, 
                                  fc, fcp,
                                  gc, gcp,
				  F_exp, G_exp
                                  );
}


int
gsl_sf_coulomb_wave_sphF_impl(double lam_min, int kmax,
                              double eta, double x,
		              double * fc,
			      double * F_exp)
{
  int k;

  if(x < 0.0 || lam_min < -0.5) {
    return GSL_EDOM;
  }
  else if(x < 1.0/DBL_MAX) {
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
    gsl_sf_coulomb_wave_impl(lam_min, kmax,
                             eta, x,
                             fc, (double *) 0,
                             (double *) 0, (double *) 0,
			     F_exp, &G_exp
                             );
    for(k=0; k<=kmax; k++) {
      fc[k] = fc[k] / x;
    }
    return GSL_SUCCESS;
  }
}


int
gsl_sf_coulomb_wave_F_e(double lam_min, int kmax,
                        double eta, double x,
                        double * fc,
			double * F_exponent
                        )
{
  int status = gsl_sf_coulomb_wave_F_impl(lam_min, kmax, x, eta,
                                          fc,
					  F_exponent
					  );
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_coulomb_wave_F_e", status);
  }
  return status;
}


int
gsl_sf_coulomb_wave_FG_e(double lam_min, int kmax,
                         double eta, double x,
                         double * fc,
			 double * gc,
			 double * F_exponent,
			 double * G_exponent
                         )
{
  int status = gsl_sf_coulomb_wave_FG_impl(lam_min, kmax,
                                           eta, x,
                                           fc, gc,
					   F_exponent, G_exponent
					   );
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_coulomb_wave_FG_e", status);
  }
  return status;
}


int
gsl_sf_coulomb_wave_FGp_e(double lam_min, int kmax,
                          double eta, double x,
                          double * fc, double * fcp,
                          double * gc, double * gcp,
                          double * F_exponent,
                          double * G_exponent
                          )
{
  int status = gsl_sf_coulomb_wave_FGp_impl(lam_min, kmax,
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
}


int
test_coulomb(void)
{
  const int kmax = 3;
  int k;

  double lam_min = 0.0;
  double eta = -50.0;
  double x = 5.0;

  double fc[kmax+1], fcp[kmax+1], gc[kmax+1], gcp[kmax+1];
  double F_e, G_e;

  int stat = gsl_sf_coulomb_wave_impl(lam_min, kmax,
                                      eta, x,
                                      fc, fcp,
				      gc, gcp,
				      &F_e, &G_e);
  
  for(k=0; k<=kmax; k++) {
    printf("%5.3g   %20.16g  %20.16g  %20.16g  %20.16g  %s\n",
           lam_min + k, fc[k], fcp[k], gc[k], gcp[k], gsl_strerror(stat)
	   );
  }

  return stat;
}
