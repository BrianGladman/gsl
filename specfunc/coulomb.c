/* Author:  G. Jungman
 * RCS:     $Id$
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
 *   eta < 20
 *   x   < 20
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
    double rk2 = 1. + ek*ek;
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
  while(fabs(df-1.) > CF1_acc);
  
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


static
int
coulomb(double lam_min, int kmax,
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

    /* Obtain F'/F and sign(F) at the upper lambda value.
     */
    coulomb_CF1(lam_min, kmax, eta, x, &F_sign_lam_max, &F_ratio_lam_max);

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


/* junk for coulfg() */
#define CFG_ABORT	 1.e5 /* 2.e4 */
#define CFG_TM30	 1.e-30

#define CFG_ACCUR	 GSL_MACH_EPS
#define CFG_ACC	         (10.*CFG_ACCUR)
#define CFG_ACC4	 (CFG_ACC*100.*100.)


/* coulfg() mode control */
#define Mode_F   3
#define Mode_FG  2
#define Mode_FGp 1


/* zero argument calculation of coulomb wave functions
 * assumes xlmin >= 0.0
 *   fc[0] .. fc[kmax] = F_{xlmin}(eta,0) .. F_{xlmin+kmax}(eta,0)
 *   etc.
 */
static
int
coulF_zero_x(double eta, double xlmin, int kmax, double * fc)
{
  int k;
  for(k=0; k<=kmax; k++) fc[k] = 0.0;
  return GSL_SUCCESS;
}
static
int
coulG_zero_x(double eta, double xlmin, int kmax, double * gc)
{
  if(xlmin == 0.0 && kmax == 0) {
    gc[0] = 1.0/sqrt(C0sq(eta));
    return GSL_SUCCESS;
  }
  else {
    return GSL_EDOM;
  }
}
static
int
coulFp_zero_x(double eta, double xlmin, int kmax, double * fcp)
{
  int k;
  if(xlmin == 0.0) {
    fcp[0] = sqrt(C0sq(eta));
  }
  else {
    fcp[0] = 0.0;
  }
  for(k=1; k<=kmax; k++) fcp[k] = 0.0;
  return GSL_SUCCESS;
}
static
int
coulGp_zero_x(double eta, double xlmin, int kmax, double * gcp)
{
  if(xlmin == 0.0 && kmax == 0) {
    gcp[0] = 0.0;
    return GSL_SUCCESS;
  }
  else {
    return GSL_EDOM;
  }
}
static
int
coulfg_zero_x(double eta, double xlmin, int kmax,
              double * fc, double * gc, double * fcp, double * gcp,
              int mode)
{
  if(mode==Mode_F) {
    return coulF_zero_x(eta, xlmin, kmax, fc);
  }
  else if(mode==Mode_FG) {
    int status_F = coulF_zero_x(eta, xlmin, kmax, fc);
    int status_G = coulG_zero_x(eta, xlmin, kmax, gc);
    return status_F | status_G;
  }
  else if(mode==Mode_FGp) {
    int status_F  =  coulF_zero_x(eta, xlmin, kmax, fc);
    int status_G  =  coulG_zero_x(eta, xlmin, kmax, gc);
    int status_Fp = coulFp_zero_x(eta, xlmin, kmax, fcp);
    int status_Gp = coulGp_zero_x(eta, xlmin, kmax, gcp);
    return status_F | status_G | status_Fp | status_Gp;
  }
  else {
    /* NOT REACHED */
    return GSL_EFAILED;
  }
}


/* small argument calculation of coulomb wave functions 
 * based on expansion in terms of spherical Bessel functions
 * [Abramowitz and Stegun 14.4.5]
 *
 * This is effectively a series in eta*x, as long as x is small.
 *
 *   F_L  ~ C_L(eta) x^(L+1)
 *   Fp_L ~ (L+1) C_L(eta) x^L = (L+1)/x F_L
 *
 *   G_L  ~ 1/(2L+1) 1/C_L(eta) 1/x^L = 1/(2L+1) x/F_L,  L > 0
 *   G_0  ~ 2eta/C_0^2(eta) F_0 ln(2x) + 1/C_0(eta) ~ x/F_0(1 + 2eta x ln(2x))
 *
 *   Gp_L ~ -L/(2L+1) 1/C_L(eta) 1/x^(L+1) = -L/x G_L, L > 0
 *   Gp_0 ~ 2eta / C_0(eta) (ln(2x) + 1) ~ 2eta x/F_0 (ln(2x) + 1)
 *
 * assumes x > 0.0
 */
static
int
coulfg_small_args(double x, double eta, double xlmin, int kmax,
                  double * fc, double * gc, double * fcp, double * gcp,
                  int mode)
{
  int k;
  double * cl = (double *)malloc((kmax+1)*sizeof(double));

  if(cl == 0){
    return GSL_ENOMEM;
  }

  gsl_sf_coulomb_CL_list(xlmin, kmax, eta, cl);
 
  /* calculate F_L */
  if(mode==Mode_F || mode==Mode_FG || mode==Mode_FGp) {
    for(k=0; k<=kmax; k++) {
      fc[k] = cl[k] * pow(x,xlmin+k+1.0);
    }
  }
  
  /* calculate G_L */
  if(mode==Mode_FG || mode==Mode_FGp) {
    if(xlmin == 0.0) {
      gc[0] = x/fc[0] * (1.0 + 2.0*eta * x * log(2.0*x));
    }
    else {
      gc[0] = 1.0/(2.0*xlmin+1.0) /cl[0] /pow(x,xlmin);
    }
    for(k=1; k<=kmax; k++) {
      gc[k] = x/(2.0*(xlmin+k)+1.) / fc[k];
    }
  }

  /* calculate F_L' and G_L' */
  if(mode==Mode_FGp) {
    for(k=0; k<=kmax; k++) {
      fcp[k] = fc[k] * (xlmin+k+1.0)/x;
    }
    if(xlmin == 0.0) {
      gcp[0] = 2.0*eta * x/fc[0] * (1.0 + log(2.0*x));
    }
    else {
      gcp[0] = gc[0] * (-xlmin)/x;
    }
    for(k=1; k<=kmax; k++) {
      gcp[k] = -gc[k] * (xlmin + k)/x;
    }
  }

  free(cl);
  return GSL_SUCCESS;
}



/* ------------------------------------------------------------ 

 The following hack job is based on some fortran code.
   There is a small function jwkb() and the main function
   coulfg(). See the comments with coulfg().
   
  ------------------------------------------------------------ */


/* overflow exponent, can be != 0 if scaling required,
   which occurs if the wkb method is invoked and
   it generates large exponents
 */
static int over_exp_ = 0;
int coul_wave_overflow_exp(void) { return over_exp_; }


/*
  COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS FOR XL.GE. 0
  AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
  
  F_xl(eta, x)

*/
#define ALOGE  0.43429448
static
int
jwkb(double x, double eta, double xl,
     double * fjwkb, double * gjwkb, int *iexp)
{

  double gh2  = x*(2.*eta - x);
  double xll1 = locMax(xl*(xl + 1.), 0.0);

  if(gh2 + xll1 > 0.0) {   /* if( x < turning point) */
    double hll = xll1 + 6.0/35.0;
    double hl  = sqrt(hll);
    double sl  = eta/hl + hl/x;
    double rl2 = 1.0 + eta*eta/hll;
    double gh  = sqrt(gh2 + hll)/x;
    double phi = x*gh - 0.5*(hl*log((gh+sl)*(gh+sl)/rl2) - log(gh));
    double phi10;

    if(eta != 0.0) phi -= eta*atan2(x*gh,x - eta);
    phi10 = -phi*ALOGE;
    *iexp = (int)phi10;

    if(*iexp > 300)  {
      /* overflow in G: scale F,F' by 10^-iexp and G,G' by 10^iexp */
      over_exp_ = *iexp;
      *gjwkb = pow(10.0, phi10 - *iexp);
      return GSL_EOVRFLW;
    }
    if(*iexp <= 300) {
      *gjwkb = exp(-phi);
      *iexp  = 0;
    }
    *fjwkb = 0.5/(gh * *gjwkb);

    return GSL_SUCCESS;
  }
  else {
    /* above turning point: INTERNAL ERROR */
    return GSL_EFAILED;
  }
}
#undef ALOGE


/*
 SUBROUTINE COULFG(X,ETA,XLMIN,XLMAX, FC,GC,FCP,GCP, MODE1,KFN,IFAIL)

C  REVISED IJT WITH L-T ALGORITHMN FOR CONTINUED FRACTIONS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								       C
C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD	       C
C								       C
C  A. R. BARNETT	   MANCHESTER  MARCH   1981		       C
C								       C
C  ORIGINAL PROGRAM 'RCWFN'	 IN    CPC  8 (1974) 377-395	       C
C		  + 'RCWFF'	 IN    CPC 11 (1976) 141-142	       C
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314	       C
C  THIS VERSION WRITTEN UP	 IN    CPC 27 (1982) 147-166	       C
C								       C
C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
C   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    C
C   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   C
C   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     C
C   THE DIRAC EQUATION ,ALSO SPHERICAL & CYLINDRICAL BESSEL EQUATIONS  C
C								       C
C  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    C
C  STARTING ARRAY ELEMENT IS M1 = MAX0(IDINT(XLMIN+ACCUR),0) + 1       C
C      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES		       C
C								       C
C  IF 'MODE' = 1  GET F,G,F',G'	  FOR INTEGER-SPACED LAMBDA VALUES     C
C	     = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    C
C	     = 3      F		      CALL TO AT LEAST LENGTH (1)      C
C  IF 'KFN'  = 0 REAL	     COULOMB FUNCTIONS ARE RETURNED	       C
C	     = 1 SPHERICAL   BESSEL	 "	"     "		       C
C	     = 2 CYLINDRICAL BESSEL	 "	"     "		       C
C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT			       C
C								       C
C  PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    C
C   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     C
C   COULFG IS CODED FOR REAL*8 ON IBM OR EQUIVALENT  ACCUR = 10**-16   C
C   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   C
C   FOR MANTISSAS OF 56 & 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) C
C   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*/

/* I turned off the KFN input, since this will only be used
   for calculating Coulomb wave functions explicitly. Of course,
   there is no reason this cannot be called with eta=0, in which
   case Bessel functions is what you get...
   
   F_L(0, x) =  x j_L(x) =  sqrt(x) sqrt(pi/2) J_{L+1/2}(x)
   G_L(0, x) = -x y_L(x) = -sqrt(x) sqrt(pi/2) Y_{L+1/2}(x)
   
   ==>
   
   J_nu(x) =  sqrt(2/(pi x)) F_{nu-1/2}(0,x)
   Y_nu(x) = -sqrt(2/(pi x)) G_{nu-1/2}(0,x)

 */

/* This troglodyte code makes some sense if you read the paper:
   Barnett, Comp. Phys. Comm., 21, 297 (1981).
   But not much.
 */

/* Behavioural change: the original fortran version would load
   the computed values into the arrays starting at an index given
   by the lowest requested l value (integer part). This is dumb.
   I changed it so that it always loads them starting at index 0.
   You can easily get the other effect by passing an offset
   pointer if you like.
 */

static
int
coulfg(double x, double eta, double xlmin, int kmax,
       double * fc, double * gc, double * fcp, double * gcp,
       int mode)
{
  double acch = sqrt(CFG_ACC);

  double xlmax    = xlmin + kmax;
  double upper_l  = xlmax;
  int i_delta_lam = kmax;

  int iexp = 1;
  double fjwkb = 0.0;
  double gjwkb = 0.0;

  if(x < 0.0 || xlmin <= -0.5 || kmax < 0) {
    return GSL_EDOM;
  }

  if(x == 0.0) {
    return coulfg_zero_x(eta, xlmin, xlmax, fc, gc, fcp, gcp, mode);
  }
  else if(x < 0.001 && fabs(eta*x) < 0.01) {
    return coulfg_small_args(x, eta, xlmin, xlmax, fc, gc, fcp, gcp, mode);
  }
  else if(x < acch) {
    /* This is the case when eta is very large.
       There must be something easy we can do here.
    */
    return GSL_EFAILED;
  }
  else {
    /* check if x is below the turning point */
    int xlturn = (x*(x - 2.*eta) < xlmin*(xlmin + 1.0) ? 1 : 0);

    double e2mm1 = eta*eta + xlmin*(xlmin + 1.0);

    double x_inv = 1./x;
    double fcl = 1.0;
    double gcl;
    double pk  = upper_l + 1.0;
    double px  = pk  + CFG_ABORT;
    double df;

    int L;
    double xl;

    double F = eta/pk + pk*x_inv;
    double D;
    double C;

    /* the real ang imaginary parts
     * of the second continued fraction
     */
    double P;
    double Q;

    double gam;
    double fcm;
    double W;
    double gpl;

    if(fabs(F) < CFG_TM30) F = CFG_TM30;
    D = 0.0;
    C = F;

    /* Compute first continued fraction. This
     * gives F'/F at the upper lambda value.
     */
    do {
      double pk1 = pk + 1.0;
      double ek  = eta / pk;
      double rk2 = 1. + ek*ek;
      double tk  = (pk + pk1)*(x_inv + ek/pk1);
      D   =  tk - rk2 * D;
      C   =  tk - rk2 / C;
      if(fabs(C) < CFG_TM30) C = CFG_TM30;
      if(fabs(D) < CFG_TM30) D = CFG_TM30;
      D = 1.0/D;
      df = D * C;
      F  = F * df;
      if(D < 0.) fcl = -fcl; /* sign of result depends on sign of denominator */
      pk = pk1;
      if( pk > px ) {
        /* first continued fraction not converging */
	return GSL_ERUNAWAY;
      }
    }
    while(fabs(df-1.) > CFG_ACC);


    /* At this point we have computed F'/F at the
     * upper lambda value and stored it in F.
     *
     * Now use downward recurrence to the minimum lambda value.
     *    F_{lam-1}  = (S_lam F_lam + F_lam') / R_lam
     *    F_{lam-1}' = (S_lam F_{lam-1} - R_lam F_lam)
     * where
     *    R_lam = sqrt(1 + (eta/lam)^2)
     *    S_lam = lam/x + eta/lam
     *
     * Combining this with the ratio computed above, we
     * obtain F and F' at the minimum lambda value.
     *
     * Use the array gc[] to store the values of R_lam on
     * the way down, which we will need later if we are
     * calculating the G functions.
     */
    if(i_delta_lam != 0) {
      int lp;
      double fpl;
      fcl *= CFG_TM30;
      fpl = fcl*F;
      if(mode == Mode_FGp) fcp[i_delta_lam] = fpl;
      fc[i_delta_lam] = fcl;
      xl  = upper_l;
      for(lp=1; lp<=i_delta_lam; lp++) {
	double el = eta/xl;
	double rl = sqrt(1. + el*el);
	double sl = el  + xl*x_inv;
	double fc_lm1;
	fc_lm1 = (fcl*sl + fpl)/rl;
	fpl    =  fc_lm1*sl - fcl*rl;
	fcl    =  fc_lm1;
	L      =  i_delta_lam - lp;
	fc[L]  =  fcl;
	if(mode == Mode_FGp) fcp[L]  = fpl;
	if(mode != Mode_F  /* && eta_ne_zero */) gc[L+1] = rl;
	--xl;
      }
      
      if(fcl == 0.) fcl = CFG_ACC;
      F	= fpl/fcl;
    }


    /* At this point we have calculated the values
     * of F and F' at the minimum lambda:  fcl, fpl.
     * The unnormalized values of fcl[] and fcp[]
     * have also been set on the way down.
     */
    
    /* If we are below the turning point, there is a danger
     * that G is very large. Evaluate F and G at the minimum lambda,
     * using the WKB approximation.
     */
    if(xlturn) jwkb(x, eta, locMax(xlmin,0.0), &fjwkb, &gjwkb, &iexp);
    
    /* if(iexp != 1) fprintf(stderr,"iexp= %d\n", iexp); */
    /* FIXME: have to do something about this overflow business */

    if(iexp > 1 || gjwkb > 1.0/(acch*100.)) {
      /* ARRIVE HERE IF G(xlmin) > 10**6 OR iexp > 70 and xlturn = true
       */
      /* We do not attempt the second continued fraction
       * in this case since it will lose accuracy. Simply
       * set the normalization as determined by the WKB
       * solution and move on.
       */
      W	 = fjwkb;
      gam = gjwkb*W;
      P	  = F;
      Q	  = 1.;
    }
    else {
      /* Evaluate the second continued fraction to 
       * obtain the ratio
       *    (G' + i F')/(G + i F)
       * at the minimum lambda value.
       */
      double ta =  2.0*CFG_ABORT;
      double wi =  2.*eta;
      double ar = -e2mm1;
      double ai =  eta;
      double br =  2.*(x - eta);
      double bi =  2.;
      double dr =  br/(br*br + bi*bi);
      double di = -bi/(br*br + bi*bi);

      double dp = -x_inv*(ar*di + ai*dr);
      double dq =  x_inv*(ar*dr - ai*di);
      
      double A;
      double B;

      xlturn = 0;
      pk =  0.0;
      P	 =  0.0;
      Q	 =  1. - eta*x_inv;
      
      do {
	P += dp;
	Q += dq;
	pk += 2.;
	ar += pk;
	ai += wi;
	bi += 2.;
	D  = ar*dr - ai*di + br;
	di = ai*dr + ar*di + bi;
	C  = 1./(D*D + di*di);
	dr =  C*D;
	di = -C*di;
	A  = br*dr - bi*di - 1.;
	B  = bi*dr + br*di;
	C  = dp*A  - dq*B;
	dq = dp*B  + dq*A;
	dp = C;
	if(pk > ta) {
	  /* second continued fraction not converging */
	  return GSL_ERUNAWAY;
	}
      }
      while(fabs(dp)+fabs(dq) >= (fabs(P)+fabs(Q))*CFG_ACC);

      gam = (F - P)/Q;
      if(Q <= CFG_ACC4*fabs(P)) {
        /*
	GSL_ERROR("coulfg: final Q < abs(P)*CFG_ACC*1e4", GSL_EFAILED);
	*/
	return GSL_EFAILED;
	/* Not sure what this condition is */
      }
      
      /* This is the normalization factor */
      W	= 1.0/sqrt((F - P)*gam + Q);
    }
    
    /* Get the right sign for the value at the minimum lambda */
    fcm   = ( fcl > 0.0 ? W : -W);  /* fcm   = DSIGN(W,fcl)*beta; beta=1*/
    fc[0] = fcm;

    /* Calculate the other values at the
     * minimum lambda point, if needed.
     */
    if(mode == Mode_FG || mode == Mode_FGp) {
      if(! xlturn)   gcl =  fcm*gam;
      if(  xlturn)   gcl =  gjwkb;
      gc[0]  = gcl;
      gpl =  gcl*(P - Q/gam);
    }
    if(mode == Mode_FGp) {
      gcp[0] = gpl;
      fcp[0] = fcm*F;
    }
    
    /* If we are calculating for one lambda
     * value only, then we are done.
     */
    if(i_delta_lam == 0 ) return GSL_SUCCESS;

    /* Now we have the correct normalization factor. */
    W = W/fabs(fcl);

    /* Perform upward recursion to obtain the G
     * values, and update the F' and G' values
     * if those are needed. Recall that the R_lam
     * values from the first (downward) recursion are
     * stored in gc[].
     */
    for(L=0; L<=i_delta_lam-1; L++) {

      double gcl1;
      
      xl += 1.0;

      if(mode == Mode_FGp || mode == Mode_FG) {
	/* if(eta_ne_zero) */ double el = eta/xl;
	/* if(eta_ne_zero) */ double rl = gc[L+1];
	double sl = el + xl*x_inv;
	gcl1  = ((sl)*gcl - gpl)/rl;
	gpl   = rl*gcl - (sl)*gcl1;
	gcl   = gcl1;
	gc[L+1] = gcl1;
      }
      if(mode == Mode_FGp) {
	gcp[L+1] = gpl;
	fcp[L+1] = W*fcp[L+1];
      }
      fc[L+1] = W * fc[L+1];
    }
    
    return GSL_SUCCESS;
  }
}



int
gsl_sf_coulomb_wave_F_impl(double x, double eta,
		           double lam_min, double lam_max,
                           double * fc)
{
  return coulfg(x, eta, lam_min, lam_max,
	        fc, (double *) 0, (double *) 0, (double *) 0,
	        Mode_F
	        );
}


int
gsl_sf_coulomb_wave_FG_impl(double x, double eta,
                            double lam_min, double lam_max,
                            double * fc, double * gc)
{
  return coulfg(x, eta, lam_min, lam_max,
	        fc, gc, (double *) 0, (double *) 0,
	        Mode_FG
	        );
}

int
gsl_sf_coulomb_wave_FGp_impl(double x, double eta,
		             double lam_min, double lam_max,
		             double * fc, double * gc,
		             double * fcp, double * gcp)

{
  return coulfg(x, eta, lam_min, lam_max,
	        fc, gc, fcp, gcp,
	        Mode_FGp
	        );
}


int
gsl_sf_coulomb_wave_sphF_impl(double x, double eta,
		              double lam_min, double lam_max,
		              double * fc)
{
  int mode;
  double delta_lam = lam_max - lam_min + CFG_ACC;
  int i_delta_lam  = (int)delta_lam;
  int i;

  if(fabs(x) <= 100.*GSL_MACH_EPS) {
    /* if x==0 then it is nonzero and finite only for l=0 */
    double ell = lam_min;
    mode = 0;
    for(i=0; i<=i_delta_lam; i++) {
      if(fabs(ell) < 1000.*GSL_MACH_EPS) {
	fc[i] = sqrt(C0sq(eta));
      }
      else {
	fc[i] = 0.;
      }
      ell += 1.;
    }
  }
  else {
    mode = 1;
    coulfg(x, eta, lam_min, lam_max,
	   fc, (double *) 0, (double *) 0, (double *) 0,
	   Mode_F
	   );
    for(i=0; i<=i_delta_lam; i++) { fc[i] = fc[i] / x; }
  }
  return GSL_SUCCESS;
}


int
test_coulomb(void)
{
  const int kmax = 100;
  int k;

  double lam_min = 0.0;
  double eta = 10.0;
  double x = 20.0;

  double fc[kmax+1], fcp[kmax+1], gc[kmax+1], gcp[kmax+1];
  double F_e, G_e;

  int stat = coulomb(lam_min, kmax, eta, x, fc, fcp, gc, gcp, &F_e, &G_e);
  
  for(k=0; k<=kmax; k++) {
    printf("%5.3g   %16.12g  %16.12g  %16.12g  %16.12g\n",
           lam_min + k, fc[k], fcp[k], gc[k], gcp[k]
	   );
  }

  return stat;
}
