/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "hyperg.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_bessel.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_hyperg.h"

#define locMAX(a,b)     ((a) > (b) ? (a) : (b))
#define locEPS          (1000.0*GSL_MACH_EPS)


/* 1F1(a,b,x)/Gamma(b) for b->0
 * [limit of Abramowitz+Stegun 13.3.7]
 */
static
int
hyperg_1F1_renorm_b0(const double a, const double x, double * result)
{
  double eta = a*x;
  if(eta > 0.0) {
    double root_eta = sqrt(eta);
    double I1_scaled;
    int stat_I = gsl_sf_bessel_I1_scaled_impl(2.0*root_eta, &I1_scaled);
    if(stat_I != GSL_SUCCESS) {
      *result = 0.0;
      return stat_I;
    }
    else {
      double lnr = 0.5*x + 0.5*log(eta) + fabs(x) + log(I1_scaled);
      return gsl_sf_exp_impl(lnr, result);
    }
  }
  else if(eta == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    double root_eta = sqrt(-eta);
    double J1;
    int stat_J = gsl_sf_bessel_J1_impl(2.0*root_eta, &J1);
    if(stat_J != GSL_SUCCESS) {
      *result = 0.0;
      return stat_J;
    }
    else {
      double lnr = 0.5*x + 0.5*log(-eta) + fabs(x) + log(J1);
      double ex;
      int stat_e = gsl_sf_exp_impl(lnr, &ex);
      *result = -ex;
      return stat_e;
    }
  }
  
}


/* 1F1'(a,b,x)/1F1(a,b,x)
 * Uses Gautschi's version of the CF.
 * [Gautschi, Math. Comp. 31, 994 (1977)]
 *
 * Supposedly this suffers from the "anomalous convergence"
 * problem when b < x. I have seen anomalous convergence
 * in several of the continued fractions associated with
 * 1F1(a,b,x). This particular CF formulation seems stable
 * for b > x. However, it does display a painful artifact
 * of the anomalous convergence; the convergence plateaus
 * unless b >>> x. For example, even for b=1000, x=1, this
 * method locks onto a ratio which is only good to about
 * 4 digits. Apparently the rest of the digits are hiding
 * way out on the plateau, but finite-precision lossage
 * means you will never get them.
 */
static
int
hyperg_1F1_CF1_p(const double a, const double b, const double x, double * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int n = 1;
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = 1.0;
  double b1 = 1.0;
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
    an = (a+n)*x/((b-x+n-1)*(b-x+n));
    bn = 1.0;
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
    
    if(fabs(del - 1.0) < 10.0*GSL_MACH_EPS) break;
  }

  *result = a/(b-x) * fn;

  if(n == maxiter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* 1F1'(a,b,x)/1F1(a,b,x)
 * Uses Gautschi's series transformation of the
 * continued fraction. This is apparently the best
 * method for getting this ratio in the stable region.
 * The convergence is monotone and supergeometric
 * when b > x.
 */
static
int
hyperg_1F1_CF1_p_ser(const double a, const double b, const double x, double * result)
{
  const int maxiter = 5000;
  double sum  = 1.0;
  double pk   = 1.0;
  double rhok = 0.0;
  int k;
  for(k=1; k<maxiter; k++) {
    double ak = (a + k)*x/((b-x+k-1.0)*(b-x+k));
    rhok = -ak*(1.0 + rhok)/(1.0 + ak*(1.0+rhok));
    pk  *= rhok;
    sum += pk;
    if(fabs(pk/sum) < 2.0*GSL_MACH_EPS) break;
  }
  *result = a/(b-x) * sum;
  if(k == maxiter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* 1F1(a+1,b,x)/1F1(a,b,x)
 *
 * I think this suffers from typical "anomalous convergence".
 * I could not find a region where it was truly useful.
 */
static
int
hyperg_1F1_CF1(const double a, const double b, const double x, double * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int n = 1;
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = b - a - 1.0;
  double b1 = b - x - 2.0*(a+1.0);
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
    an = (a + n - 1.0) * (b - a - n);
    bn = b - x - 2.0*(a+n);
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
    
    if(fabs(del - 1.0) < 10.0*GSL_MACH_EPS) break;
  }

  *result = fn;
  if(n == maxiter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* 1F1(a,b+1,x)/1F1(a,b,x)
 *
 * This seemed to suffer from "anomalous convergence".
 * However, I have no theory for this recurrence.
 */
static
int
hyperg_1F1_CF1_b(const double a, const double b, const double x, double * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int n = 1;
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = b + 1.0;
  double b1 = (b + 1.0) * (b - x);
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
    an = (b + n) * (b + n - 1.0 - a) * x;
    bn = (b + n) * (b + n - 1.0 - x);
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
    
    if(fabs(del - 1.0) < 10.0*GSL_MACH_EPS) break;
  }

  *result = fn;
  if(n == maxiter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* We define two scaled functions for general use.
 * One is used for a > 0, the other for a < 0.
 *
 * a > 0, a = ap + N: 
 *   Y(ap,N,b,x) := Gamma[ap+N]/Gamma[1+ap+N-b] 1F1(ap+N, b, x)
 *
 * a < 0, a = ap - N: 
 *   Z(ap,N,b,x) := 1/Gamma(1+ap-N-b) 1F1(ap-N,b,x)
 *
 * Recursion relations will be written in terms of
 * these. This was convenient since the discussion
 * in Temme's paper on U(a,b,x) is in terms of 
 * these recursions. We make use of some of his
 * statements about the dominant solutions of these.
 */


/* Determine 1F1(a,b,x), 1F1(a+1,b,x)
 * from Y(a,0,b,x) and Y(a,1,b,x).
 */
static
int
hyperg_1F1_from_Y(const double a, const double b, const double x,
                  const double Y0, const double Y1,
                  double * Fa, double * Fap1)
{
  double lg_a, sg_a;
  double lg_ab, sg_ab;
  int stat_a  = gsl_sf_lngamma_sgn_impl(a,       &lg_a,  &sg_a);
  int stat_ab = gsl_sf_lngamma_sgn_impl(1.0+a-b, &lg_ab, &sg_ab);
  if(stat_a == GSL_SUCCESS && stat_ab == GSL_SUCCESS) {
    /* FIXME: handle overflows properly */
    double pre = sg_a * sg_ab * exp(lg_ab - lg_a);
    *Fa   = pre * Y0;
    *Fap1 = pre * (1.0+a-b)/a * Y1;
    return GSL_SUCCESS;
  }
  else {
    *Fa   = 0.0;
    *Fap1 = 0.0;
    return GSL_EDOM;
  }
}


/* Upward recursion on the parameter 'a', which
 * a forward recursion for the function Y(a,n,b,x).
 *
 *     (a+n+1-b)Y(n+1) + (b-2a-2n-x)Y(n) + (a+n-1)Y(n-1) = 0
 *
 * Note that this will bomb if a-b+1 is a negative integer and
 * the recursion passes through that integer.
 */
static
int
hyperg_1F1_Y_recurse_up(double a, double b, double x,
                        int n, double Ynm1, double Yn,
		        int N,
		        double * YNm1, double * YN)
{
  int k;
  double Ykm1 = Ynm1;
  double Yk   = Yn;
  double Ykp1;

  for(k=n; k<N; k++) {
    double ckp1 = (a + k + 1 - b);
    double ck   = b - 2*a - 2*k - x;
    double ckm1 = a + k - 1;
    
    if(fabs(ckp1) < GSL_MACH_EPS) {
      *YN = 0.0;
      return GSL_EDOM;
    }
    Ykp1 = (-ck*Yk - ckm1*Ykm1)/ckp1;
    
    Ykm1 = Yk;
    Yk   = Ykp1;
  }
  
  *YNm1 = Ykm1;
  *YN   = Yk;
  return GSL_SUCCESS;
}


/* Evaluate Y(ap,N+1,b,x)/Y(ap,N,b,x) by Steed's continued fraction method.
 * As usual, this is used to pick out the minimal solution of the
 * associated recursion relation.
 */
static
int
hyperg_1F1_Y_CF1(const double ap, const double b, const int N, const double x, double * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int n = 1;
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = -(ap + N);
  double b1 =  (b - 2.0*ap - x - 2.0*(N+1));
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
    an = -(ap + N + n - b)*(ap + N + n - 1.0);
    bn =  (b - 2.0*ap - x - 2.0*(N+n));
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
    
    if(fabs(del - 1.0) < 10.0*GSL_MACH_EPS) break;
  }

  *result = fn;
  if(n == maxiter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* Recursion in the negative a direction. This is
 * forward recursion on the function Z(a, n, b, x).
 *
 * Z[n+1] + (b-2a-x+2n) Z[n] + (a-n)(1+a-n-b) Z[n-1] = 0
 */
static
int
hyperg_1F1_Z_recurse(const double a, const double b, const double x,
                     int n, double Znm1, double Zn,
		     int N,
		     double * ZNm1, double * ZN)
{
  int k;
  double Zkm1 = Znm1;
  double Zk   = Zn;
  double Zkp1;

  for(k=n; k<N; k++) {
    Zkp1 = -(b-2.0*a-x+2.0*k) * Zk - (a-k)*(1.0+a-k-b) * Zkm1;
    Zkm1 = Zk;
    Zk   = Zkp1;
  }
  
  *ZNm1 = Zkm1;
  *ZN   = Zk;
  return GSL_SUCCESS;
}


/* Downward recursion on the parameter 'a'. This
 * is backward recursion for the function Y(a,n,b,x).
 *
 *     (a+n+1-b)Y(n+1) + (b-2a-2n-x)Y(n) + (a+n-1)Y(n-1) = 0
 *
 * Note that this will bomb if a-b+1 is a negative integer and
 * the recursion passes through that integer.
 * It will also bomb if a=0 and N=0.
 */
static
int
hyperg_1F1_Y_recurse_down(double a, double b, double x,
                          int n, double Ynp1, double Yn,
		          int N,
		          double * YNp1, double * YN)
{
  int k;
  double Ykp1 = Ynp1;
  double Yk   = Yn;
  double Ykm1;

  for(k=n; k>N; k--) {
    double ckp1 = (a + k + 1 - b);
    double ck   = b - 2*a - 2*k - x;
    double ckm1 = a + k - 1;
    
    if(fabs(ckm1) < GSL_MACH_EPS) {
      *YN = 0.0;
      return GSL_EDOM;
    }
    Ykm1 = (-ck*Yk - ckp1*Ykp1)/ckm1;
    
    Ykp1 = Yk;
    Yk   = Ykm1;
  }

  *YNp1 = Ykp1;
  *YN   = Yk;
  return GSL_SUCCESS;
}


/* Luke's rational approximation.
 * See [Luke, Algorithms for the Computation of Mathematical Functions, p.182]
 *
 * Like the case of the 2F1 rational approximations, these are
 * probably guaranteed to converge for x < 0, barring gross
 * numerical instability in the pre-asymptotic regime.
 */
static
int
hyperg_1F1_luke(const double a, const double c, const double xin,
                double * result, double * prec)
{
  const double RECUR_BIG = 1.0e+50;
  const int nmax = 5000;
  int n = 3;
  const double x  = -xin;
  const double x3 = x*x*x;
  const double t0 = a/c;
  const double t1 = (a+1.0)/(2.0*c);
  const double t2 = (a+2.0)/(2.0*(c+1.0));
  double F = 1.0;

  double Bnm3 = 1.0;                                  /* B0 */
  double Bnm2 = 1.0 + t1 * x;                         /* B1 */
  double Bnm1 = 1.0 + t2 * x * (1.0 + t1/3.0 * x);    /* B2 */
 
  double Anm3 = 1.0;                                                      /* A0 */
  double Anm2 = Bnm2 - t0 * x;                                            /* A1 */
  double Anm1 = Bnm1 - t0*(1.0 + t2*x)*x + t0 * t1 * (c/(c+1.0)) * x*x;   /* A2 */

  while(1) {
    double npam1 = n + a - 1;
    double npcm1 = n + c - 1;
    double npam2 = n + a - 2;
    double npcm2 = n + c - 2;
    double tnm1  = 2*n - 1;
    double tnm3  = 2*n - 3;
    double tnm5  = 2*n - 5;
    double F1 =  (n-a-2) / (2*tnm3*npcm1);
    double F2 =  (n+a)*npam1 / (4*tnm1*tnm3*npcm2*npcm1);
    double F3 = -npam2*npam1*(n-a-2) / (8*tnm3*tnm3*tnm5*(n+c-3)*npcm2*npcm1);
    double E  = -npam1*(n-c-1) / (2*tnm3*npcm2*npcm1);

    double An = (1.0+F1*x)*Anm1 + (E + F2*x)*x*Anm2 + F3*x3*Anm3;
    double Bn = (1.0+F1*x)*Bnm1 + (E + F2*x)*x*Bnm2 + F3*x3*Bnm3;
    double r = An/Bn;

    *prec = fabs((F - r)/F);
    F = r;

    if(*prec < GSL_MACH_EPS || n > nmax) break;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An   /= RECUR_BIG;
      Bn   /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
      Anm3 /= RECUR_BIG;
      Bnm3 /= RECUR_BIG;
    }
    else if(fabs(An) < 1.0/RECUR_BIG || fabs(Bn) < 1.0/RECUR_BIG) {
      An   *= RECUR_BIG;
      Bn   *= RECUR_BIG;
      Anm1 *= RECUR_BIG;
      Bnm1 *= RECUR_BIG;
      Anm2 *= RECUR_BIG;
      Bnm2 *= RECUR_BIG;
      Anm3 *= RECUR_BIG;
      Bnm3 *= RECUR_BIG;
    }

    n++;
    Bnm3 = Bnm2;
    Bnm2 = Bnm1;
    Bnm1 = Bn;
    Anm3 = Anm2;
    Anm2 = Anm1;
    Anm1 = An;
  }

  *result = F;

  if(*prec > 10.0 * locEPS)
    return GSL_ELOSS;
  else
    return GSL_SUCCESS;
}


/* Asymptotic result for 1F1(a, b, x)  x -> -Infinity.
 * Assumes b-a != neg integer and b != neg integer.
 */
static
int
hyperg_1F1_asymp_negx(const double a, const double b, const double x,
                      double * result, double * prec
                      )
{
  double lg_b, sgn_b;
  double lg_bma, sgn_bma;
  double ln_pre;
  double ln_F;
  double prec_F;
  double F;
  int stat_b   = gsl_sf_lngamma_sgn_impl(b, &lg_b, &sgn_b);
  int stat_bma = gsl_sf_lngamma_sgn_impl(b-a, &lg_bma, &sgn_bma);
  
  if(stat_b == GSL_SUCCESS && stat_bma == GSL_SUCCESS) {
    gsl_sf_hyperg_2F0_series_impl(a, 1.0+a-b, -1.0/x, -1, &F, &prec_F);

    ln_pre = lg_b - a*log(-x) - lg_bma;
    ln_F = log(fabs(F));
  
    if(ln_pre + ln_F  <  GSL_LOG_DBL_MAX-1.0) {
      *result = sgn_b * sgn_bma * exp(ln_pre) * F;
      return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      return GSL_EOVRFLW;
    }
  }
  else {
    *result = 0.0;
    return GSL_EDOM;
  }
}


/* Asymptotic result for 1F1(a, b, x)  x -> +Infinity
 * Assumes b != neg integer and a != neg integer
 */
static
int
hyperg_1F1_asymp_posx(const double a, const double b, const double x,
                      double * result, double * prec
                      )
{
  double lg_b, sgn_b;
  double lg_a, sgn_a;
  double ln_pre;
  double ln_F;
  double prec_F;
  double F;

  int stat_b = gsl_sf_lngamma_sgn_impl(b, &lg_b, &sgn_b);
  int stat_a = gsl_sf_lngamma_sgn_impl(a, &lg_a, &sgn_a);

  if(stat_a == GSL_SUCCESS && stat_b == GSL_SUCCESS) {
    gsl_sf_hyperg_2F0_series_impl(b-a, 1.0-a, 1.0/x, -1, &F, &prec_F);

    ln_pre = lg_b - lg_a + x + (a-b)*log(x);

    if(ln_pre + ln_F  <  GSL_LOG_DBL_MAX) {
      *result = sgn_b * sgn_a * exp(ln_pre) * F;
      return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      return GSL_EOVRFLW;
    }
  }
  else {
    *result = 0.0;
    return GSL_EDOM;
  }
}


/* Use this when a is small.
 * We have to handle a few cases for a=0,a=-1,... to be correct.
 */
static
int
hyperg_1F1_small_a(const double a, const double b, const double x, double * result)
{
  double bma = b-a;
  double oma = 1.0-a;
  double ap1mb = 1.0+a-b;
  double abs_a = fabs(a);
  double abs_b = fabs(b);
  double abs_x = fabs(x);
  double abs_bma = fabs(bma);
  double abs_oma = fabs(oma);
  double abs_ap1mb = fabs(ap1mb);

  if(abs_a < locEPS) {
    /* a == 0 */
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(fabs(a+1.0) < locEPS) {
    /* a == -1 */
    *result = 1.0 + a/b * x;
    return GSL_SUCCESS;
  }
  else if(fabs(a+2.0) < locEPS) {
    /* a == -2 */
    *result = 1.0 + a/b*x * (1.0 + 0.5*(a+1.0)/(b+1.0)*x);
    return GSL_SUCCESS;
  }
  else if(fabs(x) < 8.0 || (b > 0.0 && abs_x < 0.7 * abs_b)) {
    /* Series is easy or is dominated and safe, though
     * it may be little slow to converge in the latter
     * case, being like Sum[(x/b)^n] in the worst case.
     */
    double prec;
    return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
  }
  else if(x < 0.0 && locMAX(abs_a,1.0)*locMAX(abs_ap1mb,1.0) < 0.99*abs_x) {
    double prec;
    return hyperg_1F1_asymp_negx(a, b, x, result, &prec);
  }
  else if(x > 0.0 && locMAX(abs_bma,1.0)*locMAX(abs_oma,1.0) < 0.99*abs_x) {
    double prec;
    return hyperg_1F1_asymp_posx(a, b, x, result, &prec);
  }
  else if(   (abs_b > 1000.0 && abs_x < 0.8 * abs_b)
          || (abs_b >  500.0 && abs_x < 0.5 * abs_b)
    ) {
    return gsl_sf_hyperg_1F1_large_b_impl(a, b, x, result);
  }
  else {
    /* We are left with small wedges around b=|x|, and
     * somewhat larger wedges around b=-|x|, as well as
     * a chunk for -500<b<0.
     */
    if(x < 0.0) {
      double prec;
      return hyperg_1F1_luke(a, b, x, result, &prec);
    }
    else {
    }
  }
}


/* Evaluate for small a:
 *   Y(a,0,b,x)  when a > 0
 *   Z(a,0,b,x)  when a < 0
 * Not defined for a = 0.
 */
static
int
hyperg_1F1_YZ_small_a(const double a, const double b, const double x, double * result)
{
  double F;
  hyperg_1F1_small_a(a, b, x, &F);

  if(F == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }

  if(a == 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(a > 0.0) {
    /* Y(a,0,b,x) := Gamma[a]/Gamma[1+a-b] 1F1(a, b, x)
     */
    double lg_ab, lg_a, sg_ab, sg_a;
    int stat_a  = gsl_sf_lngamma_sgn_impl(a,       &lg_a,  &sg_a);
    int stat_ab = gsl_sf_lngamma_sgn_impl(1.0+a-b, &lg_ab, &sg_ab);
    if(stat_a == GSL_SUCCESS && stat_ab == GSL_SUCCESS) {
      double sgF = (F > 0.0 ? 1.0 : -1.0);
      double lnF = log(fabs(F));
      double lnr = lnF + lg_a - lg_ab;
      if(lnr < GSL_LOG_DBL_MAX) {
        *result = sg_a * sg_ab * sgF * exp(lnr);
	return GSL_SUCCESS;
      }
      else {
        *result = 0.0;
	return GSL_EOVRFLW;
      }
    }
    else {
      *result = 0.0;
      return GSL_EDOM;
    }
  }
  else {
    /* Z(a,0,b,x) := 1/Gamma(1+a-b) 1F1(a,b,x)
     */
    double lg_ab, sg_ab;
    int stat_ab = gsl_sf_lngamma_sgn_impl(1.0+a-b, &lg_ab, &sg_ab);
    if(stat_ab == GSL_SUCCESS) {
      double sgF = (F > 0.0 ? 1.0 : -1.0);
      double lnF = log(fabs(F));
      double lnr = lnF - lg_ab;
      if(lnr < GSL_LOG_DBL_MAX) {
        *result = sg_ab * sgF * exp(lnr);
	return GSL_SUCCESS;
      }
      else {
        *result = 0.0;
	return GSL_EOVRFLW;
      }
    }
  }
}




/* Evaluates 1F1(a+N,b,x) and 1F1(a+N+1,b,x) by recursing
 * up from 1F1(a,b,x) and 1F1(a+1,b,x).
 *
 * Assumes a is small and positive.
 */
static
int
hyperg_1F1_recurse_posa(double a, double b, double x, int N, double * FN, double * FNp1)
{
  double F0, F1;
  int stat_0 = hyperg_1F1_small_a(a,	 b, x, &F0);
  int stat_1 = hyperg_1F1_small_a(a+1.0, b, x, &F1);

  if(stat_0 == GSL_EOVRFLW || stat_1 == GSL_EOVRFLW) {
    *FN   = 0.0;
    *FNp1 = 0.0;
    return GSL_EOVRFLW;
  }
  else {
    double lg_a0;   /* log(Gamma(a))     */
    double lg_ab0;  /* log(Gamma(1+a-b)) */
    double sg_ab0;
    int stat_lg_a0  = gsl_sf_lngamma_impl(a, &lg_a0);
    int stat_lg_ab0 = gsl_sf_lngamma_sgn_impl(1+a-b, &lg_ab0, &sg_ab0);

    double ln_Y0 = lg_a0 - lg_ab0 + log(fabs(F0));
    double ln_Y1 = lg_a0 - lg_ab0 + log(fabs(a/(1.0+a-b)*F1));

    if(ln_Y0 > GSL_LOG_DBL_MAX || ln_Y1 > GSL_LOG_DBL_MAX) {
      *FN   = 0.0;
      *FNp1 = 0.0;
      return GSL_EOVRFLW;
    }
    else {
      double lg_aN;       /* log(Gamma(a+N))     */
      double lg_abN;      /* log(Gamma(1+a+N-b)) */
      double sg_abN;
      int stat_lg_aN  = gsl_sf_lngamma_impl(a+N, &lg_aN);
      int stat_lg_abN = gsl_sf_lngamma_sgn_impl(1+a+N-b, &lg_abN, &sg_abN);

      double e0 = sg_ab0 * exp(lg_a0 - lg_ab0);
      double Y0 = e0 * F0;
      double Y1 = a/(1.0+a-b) * e0 * F1;
      double ln_FN, ln_FNp1;
      double YN, YNp1;

      hyperg_1F1_Y_recurse_up(a, b, x, 1, Y0, Y1, N+1, &YN, &YNp1);

      ln_FN   = -(lg_aN - lg_abN) + log(fabs(YN));
      ln_FNp1 = -(lg_aN - lg_abN) + log(fabs((1.0+a+N-b)/(a+N)*YNp1));
      
      if(ln_FN > GSL_LOG_DBL_MAX || ln_FNp1 > GSL_LOG_DBL_MAX) {
        *FN   = 0.0;
	*FNp1 = 0.0;
	return GSL_EOVRFLW;
      }
      else {
        double eN = exp(lg_aN - lg_abN);
        *FN   = YN / eN;
        *FNp1 = (1.0+a+N-b)/(a+N) / eN * YNp1;
        return GSL_SUCCESS;
      }
    }
  }
}




/* Evaluates 1F1(a-N,b,x) and 1F1(a-N-1,b,x) by recursing
 * down from 1F1(a,b,x) and 1F1(a-1,b,x).
 *
 */
static
int
hyperg_1F1_recurse_nega(double a, double b, double x, int N, double * FN, double * FNp1)
{
  double prec;
  double F0, F1;
  int stat_0 = hyperg_1F1_small_a(a,	 b, x, &F0);
  int stat_1 = hyperg_1F1_small_a(a-1.0, b, x, &F1);

  if(stat_0 == GSL_EOVRFLW || stat_1 == GSL_EOVRFLW) {
    *FN   = 0.0;
    *FNp1 = 0.0;
    return GSL_EOVRFLW;
  }
  else {
    double lg_ab0;  /* log(Gamma(1+a-b)) */
    double sg_ab0;
    int stat_lg_ab0 = gsl_sf_lngamma_sgn_impl(1+a-b, &lg_ab0, &sg_ab0);

    double ln_Y0 = -lg_ab0 + log(fabs(F0));
    double ln_Y1 = -lg_ab0 + log(fabs((a-b)*F1));

    if(ln_Y0 > GSL_LOG_DBL_MAX || ln_Y1 > GSL_LOG_DBL_MAX) {
      *FN   = 0.0;
      *FNp1 = 0.0;
      return GSL_EOVRFLW;
    }
    else {
      double lg_abN;      /* log(Gamma(1+a-N-b)) */
      double sg_abN;
      int stat_lg_abN = gsl_sf_lngamma_sgn_impl(1+a-N-b, &lg_abN, &sg_abN);

      double e0 = sg_ab0 * exp(-lg_ab0);
      double Y0 = e0 * F0;
      double Y1 = (a-b) * e0 * F1;
      double ln_FN, ln_FNp1;
      double YN, YNp1;

      hyperg_1F1_Z_recurse(a, b, x, 1, Y0, Y1, N+1, &YN, &YNp1);

      ln_FN   = lg_abN + log(fabs(YN));
      ln_FNp1 = lg_abN + log(fabs(YNp1/(a-N-b)));
      
      if(ln_FN > GSL_LOG_DBL_MAX || ln_FNp1 > GSL_LOG_DBL_MAX) {
        *FN   = 0.0;
	*FNp1 = 0.0;
	return GSL_EOVRFLW;
      }
      else {
        double eN = sg_abN * exp(lg_abN);
        *FN   = eN * YN;
        *FNp1 = eN * YNp1 / (a-N-b);
        return GSL_SUCCESS;
      }
    }
  }
}


/* 1F1(a,2a,x) = Gamma(a + 1/2) E(x) (|x|/4)^(-a+1/2) scaled_I(a-1/2,|x|/2)
 *
 * E(x) = exp(x) x > 0
 *      = 1      x < 0
 *
 * a >= 1/2
 */
static
int
hyperg_1F1_beq2a_pos(const double a, const double x, double * result)
{
  if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else {
    double I;
    int stat_I = gsl_sf_bessel_Inu_scaled_impl(a-0.5, 0.5*fabs(x), &I);
    if(stat_I == GSL_SUCCESS) {
      double lg;
      double lr;
      gsl_sf_lngamma_impl(a + 0.5, &lg);
      lr = lg + locMAX(x,0.0) + (0.5-a)*log(0.25*fabs(x)) + log(fabs(I));
      return gsl_sf_exp_impl(lr, result);
    }
    else {
      *result = 0.0;
      return stat_I;
    }
  }
}


/* Determine middle parts of diagonal recursion along b=2a
 * from two endpoints, i.e.
 *
 * given:  M(a,b)      and  M(a+1,b+2)
 * get:    M(a+1,b+1)  and  M(a,b+1)
 */
#ifdef HAVE_INLINE
inline
#endif
static
int
hyperg_1F1_diag_step(const double a, const double b, const double x,
                     const double Mab, const double Map1bp2,
                     double * Map1bp1, double * Mabp1)
{
  if(a == b) {
    *Map1bp1 = Mab;
    *Mabp1   = Mab - x/(b+1.0) * Map1bp2;
  }
  else {
    *Map1bp1 = Mab - x * (a-b)/(b*(b+1.0)) * Map1bp2;
    *Mabp1   = (a * *Map1bp1 - b * Mab)/(a-b);
  }
  return GSL_SUCCESS;
}


/* Determine endpoint of diagonal recursion.
 *
 * given:  M(a,b)    and  M(a+1,b+2)
 * get:    M(a+1,b)  and  M(a+1,b+1)
 */
#ifdef HAVE_INLINE
inline
#endif
static
int
hyperg_1F1_diag_end_step(const double a, const double b, const double x,
                         const double Mab, const double Map1bp2,
                         double * Map1b, double * Map1bp1)
{
  *Map1bp1 = Mab - x * (a-b)/(b*(b+1.0)) * Map1bp2;
  *Map1b   = Mab + x/b * *Map1bp1;
  return GSL_SUCCESS;
}


/* Handle the case of a and b both positive integers.
 * Assumes a > 0 and b > 0.
 */
static
int
hyperg_1F1_ab_posint(const int a, const int b, const double x, double * result)
{
  double ax = fabs(x);

  if(a == b) {
    return gsl_sf_exp_impl(x, result);             /* 1F1(a,a,x) */
  }
  else if(a == 1) {
    return gsl_sf_exprel_n_impl(b-1, x, result);   /* 1F1(1,b,x) */
  }
  else if(b == a + 1) {
    double K;
    int stat_K = gsl_sf_exprel_n_impl(a, -x, &K);  /* 1F1(1,1+a,-x) */
    if(K == 0.0) {
      *result = 0.0;
      return stat_K;
    }
    if(stat_K == GSL_SUCCESS) {
      double lK = log(fabs(K));
      double lr = lK + x;
      return gsl_sf_exp_sgn_impl(lr, K, result);
    }
    else {
      *result = 0.0;
      return stat_K;
    }
  }
  else if(a == b + 1) {
    *result = exp(x) * (1.0 + x/b);
    return GSL_SUCCESS;
  }
  else if(a == b + 2) {
    *result = exp(x) * (1.0 + x/b*(2.0 + x/(b+1)));
    return GSL_SUCCESS;
  }
  else if(b == 2*a) {
    return hyperg_1F1_beq2a_pos(a, x, result);  /* 1F1(a,2a,x) */
  }
  else if(   ( b < 10 && a < 10 && ax < 5.0 )
          || ( b > a*ax )
	  || ( b > a && ax < 5.0 )
    ) {
    double prec;
    return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
  }
  else if(a > b) {
    /* Forward recursion from a=b.
     * Note that a > b + 1 as well, since we already tried a = b + 1.
     */
    if(x + log(fabs(x/b)) < GSL_LOG_DBL_MAX-2.0) {
      double ex = exp(x);
      int n;
      double Mnm1 = ex; 		/* 1F1(b,b,x)   */
      double Mn   = ex * (1.0 + x/b);   /* 1F1(b+1,b,x) */
      double Mnp1;
      for(n=b+1; n<a; n++) {
        Mnp1 = ((b-n)*Mnm1 + (2*n-b+x)*Mn)/n;
	Mnm1 = Mn;
	Mn   = Mnp1;
      }
      *result = Mn;
      return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      return GSL_EOVRFLW;
    }
  }
  else if(2*a > b){
    double Mnm1;
    double Mn;
    double Mnp1;
    int a_start;
    int s = 0;

    if(GSL_IS_EVEN(b)) {
      /* Forward recursion from b/2+1 and b/2.
       */
      a_start = b/2 + 1;

      if(b == 2) {
        s += gsl_sf_exprel_impl(x, &Mnm1);  /* 1F1(1,2,x) =(e^x-1)/x */
        s += gsl_sf_exp_impl(x, &Mn);       /* 1F1(2,2,x) = e^x      */
      }
      else {
        double M12, M11;
        s += hyperg_1F1_beq2a_pos(b/2,   x, &Mnm1);  /* 1F1(b/2,b,x)     */
        s += hyperg_1F1_beq2a_pos(b/2+1, x, &M12);   /* 1F1(b/2+1,b+2,x) */
        s += hyperg_1F1_diag_end_step(b/2, b, x, Mnm1, M12, &Mn, &M11);
      }
    }
    else {
      /* Forward recursion from (b+1)/2 and (b-1)/2.
       */
      a_start = (b+1)/2;

      if(b == 1) {
        Mnm1 = 0.0;                     /* 1F1(0,1,x) */
        s += gsl_sf_exp_impl(x, &Mn);   /* 1F1(1,1,x) */
      }
      else {
        double M00, M12;
        s += hyperg_1F1_beq2a_pos((b-1)/2, x, &M00);
        s += hyperg_1F1_beq2a_pos((b+1)/2, x, &M12);
        s += hyperg_1F1_diag_step((b-1)/2, b-1, x, M00, M12, &Mnm1, &Mn);
      }
    }

    if(s == 0) {
      int n;
      for(n=a_start; n<a; n++) {
        Mnp1 = ((b-n)*Mnm1 + (2*n-b+x)*Mn)/n;
	Mnm1 = Mn;
	Mn   = Mnp1;
      }
      *result = Mn;
      return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      return GSL_EFAILED;
    }
  }
  else if(b > 2*a + x) {
    /* b > x, so use Gautschi series representation of
     * continued fraction. Then recurse backward since
     * we are in the stable region for that as well.
     */
    double rap;
    int stat_CF1 = hyperg_1F1_CF1_p_ser(a, b, x, &rap);
    double ra = 1.0 + x/a * rap;

    if(stat_CF1 == GSL_SUCCESS || stat_CF1 == GSL_EMAXITER) {
      double Ma   = GSL_SQRT_DBL_MIN;
      double Map1 = ra * Ma;
      double Mnp1 = Map1;
      double Mn   = Ma;
      double Mnm1;
      int n;
      for(n=a; n>0; n--) {
        Mnm1 = (n * Mnp1 - (2*n-b+x) * Mn) / (b-n);
        Mnp1 = Mn;
        Mn   = Mnm1;
      }
      *result = Ma/Mn;
      return stat_CF1;
    }
    else {
      *result = 0.0;
      return stat_CF1;
    }
  }
  else {
    /* 2a + x > b > 2a
     */
    double ra;
    int stat_CF1;
    int n;
    double Ma;
    double Mnm1;
    double Mn;
    double Mnp1;

    if(b > x) {
      /* Gautschi stable region for continued fraction.
       */
      double rap;
      stat_CF1 = hyperg_1F1_CF1_p_ser(a, b, x, &rap);
      ra = 1.0 + x/a * rap;
    }
    else {
      /* This is the "anomalous convergence" region.
       * Direct application of any method related
       * to the continued fraction will not work.
       * However, we can make use of the relation
       *
       * M(a+1,b,x)/M(a,b,x) = M(b-a-1,b,-x)/M(b-a,b,-x),
       *
       * and the second ratio is in the Gautschi stable
       * region and can be evaluated. [The Kummer transform
       * is a reflection about b=2a].
       */
      double rap_Kummer;
      double ra_Kummer;
      stat_CF1 = hyperg_1F1_CF1_p_ser(b-a-1, b, -x, &rap_Kummer);
      ra_Kummer = 1.0 + (-x/(b-a-1)) * rap_Kummer;
      ra = 1.0/ra_Kummer;
    }

    /* Recurse forward to a=b to determine normalization.
     * Since b < 2a + x, this is stable.
     */
    Ma   = GSL_SQRT_DBL_MIN;
    Mnm1 = Ma;
    Mn   = ra * Mnm1;

    for(n=a+1; n<b; n++) {
      Mnp1 = ((b-n)*Mnm1 + (2*n-b+x)*Mn)/n;
      Mnm1 = Mn;
      Mn   = Mnp1;
    }
    *result = Ma/Mn * exp(x);
    return GSL_SUCCESS;
  }
}


static
int
hyperg_1F1_ab_negint(const double a, const double b, const double x, double * result)
{
  /* Safe to recurse backward with Z function.
   * Z(ap,N,b,x) := 1/Gamma(1+ap-N-b) 1F1(ap-N,b,x)
   * Z[n+1] + (b-2a-x+2n) Z[n] + (a-n)(1+a-n-b) Z[n-1] = 0
   */
  double Znm1 = 1.0;	     /* Z(0,0,b,x) */
  double Zn   = 1.0 - x/b;   /* Z(0,1,b,x) */
  double Znp1;
  int n;
  for(n=2; n<a; n++) {
    Znp1 = -(b-x-2*n) * Zn + n*(1-n-b) * Znm1;
    Znm1 = Zn;
    Zn   = Znp1;
  }
  if(Zn == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    double lg_ab;
    gsl_sf_lngamma_impl(1.0+a-b, &lg_ab);  /* 1+a-b > 0  here */
    return gsl_sf_exp_sgn_impl(lg_ab + log(fabs(Zn)), Zn, result);
  }
}


/* Handle a = positive integer cases.
 * Because of the way we use the Z and Y
 * functions, this is actually useful.
 * It is also somewhat clearer to separate
 * this case explicitly.
 */
static
int
hyperg_1F1_a_posint(const int a, const double b, const double x);


/* Handle a = negative integer cases.
 * Assumes b is such that the numerator terminates
 * before the denominator, which is of course required
 * in any case. Also assumes a <= 0.
 */
static
int
hyperg_1F1_a_negint(const int a, const double b, const double x, double * result)
{
  if(a == 0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(a == -1) {
    *result = 1.0 + a/b * x;
    return GSL_SUCCESS;
  }
  else if(a == -2) {
    *result = 1.0 + a/b*x * (1.0 + 0.5*(a+1.0)/(b+1.0)*x);
    return GSL_SUCCESS;
  }
  else {
    int N = -a;
    double poly = 1.0;
    int k;
    for(k=N-1; k>=0; k--) {
      double t = (a+k)/(b+k) * (x/(k+1));
      poly = 1.0 + poly * t;
    }
    *result = poly;
    return GSL_SUCCESS;
  }
}


/* Evaluate Y(a, 0, b, x), Y(a, 1, b, x) for b > 2a, a > 0.
 * Since we cannot recurse forward on a in this regime,
 * we use Steed's continued fraction, then recurse backward
 * and evaluate at small a to obtain the normalization.
 */
static
int
hyperg_1F1_Y_agt0_bgt2a_recurse(const double a, const double b, const double x,
                                double * Ya, double * Yap1)
{
  int N     = floor(a);
  double ap = a - N;
  double yapp1, yap;
  double yNp1, yN;
  double fN;
  double yap_true;
  double scale;
  double lg_a, sg_a, lg_ab, sg_ab;
  double pre;
  int stat_CF1;
  int stat_rec;

  /* We cannot allow the lower a value to be zero since
   * we are working in terms of Y.
   */
  if(ap == 0.0) {
    ap += 1.0;
    N  -= 1;
  }

  stat_CF1 = hyperg_1F1_Y_CF1(ap, b, N, x, &fN);  /* fN = Y_{N+1}/Y_N */

  yN   = GSL_SQRT_DBL_MIN;
  yNp1 = fN * yN;

  stat_rec = hyperg_1F1_Y_recurse_down(ap, b, x, N, yNp1, yN, 0, &yapp1, &yap);

  hyperg_1F1_YZ_small_a(ap, b, x, &yap_true);
  scale = yap_true / yap;

  *Ya   = scale * yN;
  *Yap1 = scale * yNp1;

  return GSL_SUCCESS;
}



/* Evaluate 1F1(a, b, x) for b < 2a, a > 0.
 * We evaluate near b=2a and then recurse forward.
 */
static
int
hyperg_1F1_Y_agt0_blt2a_recurse(const double a, const double b, const double x,
                                double * Ya, double * Yap1)
{
  /* a = N + ap
   * a_start = M + ap
   * M <= N
   */
  int N = floor(a);
  int M = ( b > 0.0 ? floor(0.5*b) : 0 );
  double ap = a - N;
  double aM = ap + M + 1.0;
  double YM, YMp1;
  double YN, YNm1;

  int stat_1 = hyperg_1F1_Y_agt0_bgt2a_recurse(aM, b, x, &YM, &YMp1);

  int stat_2 = hyperg_1F1_Y_recurse_up(aM, b, x, 0, YM, YMp1, N-M, &YNm1, &YN);

  *Ya   = YNm1;
  *Yap1 = YN;
  if(stat_1 != GSL_SUCCESS)
    return stat_1;
  else if(stat_2 != GSL_SUCCESS)
    return stat_2;
  else
    return GSL_SUCCESS;
}


int
gsl_sf_hyperg_1F1_int_impl(const int a, const int b, const double x, double * result)
{
  if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(a == b) {
    return gsl_sf_exp_impl(x, result);
  }
  else if(b == 0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(a == 0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(b < 0 && (a < b || a > 0)) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(a == -1) {
    *result = 1.0 + (double)a/b * x;
    return GSL_SUCCESS;
  }
  else if(a == -2) {
    *result = 1.0 + (double)a/b*x * (1.0 + 0.5*(a+1.0)/(b+1.0)*x);
    return GSL_SUCCESS;
  }
  else if(a < 0 && b < 0) {
    return hyperg_1F1_ab_negint(a, b, x, result);
  }
  else if(a < 0 && b > 0) {
    /* Use Kummer to reduce it to the positive integer case.
     * Note that b > a, strictly, since we already trapped b = a.
     */
    double Kummer_1F1;
    int stat_K = hyperg_1F1_ab_posint(b-a, b, -x, &Kummer_1F1);
    if(stat_K == GSL_SUCCESS) {
      if(Kummer_1F1 == 0.0) {
        *result = 0.0;
	return GSL_SUCCESS;
      }
      else {
        double lnr = log(fabs(Kummer_1F1)) + x;
        return gsl_sf_exp_sgn_impl(lnr, Kummer_1F1, result); 
      }
    }
    else {
      *result = 0.0;
      return stat_K;
    }
  }
  else {
    /* a > 0 and b > 0 */
    return hyperg_1F1_ab_posint(a, b, x, result);
  }
}


int
gsl_sf_hyperg_1F1_impl(const double a, const double b, const double x,
                       double * result
                       )
{
  const double bma = b - a;
  const double amb = a - b; 
 
 /* Note that we draw a distinction between the negative
  * integer cases and the zero cases. This makes sense since
  * it may be possible to specify, for instance, b as a very
  * small number and mean it, whereas the ability to specify
  * a number close to a given nonzero negative integer is
  * limited by the precision of the FP representation.
  * For the negative integer test, we resort to an
  * (arbitrary) nearness measure.
  */
  const int a_integer = ( fabs(a - rint(a)) < locEPS );
  const int b_integer = ( fabs(b - rint(b)) < locEPS );
  const int a_neg_integer = ( a < -0.1 &&  a_integer );
  const int b_neg_integer = ( b < -0.1 &&  b_integer );
  const int bma_neg_integer = ( bma < -0.1  &&  fabs(bma - rint(bma)) < locEPS );
  const int amb_neg_integer = ( amb < -0.1  &&  fabs(amb - rint(amb)) < locEPS );

  /* Testing for this before testing a and b
   * is somewhat arbitrary. The result is that
   * we can have 1F1(0,0,0) = 1. Whatever.
   */
  if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }

  if(b == 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }

  if(a == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }

  /* case: a==b; exp(x)
   * It's good to test exact equality now.
   * We also test approximate equality later.
   * Note that we need to test this before we
   * test for the domain error below (denominator
   * zero before numerator), because a and b
   * might be negative integers with a==b, which would
   * wrongly trip the test below.
   */
  if(a == b) {
    return gsl_sf_exp_impl(x, result);
  }

  /* case: denominator zeroes before numerator
   * Note that we draw a distinction between the
   * negative integers and zero. Approximate zero
   * is handled in the cases below, and is not a
   * negative integer. So we can test the negative
   * integer case now since there is no interference
   * between the cases.
   */
  if(b_neg_integer && !(a_neg_integer && a > b + 0.1)) {
    *result = 0.0;
    return GSL_EDOM;
  }

  /* case: a and b very small; 1 + a/b (exp(x)-1)
   * Note that neither a nor b is zero, since
   * we eliminated that with the above tests.
   */
  if(fabs(a) < 10.0*locEPS && fabs(b) < 10.0*locEPS) {
    double expm1;
    int stat_e = gsl_sf_expm1_impl(x, &expm1);
    if(stat_e == GSL_SUCCESS) {
      double sa = ( a > 0.0 ? 1.0 : -1.0 );
      double sb = ( b > 0.0 ? 1.0 : -1.0 );
      double se = ( expm1 > 0.0 ? 1.0 : -1.0 );
      double lnr = log(fabs(a)) - log(fabs(b)) + log(fabs(expm1));
      double hx;
      int stat_hx = gsl_sf_exp_sgn_impl(lnr, sa * sb * se, &hx);
      if(stat_hx == GSL_SUCCESS) {
        *result = 1.0 + hx;
	return GSL_SUCCESS;
      }
      else {
        *result = 0.0;
        return stat_hx;
      }
    }
    else {
      *result = 0.0;
      return stat_e;
    }
  }

  /* case: b very small
   * (and a not very small, since we did not trigger the above test)
   */
  if(fabs(b) < 10.0*locEPS) {
    double F_renorm;
    int stat_F = hyperg_1F1_renorm_b0(a, x, &F_renorm);
    if(F_renorm == 0.0) {
      /* It is possible to get zero, since we might
       * hit a zero of the Bessel function.
       */
      *result = 0.0;
      return stat_F;
    }
    else {
      double sF = ( F_renorm > 0.0 ? 1.0 : -1.0 );
      double sb = ( b > 0.0 ? 1.0 : -1.0 );
      double lnr = log(fabs(F_renorm)) - log(fabs(b));
      return gsl_sf_exp_sgn_impl(lnr, sF * sb, result);
    }
  }

  /* Now we make the arbitrary "near an integer" test.
   * Note that we deal not only with the negative
   * integers, but all the integers. Basically, the
   * integer cases are always a problem for the
   * recursions in terms of the Y() or Z() functions.
   */
  if(a_integer && b_integer) {
    int inta = floor(a + 0.1);
    int intb = floor(b + 0.1);
    return gsl_sf_hyperg_1F1_int_impl(inta, intb, x, result);
  }

  /* case: approximate a==b;
   *
   * 1F1(a,a+eps,x) = exp(ax/b) (1 + eps x^2 (v2 + v3 x + ...) + ...)
   *
   *   v2 = a/(2b^2(b+1))
   *   v3 = a(b-2a)/(3b^3(b+1)(b+2))
   *   ...
   *
   * See [Luke, Mathematical Functions and Their Approximations, p.292]
   *
   * This cannot be used for b near a negative integer or zero.
   * Also, if x/b is large the deviation from exp(x) behaviour grows.
   */
  if(b > locEPS && fabs(bma) < GSL_SQRT_MACH_EPS && fabs(b) > fabs(x) ) {
    double eps = bma;
    double exab;
    int stat_e = gsl_sf_exp_impl(a*x/b, &exab);
    if(stat_e == GSL_SUCCESS) {
      double v2 = a/(2.0*b*b*(b+1.0));
      double v3 = a*(b-2.0*a)/(3.0*b*b*b*(b+1.0)*(b+2.0));
      double v  = v2 + v3 * x;
      *result = exab * (1.0 + eps*x*x*(v2 + v3 * x));
      return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      return stat_e;
    }
  }

  /* case: a = neg integer; use explicit evaluation */
  if(a_neg_integer) {
    int inta = floor(a + 0.1);
    return hyperg_1F1_a_negint(inta, b, x, result);
  }

  /* b-a = negative integer; Kummer and explicit evaluation */
  if(bma_neg_integer) {
    int    intbma = floor(bma + 0.1);
    double Kummer_1F1;
    int    stat_K = hyperg_1F1_a_negint(intbma, b, -x, &Kummer_1F1);
    if(Kummer_1F1 == 0.0) {
      *result = 0.0;
      return stat_K;
    }
    else {
      double lnr = log(fabs(Kummer_1F1)) + x;
      double sK  = ( Kummer_1F1 > 0.0 ? 1.0 : - 1.0 );
      if(lnr < GSL_LOG_DBL_MAX) {
        *result = sK * exp(lnr);
        return stat_K;
      }
      else {
        *result = 0.0;  /* FIXME: should be Inf */
        return GSL_EOVRFLW;
      }
    }
  }

  /* Trap the generic cases where some form
   * of series evaluation will work.
   */
  if(fabs(x) < 10.0) {

    if( (fabs(a) < 20.0 && fabs(b) < 20.0) || (b >= fabs(a)) ) {
      /* Arguments small enough to evaluate series directly
       * or series is dominated and safe.
       */
      double prec;
      return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
    }

    if( (fabs(bma) < 20.0 && fabs(b) < 20.0) || (b >= fabs(bma)) ) {
      /* Use Kummer transformation to render series safe.
       * We do not have to worry about overflow in
       * exp(x) * Kummer_1F1, because neither term can be very large.
       */
      double prec;
      double Kummer_1F1;
      double Ex = exp(x);
      int stat_K = gsl_sf_hyperg_1F1_series_impl(bma, b, -x, &Kummer_1F1, &prec);
      *result = Ex * Kummer_1F1;
      return stat_K;
    }
  }

  /* When b is positive, the series can be dominated,
   * as long as x is not relatively too large.
   */
  if(b > locMAX(fabs(a),1.0)*fabs(x)) {
    double prec;
    return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
  }

  /* Kummer transformed version of above.
   */
  if(b > locMAX(fabs(bma),1.0)*fabs(x)) {
    double prec;
    double Kummer_1F1;
    int stat_K = gsl_sf_hyperg_1F1_series_impl(bma, b, -x, &Kummer_1F1, &prec);
    if(Kummer_1F1 == 0.0) {
      *result = 0.0;
      return stat_K;
    }
    else {
      double lnr = log(fabs(Kummer_1F1)) + x;
      double sgK = (Kummer_1F1 > 0.0 ? 1.0 : -1.0);
      if(lnr < GSL_LOG_DBL_MAX) {
        *result = sgK * exp(lnr);
        return stat_K;
      }
      else {
        *result = 0.0; /* FIXME: should be Inf */
        return GSL_EOVRFLW;
      }
    }
  }

  /* Handle small a explicitly. This helps clarify
   * the thinking for the recursize cases below.
   * It is also a minor premature optimization.
   */
  if(fabs(a) < 2.0) {
    return hyperg_1F1_small_a(a, b, x, result);
  }


  if(a > 0.0 && b > 0.0) {
    int N     = floor(a);
    double ap = a - N;
    int n     = 2;
    double Ynm1, Yn;
    double Ya, Yap1;
    double Fa, Fap1;
    hyperg_1F1_YZ_small_a(ap + n - 1, b, x, &Ynm1);
    hyperg_1F1_YZ_small_a(ap + n,     b, x, &Yn);
    hyperg_1F1_Y_recurse_up(ap, b, x, n, Ynm1, Yn, N+1, &Ya, &Yap1);
    hyperg_1F1_from_Y(a, b, x, Ya, Yap1, &Fa, &Fap1);
    *result = Fa;
    return GSL_SUCCESS;
  }

  /* Handle the region where we cannot recurse forward.
   */
  if(a > 0.0 && b >= 2.0*a) {
    double Ya, Yap1;
    double Fap1;
    int stat_r = hyperg_1F1_Y_agt0_bgt2a_recurse(a, b, x, &Ya, &Yap1);
    return hyperg_1F1_from_Y(a, b, x, Ya, Yap1, result, &Fap1);
  }


  /* Handle the region involving stable forward recursion.
   */
  if(a > 0.0 && b < 2.0*a) {
    double Ya, Yap1;
    double Fap1;
    int stat_r = hyperg_1F1_Y_agt0_blt2a_recurse(a, b, x, &Ya, &Yap1);
    return hyperg_1F1_from_Y(a, b, x, Ya, Yap1, result, &Fap1);
  }




  /* Large negative x asymptotic.
   */
  if(   x < -10.0
     && locMAX(fabs(a),1.0)*locMAX(fabs(1.0+a-b),1.0) < 1.2*fabs(x)
     && !b_neg_integer
     && !bma_neg_integer
    ) {
    double prec;
    return hyperg_1F1_asymp_negx(a, b, x, result, &prec);
  }


  /* Large positive x asymptotic.
   */
  if(   x > 10.0
     && locMAX(fabs(bma),1.0)*locMAX(fabs(1.0-a),1.0) < 1.2*fabs(x)
     && !b_neg_integer
     && !a_neg_integer
    ) {
    double prec;
    return hyperg_1F1_asymp_posx(a, b, x, result, &prec);
  }


  /* Large positive a obtained by recursion.
   */
  if(a > 20.0 && a < INT_MAX-2 && b > 0.0 && x > 0.0) {
    double FN, FNp1;
    int N     = floor(a);
    double ap = a - N;
    int stat_r;
    if(ap == 0.0) {
      ap += 1.0;
      N  -= 1;
    }
    stat_r = hyperg_1F1_recurse_posa(ap, b, x, N, &FN, &FNp1);
    *result = FN;
    return stat_r;
  }


  /* Large negative a by recursion.
   */
  if(a < -20.0 && a > INT_MIN+1 && b > 0.0 && x > 0.0) {
    double FN, FNp1;
    int N     = -floor(a);
    double ap = a + N;
    int stat_r;
    if(ap == 0.0) {
      ap -= 1.0;
      N  -= 1;
    }
    stat_r = hyperg_1F1_recurse_nega(ap, b, x, N, &FN, &FNp1);
    *result = FN;
    return stat_r;
  }


  /* Luke in the canonical case.
   */
  if(x < 0.0 && !a_neg_integer && !bma_neg_integer) {
    double prec;
    return hyperg_1F1_luke(a, b, x, result, &prec);
  }


  /* Luke with Kummer transformation.
   */
  if(x > 0.0 && !a_neg_integer && !bma_neg_integer) {
    double prec;
    double Kummer_1F1;
    double ex;
    int stat_F = hyperg_1F1_luke(b-a, b, -x, &Kummer_1F1, &prec);
    int stat_e = gsl_sf_exp_impl(x, &ex);
    if(stat_F == GSL_SUCCESS && stat_e == GSL_SUCCESS) {
      double lnr = log(fabs(Kummer_1F1)) + x;
      if(lnr < GSL_LOG_DBL_MAX) {
        *result = ex * Kummer_1F1;
	return GSL_SUCCESS;
      }
      else {
        *result = 0.0;  /* FIXME: should be Inf */
	return GSL_EOVRFLW;
      }
    }
    else if(stat_F != GSL_SUCCESS) {
      *result = 0.0;
      return stat_F;
    }
    else {
      *result = 0.0;
      return stat_e;
    }
  }



  /* At this point we have no more tricks. Instead we must
   * proceed systematically. If a>0 we reduce it to 0<a<1
   * by backward recursion, then use the small-a evaluation
   * to normalize. If a<0, we can recurse backward directly.
   * See [Temme, Numer. Math. 41, 63 (1983)].
   *
   * Huh? I thought forward on 'a' was stable?? That's what
   * Temme says... p.65
   */
  *result = 0.0;
  return GSL_EUNIMPL;
}


int
gsl_sf_hyperg_1F1_int_e(const int m, const int n, double x, double * result)
{
  int status = gsl_sf_hyperg_1F1_int_impl(m, n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hyperg_1F1_int_e", status);
  }
  return status;
}


int
gsl_sf_hyperg_1F1_e(double a, double b, double x, double * result)
{
  int status = gsl_sf_hyperg_1F1_impl(a, b, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hyperg_1F1_e", status);
  }
  return status;
}


double
gsl_sf_hyperg_1F1_int(int m, int n, double x)
{
  double y;
  int status = gsl_sf_hyperg_1F1_int_impl(m, n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hyperg_1F1_int", status);
  }
  return y;
}


double
gsl_sf_hyperg_1F1(double a, double b, double x)
{
  double y;
  int status = gsl_sf_hyperg_1F1_impl(a, b, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hyperg_1F1", status);
  }
  return y;
}
