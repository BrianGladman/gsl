/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "hyperg.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_hyperg.h"

#define locMAX(a,b)     ((a) > (b) ? (a) : (b))
#define locEPS          (1000.0*GSL_MACH_EPS)


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


/* Assumes b-a != neg integer and b != neg integer.
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


/* Assumes b != neg integer and a != neg integer
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


/* Use this as it is used below. Basically, it assumes
 * that a is small and that either x or b is large since,
 * if they were all small, the series would have been
 * evaluated directly.
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
    *result = 1.0 + a/b * x + 0.5*a*(a+1.0)/(b*(b+1.0)) *x*x;
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
    double prec;
    return hyperg_1F1_luke(a, b, x, result, &prec);
  }
}


/* Do the (stable) upward recursion on the parameter 'a'.
 * Work in terms of the function
 *     Y(n) := Gamma(a+n)/Gamma(1+a+n-b) 1F1(a+n;b;x)
 *
 *     (a+n+1-b)Y(n+1) + (b-2a-2n-x)Y(n) + (a+n-1)Y(n-1) = 0
 *
 * Note that this will bomb if a-b+1 is a negative integer and
 * the recursion passes through that integer.
 */
static
int
hyperg_1F1_Y_recurse_posa(double a, double b, double x,
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


/* Evaluates 1F1(a+N,b,x) and 1F1(a+N+1,b,x) by recursing
 * up from 1F1(a,b,x) and 1F1(a+1,b,x).
 *
 * Assumes a is small and positive.
 */
static
int
hyperg_1F1_recurse_posa(double a, double b, double x, int N, double * FN, double * FNp1)
{
  double prec;
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

      hyperg_1F1_Y_recurse_posa(a, b, x, 1, Y0, Y1, N+1, &YN, &YNp1);

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


/* Recursion in the negative a direction.
 *
 * Y[n] := 1F1(a-n,b,x)/Gamma(1+a-n-b)
 *
 * Y[n+1] + (b-2a-x+2n) Y[n] + (a-n)(1+a-n-b) Y[n-1] = 0
 */
static
int
hyperg_1F1_Y_recurse_nega(const double a, const double b, const double x,
                          int n, double Ynm1, double Yn,
		          int N,
		          double * YNm1, double * YN)
{
  int k;
  double Ykm1 = Ynm1;
  double Yk   = Yn;
  double Ykp1;

  for(k=n; k<N; k++) {
    Ykp1 = -(b-2.0*a-x+2.0*k) * Yk - (a-k)*(1.0+a-k-b) * Ykm1;
    Ykm1 = Yk;
    Yk   = Ykp1;
  }
  
  *YNm1 = Ykm1;
  *YN   = Yk;
  return GSL_SUCCESS;
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

      hyperg_1F1_Y_recurse_nega(a, b, x, 1, Y0, Y1, N+1, &YN, &YNp1);

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


/* Handle a = neg integer cases.
 * Assumes b is such that the numerator terminates
 * before the denominator, which is of course required
 * in any case.
 *
 * FIXME: how does this behave for b < 0 and/or x < 0 ???
 * The question is whether or not the recursion is stable.
 */
static
int
hyperg_1F1_a_negint(int a, const double b, const double x, double * result)
{
  if(a == -1) {
    *result = 1.0 + a/b * x;
    return GSL_SUCCESS;
  }
  else if(a == -2) {
    *result = 1.0 + a/b * x + 0.5*a*(a+1.0)/(b*(b+1.0)) *x*x;
    return GSL_SUCCESS;
  }
  else {
    double FN, FNp1;
    int N      = -a;
    double ap  = 0.0;
    int stat_r = hyperg_1F1_recurse_nega(ap, b, x, N, &FN, &FNp1);
    *result = FN;
    return stat_r;
  }
}


int
gsl_sf_hyperg_1F1_impl(const double a, const double b, const double x,
                       double * result
                       )
{
  int a_neg_integer;    /*  a   negative integer  */
  int b_neg_integer;    /*  b   negative integer  */
  int bma_neg_integer;  /*  b-a negative integer  */
  int amb_neg_integer;  /*  a-b negative integer  */

  const double bma = b - a;
  const double amb = a - b;

/*
double prec;
return hyperg_1F1_luke(a, b, x, result, &prec);
*/
  a_neg_integer = ( a < 0.0  &&  fabs(a - rint(a)) < locEPS );
  b_neg_integer = ( b < 0.0  &&  fabs(b - rint(b)) < locEPS );
  bma_neg_integer = ( bma < 0.0  &&  fabs(bma - rint(bma)) < locEPS );
  amb_neg_integer = ( amb < 0.0  &&  fabs(amb - rint(amb)) < locEPS );
  
  if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }

  /* case: a==0 */
  if(fabs(a) < locEPS) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  
  /* case: a==b,  exp(x) */
  if(fabs(bma) < locEPS) {
    return gsl_sf_exp_impl(x, result);
  }

  /* case: denominator zeroes before numerator */
  if(b_neg_integer && !(a_neg_integer && a > b + 0.1)) {
    *result = 0.0;
    return GSL_EDOM;
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
      int stat_K = gsl_sf_hyperg_1F1_series_impl(bma, b, x, &Kummer_1F1, &prec);
      *result = Ex * Kummer_1F1;
      return stat_K;
    }
  }


  /* a = negative integer
   */
  if(a_neg_integer) {
    int inta = floor(a + 0.1);
    return hyperg_1F1_a_negint(inta, b, x, result);
  }


  /* b-a = negative integer
   */
  if(bma_neg_integer) {
    int    intbma = floor(bma + 0.1);
    double Kummer_1F1;
    int    stat_K = hyperg_1F1_a_negint(intbma, b, -x, &Kummer_1F1);
    double lnr    = log(fabs(Kummer_1F1)) + x;
    if(lnr < GSL_LOG_DBL_MAX) {
      *result = exp(x) * Kummer_1F1;
      return stat_K;
    }
    else {
      *result = 0.0;  /* FIXME: should be Inf */
      return GSL_EOVRFLW;
    }
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
gsl_sf_hyperg_1F1_e(double a, double b, double x, double * result)
{
  int status = gsl_sf_hyperg_1F1_impl(a, b, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hyperg_1F1_e", status);
  }
  return status;
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
