/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_psi.h"
#include "gsl_sf_hyperg.h"

#define locMAX(a, b)    ((a) > (b) ? (a) : (b))
#define locEPS          (1000.0*GSL_MACH_EPS)


static int hyperg_2F1_series(const double a, const double b, const double c,
                             const double x, 
			     double * result, double * prec
			     )
{
  double sum = 1.0;
  double del = 1.0;
  double delmax = 0.0;
  double k = 0.0;
  int i = 0;

  if(fabs(c) < locEPS) {
    *prec   = 1.0;
    *result = 0.0; /* FIXME: ?? */
    return GSL_ELOSS;
  }

  do {
    if(++i > 20000) {
      *prec   = 1.0;
      *result = sum;
      return GSL_ELOSS;
    }
    del *= (a+k)*(b+k) * x / ((c+k) * (k+1.0));  /* Gauss series */
    sum += del;
    delmax = locMAX(fabs(del), delmax);
    k += 1.0;
  } while(fabs(del/sum) > GSL_MACH_EPS);

  *prec   = (GSL_MACH_EPS*delmax)/fabs(sum) + i*GSL_MACH_EPS;
  *result = sum;
  return GSL_SUCCESS;
}


/* a = aR + i aI, b = aR - i aI */
static
int
hyperg_2F1_conj_series(const double aR, const double aI, const double c,
                       double x, double * result, double * prec)
{
  if(fabs(c) < locEPS) {
    *prec   = 1.0;
    *result = 0.0; /* FIXME: ?? */
    return GSL_ELOSS;
  }
  else {
    double sum    = 1.0;
    double del    = 1.0;
    double delmax = 1.0;
    double absdel;
    int k = 0;
    do {
      del *= ((aR+k)*(aR+k) + aI*aI)/((k+1.0)*(c+k)) * x;
      sum += del;
      absdel = fabs(del);
      delmax = locMAX(absdel, delmax);
      if(++k > 20000) {
        *prec   = 1.0;
        *result = sum;
        return GSL_ELOSS;
      }
    } while(fabs(del/sum) > 0.1*locEPS);
    *result = sum;
    *prec   = (GSL_MACH_EPS*delmax)/fabs(sum) + k*GSL_MACH_EPS;
    return GSL_SUCCESS;
  }
}


/* Luke's rational approximation. The most accesible
 * discussion is in [Kolbig, CPC 23, 51 (1981)].
 * The convergence is supposedly guaranteed for x < 0.
 * You have to read Luke's books to see this and other
 * results. Unfortunately, the stability is not so
 * clear to me, although is seems very efficient when
 * it works.
 */
static
int
hyperg_2F1_luke(const double a, const double b, const double c,
                const double xin, 
                double * result, double * prec)
{
  const double RECUR_BIG = 1.0e+50;
  const int nmax = 20000;
  int n = 3;
  const double x  = -xin;
  const double x3 = x*x*x;
  const double t0 = a*b/c;
  const double t1 = (a+1.0)*(b+1.0)/(2.0*c);
  const double t2 = (a+2.0)*(b+2.0)/(2.0*(c+1.0));
  double F = 1.0;

  double Bnm3 = 1.0;                                  /* B0 */
  double Bnm2 = 1.0 + t1 * x;                         /* B1 */
  double Bnm1 = 1.0 + t2 * x * (1.0 + t1/3.0 * x);    /* B2 */
 
  double Anm3 = 1.0;                                                      /* A0 */
  double Anm2 = Bnm2 - t0 * x;                                            /* A1 */
  double Anm1 = Bnm1 - t0*(1.0 + t2*x)*x + t0 * t1 * (c/(c+1.0)) * x*x;   /* A2 */

  while(1) {
    double npam1 = n + a - 1;
    double npbm1 = n + b - 1;
    double npcm1 = n + c - 1;
    double npam2 = n + a - 2;
    double npbm2 = n + b - 2;
    double npcm2 = n + c - 2;
    double tnm1  = 2*n - 1;
    double tnm3  = 2*n - 3;
    double tnm5  = 2*n - 5;
    double n2 = n*n;
    double F1 =  (3.0*n2 + (a+b-6)*n + 2 - a*b - 2*(a+b)) / (2*tnm3*npcm1);
    double F2 = -(3.0*n2 - (a+b+6)*n + 2 - a*b)*npam1*npbm1/(4*tnm1*tnm3*npcm2*npcm1);
    double F3 = (npam2*npam1*npbm2*npbm1*(n-a-2)*(n-b-2)) / (8*tnm3*tnm3*tnm5*(n+c-3)*npcm2*npcm1);
    double E  = -npam1*npbm1*(n-c-1) / (2*tnm3*npcm2*npcm1);

    double An = (1.0+F1*x)*Anm1 + (E + F2*x)*x*Anm2 + F3*x3*Anm3;
    double Bn = (1.0+F1*x)*Bnm1 + (E + F2*x)*x*Bnm2 + F3*x3*Bnm3;
    double r = An/Bn;

    *prec = fabs((F - r)/F);
    F = r;

    if(*prec < 0.1*locEPS || n > nmax) break;

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


/* Luke's rational approximation for the
 * case a = aR + i aI, b = aR - i aI.
 */
static
int
hyperg_2F1_conj_luke(const double aR, const double aI, const double c,
                     const double xin, 
                     double * result, double * prec)
{
  const double RECUR_BIG = 1.0e+50;
  const int nmax = 10000;
  int n = 3;
  const double x = -xin;
  const double x3 = x*x*x;
  const double atimesb = aR*aR + aI*aI;
  const double apb     = 2.0*aR;
  const double t0 = atimesb/c;
  const double t1 = (atimesb +     apb + 1.0)/(2.0*c);
  const double t2 = (atimesb + 2.0*apb + 4.0)/(2.0*(c+1.0));
  double F = 1.0;

  double Bnm3 = 1.0;                                  /* B0 */
  double Bnm2 = 1.0 + t1 * x;                         /* B1 */
  double Bnm1 = 1.0 + t2 * x * (1.0 + t1/3.0 * x);    /* B2 */
 
  double Anm3 = 1.0;                                                      /* A0 */
  double Anm2 = Bnm2 - t0 * x;                                            /* A1 */
  double Anm1 = Bnm1 - t0*(1.0 + t2*x)*x + t0 * t1 * (c/(c+1.0)) * x*x;   /* A2 */

  while(1) {
    double nm1 = n - 1;
    double nm2 = n - 2;
    double npam1_npbm1 = atimesb + nm1*apb + nm1*nm1;
    double npam2_npbm2 = atimesb + nm2*apb + nm2*nm2;
    double npcm1 = nm1 + c;
    double npcm2 = nm2 + c;
    double tnm1  = 2*n - 1;
    double tnm3  = 2*n - 3;
    double tnm5  = 2*n - 5;
    double n2 = n*n;
    double F1 =  (3.0*n2 + (apb-6)*n + 2 - atimesb - 2*apb) / (2*tnm3*npcm1);
    double F2 = -(3.0*n2 - (apb+6)*n + 2 - atimesb)*npam1_npbm1/(4*tnm1*tnm3*npcm2*npcm1);
    double F3 = (npam2_npbm2*npam1_npbm1*(nm2*nm2 - nm2*apb + atimesb)) / (8*tnm3*tnm3*tnm5*(n+c-3)*npcm2*npcm1);
    double E  = -npam1_npbm1*(n-c-1) / (2*tnm3*npcm2*npcm1);

    double An = (1.0+F1*x)*Anm1 + (E + F2*x)*x*Anm2 + F3*x3*Anm3;
    double Bn = (1.0+F1*x)*Bnm1 + (E + F2*x)*x*Bnm2 + F3*x3*Bnm3;
    double r = An/Bn;

    *prec = fabs(F - r)/fabs(F);
    F = r;

    if(*prec < 0.1*locEPS || n > nmax) break;

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


/* Do the reflection described in [Moshier, p. 334].
 * Assumes a,b,c != neg integer.
 */
static
int
hyperg_2F1_reflect(const double a, const double b, const double c,
                   const double x, double * result)
{
  const double d = c - a - b;
  const int intd  = floor(d+0.5);
  const int d_integer = ( fabs(d - intd) < locEPS );

  if(d_integer) {
    const double ln_omx = log(1.0 - x);
    const double ad = fabs(d);
    int stat_F2 = GSL_SUCCESS;
    double sgn_2;
    double F1, F2;
    double d1, d2;
    double lng_c;
    double lng_ad2;
    double lng_bd2;
    int stat_c;
    int stat_ad2;
    int stat_bd2;

    if(d >= 0.0) {
      d1 = d;
      d2 = 0.0;
    }
    else {
      d1 = 0.0;
      d2 = d;
    }

    stat_ad2 = gsl_sf_lngamma_impl(a+d2, &lng_ad2);
    stat_bd2 = gsl_sf_lngamma_impl(b+d2, &lng_bd2);
    stat_c   = gsl_sf_lngamma_impl(c,    &lng_c);

    /* Evaluate F1.
     */
    if(ad < locEPS) {
      /* d = 0 */
      F1 = 0.0;
    }
    else {
      double lng_ad;
      double lng_ad1;
      double lng_bd1;
      int stat_ad  = gsl_sf_lngamma_impl(ad,   &lng_ad);
      int stat_ad1 = gsl_sf_lngamma_impl(a+d1, &lng_ad1);
      int stat_bd1 = gsl_sf_lngamma_impl(b+d1, &lng_bd1);

      if(stat_ad1 == GSL_SUCCESS && stat_bd1 == GSL_SUCCESS) {
        /* Gamma functions in the denominator are ok.
	 * Proceed with evaluation.
	 */
	int i;
        double sum1 = 1.0;
        double term = 1.0;
        double ln_pre1 = lng_ad + lng_c + d2*ln_omx - lng_ad1 - lng_bd1;

        /* Do F1 sum.
         */
        for(i=1; i<ad; i++) {
	  int j = i-1;
	  term *= (a + d2 + j) * (b + d2 + j) / (1.0 + d2 + j) / i * (1.0-x);
          sum1 += term;
        }

        if(ln_pre1 + log(fabs(sum1)) < GSL_LOG_DBL_MAX) {
          F1 = exp(ln_pre1) * sum1;
        }
        else {
          *result = 0.0; /* FIXME: should be Inf */
          return GSL_EOVRFLW;
        }
      }
      else {
        /* Gamma functions in the denominator were not ok.
	 * So the F1 term is zero.
	 */
        F1 = 0.0;
      }
    } /* end F1 evaluation */


    /* Evaluate F2.
     */
    if(stat_ad2 == GSL_SUCCESS && stat_bd2 == GSL_SUCCESS) {
      /* Gamma functions in the denominator are ok.
       * Proceed with evaluation.
       */
      const int maxiter = 1000;
      int i;
      double psi_1 = -M_EULER;
      double psi_1pd; 
      double psi_apd1;
      double psi_bpd1;
      int stat_1pd  = gsl_sf_psi_impl(1.0 + ad, &psi_1pd);
      int stat_apd1 = gsl_sf_psi_impl(a + d1,   &psi_apd1);
      int stat_bpd1 = gsl_sf_psi_impl(b + d1,   &psi_bpd1);
      double psi  = psi_1 + psi_1pd - psi_apd1 - psi_bpd1 - ln_omx;
      double fact = 1.0;
      double sum2 = psi;
      double ln_pre2 = lng_c + d1*ln_omx - lng_ad2 - lng_bd2;
      
      /* Do F2 sum.
       */
      for(i=1; i<maxiter; i++) {
        int j = i-1;
        psi  += 1.0/(1.0+j) + 1.0/(1.0+ad+j) - 1.0/(a+d1+j) - 1.0/(b+d1+j);
        fact *= (a+d1+j)*(b+d1+j)/(ad+j)/i * (1.0-x);
	sum2 += fact * psi;
      }

      if(i == maxiter) stat_F2 = GSL_EMAXITER;

      if(ln_pre2 + log(fabs(sum2)) < GSL_LOG_DBL_MAX) {
        F2 = exp(ln_pre2) * sum2;
      }
      else {
        *result = 0.0; /* FIXME: should be Inf */
        return GSL_EOVRFLW;
      }
    }
    else {
      /* Gamma functions on the denominator not ok.
       * So the F2 term is zero.
       */
      F2 = 0.0;
    } /* end F2 evaluation */

    sgn_2 = ( GSL_IS_ODD(intd) ? -1.0 : 1.0 );
    *result = F1 + sgn_2 * F2;
    return stat_F2;
  }
  else {
    /* d not an integer */

    double pre1, pre2;
    double sgn1, sgn2;
    double F1, F2;
    double prec_F1, prec_F2;
    int status_F1, status_F2;

    /* These gamma functions appear in the denominator, so we
     * catch their harmless domain errors and set the terms to zero.
     */
    double ln_g1ca,  ln_g1cb,  ln_g2a,  ln_g2b;
    double sgn_g1ca, sgn_g1cb, sgn_g2a, sgn_g2b;
    int stat_1ca = gsl_sf_lngamma_sgn_impl(c-a, &ln_g1ca, &sgn_g1ca);
    int stat_1cb = gsl_sf_lngamma_sgn_impl(c-b, &ln_g1cb, &sgn_g1cb);
    int stat_2a  = gsl_sf_lngamma_sgn_impl(a, &ln_g2a, &sgn_g2a);
    int stat_2b  = gsl_sf_lngamma_sgn_impl(b, &ln_g2b, &sgn_g2b);
    int ok1 = (stat_1ca == GSL_SUCCESS && stat_1cb == GSL_SUCCESS);
    int ok2 = (stat_2a  == GSL_SUCCESS && stat_2b  == GSL_SUCCESS);
    
    double ln_gc,  ln_gd,  ln_gmd;
    double sgn_gc, sgn_gd, sgn_gmd;
    gsl_sf_lngamma_sgn_impl( c, &ln_gc,  &sgn_gc);
    gsl_sf_lngamma_sgn_impl( d, &ln_gd,  &sgn_gd);
    gsl_sf_lngamma_sgn_impl(-d, &ln_gmd, &sgn_gmd);
    
    sgn1 = sgn_gc * sgn_gd  * sgn_g1ca * sgn_g1cb;
    sgn2 = sgn_gc * sgn_gmd * sgn_g2a  * sgn_g2b;

    if(ok1 && ok2) {
      double ln_pre1 = ln_gc + ln_gd  - ln_g1ca - ln_g1cb;
      double ln_pre2 = ln_gc + ln_gmd - ln_g2a  - ln_g2b + d*log(1.0-x);
      if(ln_pre1 < GSL_LOG_DBL_MAX && ln_pre2 < GSL_LOG_DBL_MAX) {
        pre1 = sgn1 * exp(ln_pre1);
        pre2 = sgn2 * exp(ln_pre2);
      }
      else {
   	*result = 0.0; /* FIXME: should be Inf */
   	return GSL_EOVRFLW;
      }
    }
    else if(ok1 && !ok2) {
      double ln_pre1 = ln_gc + ln_gd  - ln_g1ca - ln_g1cb;
      if(ln_pre1 < GSL_LOG_DBL_MAX) {
        pre1 = sgn1 * exp(ln_pre1);
        pre2 = 0.0;
      }
      else {
        *result = 0.0; /* FIXME: should be Inf */
   	return GSL_EOVRFLW;
      }
    }
    else if(!ok1 && ok2) {
      double ln_pre2 = ln_gc + ln_gmd - ln_g2a  - ln_g2b + d*log(1.0-x);
      if(ln_pre2 < GSL_LOG_DBL_MAX) {
        pre1 = 0.0;
        pre2 = sgn2 * exp(ln_pre2);
      }
      else {
        *result = 0.0; /* FIXME: should be Inf */
   	return GSL_EOVRFLW;
      }
    }
    else {
      pre1 = 0.0;
      pre2 = 0.0;
      *result = 0.0;
      return GSL_EUNDRFLW;
    }

    status_F1 = hyperg_2F1_series(  a,   b, 1.0-d, 1.0-x, &F1, &prec_F1);
    status_F2 = hyperg_2F1_series(c-a, c-b, 1.0+d, 1.0-x, &F2, &prec_F2);

    *result = pre1*F1 + pre2*F2;

    if(prec_F1 > 10.0 * locEPS || prec_F2 > 10.0 * locEPS)
      return GSL_ELOSS;
    else
      return GSL_SUCCESS;
  }
}


static int pow_omx(const double x, const double p, double * result)
{
  double ln_omx;
  double ln_result;
  if(fabs(x) < GSL_ROOT3_MACH_EPS) {
    ln_omx = -x*(1.0 + 0.5*x + x*x/3.0);
  }
  else {
    ln_omx = log(1.0-x);
  }
  ln_result = p * ln_omx;
  if(ln_result > GSL_LOG_DBL_MAX) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(ln_result < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    *result = exp(ln_result);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_hyperg_2F1_impl(double a, double b, const double c,
                       const double x,
                       double * result)
{
  const double d = c - a - b;
  const double rinta = floor(a + 0.5);
  const double rintb = floor(b + 0.5);
  const double rintc = floor(c + 0.5);
  int a_neg_integer;
  int b_neg_integer;
  int c_neg_integer;

  if(x < -1.0 || 1.0 <= x) return GSL_EDOM;

  a_neg_integer = ( a < 0.0  &&  fabs(a - rinta) < locEPS );
  b_neg_integer = ( b < 0.0  &&  fabs(b - rintb) < locEPS );
  c_neg_integer = ( c < 0.0  &&  fabs(c - rintc) < locEPS );

  if(c_neg_integer) {
    if(! (a_neg_integer && a > c + 0.1)) return GSL_EDOM;
    if(! (b_neg_integer && b > c + 0.1)) return GSL_EDOM;
  }

  if(fabs(c-b) < locEPS || fabs(c-a) < locEPS) {
    return pow_omx(x, d, result);  /* (1-x)^(c-a-b) */
  }

  if(fabs(a) < 10.0 && fabs(b) < 10.0) {
    /* a and b are not too large, so we attempt
     * variations on the series summation.
     */
    if(a_neg_integer) {
      double prec;
      return hyperg_2F1_series(rinta, b, c, x, result, &prec);
    }
    if(b_neg_integer) {
      double prec;
      return hyperg_2F1_series(a, rintb, c, x, result, &prec);
    }

    if(x < -0.25) {
      double prec;
      return hyperg_2F1_luke(a, b, c, x, result, &prec);
    }
    else if(x < 0.5) {
      double prec;
      return hyperg_2F1_series(a, b, c, x, result, &prec);
    }
    else {
      if(fabs(c) > 10.0) {
        double prec;
        return hyperg_2F1_series(a, b, c, x, result, &prec);
      }
      else {
        return hyperg_2F1_reflect(a, b, c, x, result);
      }
    }
  }
  else {
    /* Either a or b or both large.
     * Introduce some new variables ap,bp so that bp is
     * the larger in magnitude.
     */
    double ap, bp; 
    if(fabs(a) > fabs(b)) {
      bp = a;
      ap = b;
    }
    else {
      bp = b;
      ap = a;
    }

    if(x < 0.0) {
      /* What the hell, maybe Luke will converge.
       */
      double prec;
      return hyperg_2F1_luke(a, b, c, x, result, &prec);
    }

    if(locMAX(fabs(a),1.0)*fabs(bp)*fabs(x) < 2.0*fabs(c)) {
      /* If c is large enough or x is small enough,
       * we can attempt the series anyway.
       */
      double prec;
      return hyperg_2F1_series(a, b, c, x, result, &prec);
    }

    if(fabs(bp*bp*x*x) < 0.001*fabs(bp) && fabs(a) < 10.0) {
      /* The famous but nearly worthless "large b" asymptotic.
       */
      return gsl_sf_hyperg_1F1_impl(a, c, bp*x, result);
    }

    /* We give up. */
    *result = 0.0;
    return GSL_EUNIMPL;
  }
}


int gsl_sf_hyperg_2F1_conj_impl(const double aR, const double aI, const double c,
                                const double x,
				double * result)
{
  const double ax = fabs(x);
  const double rintc = floor(c + 0.5);
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rintc) < locEPS );

  if(ax >= 1.0) return GSL_EDOM;
  if(c_neg_integer) return GSL_EDOM;

  if(ax < 0.25 && fabs(aR) < 20.0 && fabs(aI) < 20.0) {
    double prec;
    int status = hyperg_2F1_conj_series(aR, aI, c, x, result, &prec);
    return status;
  }
  else if(fabs(aR) < 10.0 && fabs(aI) < 10.0) {
    if(x < -0.25) {
      double prec;
      return hyperg_2F1_conj_luke(aR, aI, c, x, result, &prec);
    }
    else {
      double prec;
      return hyperg_2F1_conj_series(aR, aI, c, x, result, &prec);
    }
  }
  else {
  
    if(x < 0.0) {
      /* What the hell, maybe Luke will converge.
       */
      double prec;
      return hyperg_2F1_conj_luke(aR, aI, c, x, result, &prec); 
    }

    /* Give up. */
    *result = 0.0;
    return GSL_EUNIMPL;
  }
}


int gsl_sf_hyperg_2F1_renorm_impl(const double a, const double b, const double c,
                                  const double x,
			          double * result
			          )
{
  const double rinta = floor(a + 0.5);
  const double rintb = floor(b + 0.5);
  const double rintc = floor(c + 0.5);
  const int a_neg_integer = ( a < 0.0  &&  fabs(a - rinta) < locEPS );
  const int b_neg_integer = ( b < 0.0  &&  fabs(b - rintb) < locEPS );
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rintc) < locEPS );
  
  if(c_neg_integer) {
    if((a_neg_integer && a > c+0.1) || (b_neg_integer && b > c+0.1)) {
      /* 2F1 terminates early */
      *result = 0.0;
      return GSL_SUCCESS;
    }
    else {
      /* 2F1 does not terminate early enough, so something survives */
      /* [Abramowitz+Stegun, 15.1.2] */
      double g1, g2, g3, g4, g5;
      double s1, s2, s3, s4, s5;
      int stat = 0;
      stat += gsl_sf_lngamma_sgn_impl(a-c+1, &g1, &s1);
      stat += gsl_sf_lngamma_sgn_impl(b-c+1, &g2, &s2);
      stat += gsl_sf_lngamma_sgn_impl(a, &g3, &s3);
      stat += gsl_sf_lngamma_sgn_impl(b, &g4, &s4);
      stat += gsl_sf_lngamma_sgn_impl(-c+2, &g5, &s5);
      if(stat != 0) {
        *result = 0.0;
        return GSL_EDOM;
      }
      else {
        double F;
        int stat_F = gsl_sf_hyperg_2F1_impl(a-c+1, b-c+1, -c+2, x, &F);
        double ln_pre = g1 + g2 - g3 - g4 - g5;
	double sg = s1 * s2 * s3 * s4 * s5;
	double sF = ( F > 0.0 ? 1.0 : -1.0 );
	double lnr = ln_pre + (1.0-c)*log(x) + log(fabs(F));
	if(lnr < GSL_LOG_DBL_MAX) {
          *result = sg * sF * exp(lnr);
          return GSL_SUCCESS;
	}
	else {
	  *result = 0.0;
	  return GSL_EOVRFLW;
	}
      }
    }
  }
  else {
    /* generic c */
    double F;
    double lng;
    int stat_g = gsl_sf_lngamma_impl(c, &lng);
    int stat_F = gsl_sf_hyperg_2F1_impl(a, b, c, x, &F);
    if(stat_F == GSL_SUCCESS && stat_g == GSL_SUCCESS) {
      double ln_absF   = log(fabs(F));
      double ln_result = ln_absF - lng;
      if(ln_result < GSL_LOG_DBL_MIN) {
        *result = 0.0;
        return GSL_EUNDRFLW;
      }
      else {
        double sgn = (F > 0.0 ? 1.0 : -1.0);
        *result = sgn * exp(ln_result);
        return GSL_SUCCESS;
      }
    }
    else {
      *result = 0.0;
      return GSL_EFAILED;
    }
  }
}

int gsl_sf_hyperg_2F1_conj_renorm_impl(const double aR, const double aI, const double c,
                                       const double x,
			               double * result
			               )
{
/* FIXME: copy above, when it works right */
}


int
gsl_sf_hyperg_2F1_e(double a, double b, double c, double x, double * result)
{
  int status = gsl_sf_hyperg_2F1_impl(a, b, c, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hyperg_2F1_e", status);
  }
  return status;
}

int
gsl_sf_hyperg_2F1_conj_e(double aR, double aI, double c, double x, double * result)
{
  int status = gsl_sf_hyperg_2F1_conj_impl(aR, aI, c, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hyperg_2F1_conj_e", status);
  }
  return status;
}

double
gsl_sf_hyperg_2F1(double a, double b, double c, double x)
{
  double y;
  int status = gsl_sf_hyperg_2F1_impl(a, b, c, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hyperg_2F1", status);
  }
  return y;
}

double
gsl_sf_hyperg_2F1_conj(double aR, double aI, double c, double x)
{
  double y;
  int status = gsl_sf_hyperg_2F1_conj_impl(aR, aI, c, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hyperg_2F1_conj", status);
  }
  return y;
}
