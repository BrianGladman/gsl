/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
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
static int hyperg_2F1_conj_series(const double aR, const double aI, const double c,
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
    } while(fabs(del/sum) > locEPS);
    *result = sum;
    *prec   = (GSL_MACH_EPS*delmax)/fabs(sum) + k*GSL_MACH_EPS;
    return GSL_SUCCESS;
  }
}

static int hyperg_2F1_luke(const double a, const double b, const double c,
                           const double x, 
			   double * result, double * prec)
{
  const int nmax = 10000;
  int n = 3;
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
    double E  = npam1*npbm1*(n-c-1) / (2*tnm3*npcm2*npcm1);

    double An = (1.0+F1*x)*Anm1 + (E + F2*x)*x*Anm2 + F3*x3*Anm3;
    double Bn = (1.0+F1*x)*Bnm1 + (E + F2*x)*x*Bnm2 + F3*x3*Bnm3;
    double r = An/Bn;

    *prec = fabs(F - r)/F;
    F = r;

    if(*prec < locEPS || n > nmax) break;

    n++;
    Bnm3 = Bnm2;
    Bnm2 = Bnm1;
    Bnm1 = Bn;
    Anm3 = Anm2;
    Anm2 = Anm1;
    Anm1 = An;
  }
  
  *result = F;

  if(*prec > locEPS)
    return GSL_ELOSS;
  else
    return GSL_SUCCESS;
}

/* a = aR + i aI, b = aR - i aI */
static int hyperg_2F1_conj_luke(const double aR, const double aI, const double c,
                                const double x, 
			        double * result, double * prec)
{
  const int nmax = 10000;
  int n = 3;
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
    double E  = npam1_npbm1*(n-c-1) / (2*tnm3*npcm2*npcm1);

    double An = (1.0+F1*x)*Anm1 + (E + F2*x)*x*Anm2 + F3*x3*Anm3;
    double Bn = (1.0+F1*x)*Bnm1 + (E + F2*x)*x*Bnm2 + F3*x3*Bnm3;
    double r = An/Bn;

    *prec = fabs(F - r)/F;
    F = r;

    if(*prec < locEPS || n > nmax) break;

    n++;
    Bnm3 = Bnm2;
    Bnm2 = Bnm1;
    Bnm1 = Bn;
    Anm3 = Anm2;
    Anm2 = Anm1;
    Anm1 = An;
  }
  
  *result = F;

  if(*prec > locEPS)
    return GSL_ELOSS;
  else
    return GSL_SUCCESS;
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

int gsl_sf_hyperg_2F1_impl(const double a, const double b, const double c,
                           const double x,
			   double * result
			   )
{
  const double ax = fabs(x);
  const double d  = c - a - b;
  const int a_neg_integer = ( a < 0.0  &&  fabs(a - rint(a)) < locEPS );
  const int b_neg_integer = ( b < 0.0  &&  fabs(b - rint(b)) < locEPS );
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rint(c)) < locEPS );
  const int d_integer     = ( fabs(d - rint(d)) < locEPS );
  
  if(ax >= 1.0) return GSL_EDOM;

  if(c_neg_integer) {
    if(! (a_neg_integer && a > c)) return GSL_EDOM;
    if(! (b_neg_integer && b > c)) return GSL_EDOM;
  }
  
  if(fabs(c-b) < locEPS || fabs(c-a) < locEPS) {
    return pow_omx(x, c-a-b, result);  /* (1-x)^(c-a-b) */
  }
  
  if(x < -0.5 && b > 0.0) {
    double F, F_prec, p;
    int status_F = hyperg_2F1_series(a, c-b, c, -x/(1.0-x), &F, &F_prec);
    int status_p = pow_omx(x, -a, &p);
    if(status_p == GSL_SUCCESS && status_F == GSL_SUCCESS) {
      *result = p * F;
      return GSL_SUCCESS; 
    }
    else if(status_p == GSL_EUNDRFLW) {
      *result = 0.0;
      return GSL_EUNDRFLW;
    }
    else if(status_p == GSL_EOVRFLW) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EOVRFLW;
    }
    else if(F_prec > locEPS) {
      *result = p * F;
      return GSL_ELOSS;
    }
    else {
      *result = 0.0;
      return GSL_EFAILED;
    }
  }
  
  if(x < -0.5 && b <= 0.0) {
    double F, F_prec, p;
    int status_F = hyperg_2F1_series(c-a, b, c, -x/(1.0-x), &F, &F_prec);
    int status_p = pow_omx(x, -b, &p);
    if(status_p == GSL_SUCCESS && status_F == GSL_SUCCESS) {
      *result = p * F;
      return GSL_SUCCESS; 
    }
    else if(status_p == GSL_EUNDRFLW) {
      *result = 0.0;
      return GSL_EUNDRFLW;
    }
    else if(status_p == GSL_EOVRFLW) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EOVRFLW;
    }
    else if(F_prec > locEPS) {
      *result = p * F;
      return GSL_ELOSS;
    }
    else {
      *result = 0.0;
      return GSL_EFAILED;
    }
  }
  
  if(x > 0.75 && !d_integer) {
    double pre1, pre2, F1, F2;
    double status_F1, status_F2;
    double prec_F1, prec_F2;
    double lng_c   = gsl_sf_lngamma(c);
    double ln_pre1 = lng_c + gsl_sf_ln_gamma( d) - gsl_sf_lngamma(c-a) - gsl_sf_lngamma(c-b);
    double ln_pre2 = lng_c + gsl_sf_ln_gamma(-d) - gsl_sf_lngamma(a)   - gsl_sf_lngamma(b) + d*log(1.0-x);
    if(ln_pre1 > GSL_LOG_DBL_MAX || ln_pre2 > GSL_LOG_DBL_MAX) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EOVRFLW;
    }
    pre1 = exp(ln_pre1);
    pre2 = exp(ln_pre2);
    status_F1 = hyperg_2F1_series(  a,   b, 1.0-d, 1.0-x, &F1, &prec_F1);
    status_F2 = hyperg_2F1_series(c-a, c-b, 1.0+d, 1.0-x, &F2, &prec_F2);
    *result = pre1*F1 + pre2*F2;
    if(prec_F1 > locEPS || prec_F2 > locEPS) {
      return GSL_ELOSS;
    }
    else {
      return GSL_SUCCESS;
    }
  }
  
  if(x > 0.75 && d_integer) {
    double ad = fabs(d);
    double ln_pre1, ln_pre2, pre1, pre2;
    double F1, F2;
    double d1, d2;
    if(d >= 0.0) {
      d1 = d;
      d2 = 0.0;
    }
    else {
      d1 = 0.0;
      d2 = d;
    }
    ln_pre1 = gsl_sf_lngamma(ad) + gsl_sf_lngamma(c) + d2*log(1.0-x) - gsl_sf_lngamma(a+d1) - gsl_sf_lngamma(b+d1);
    ln_pre2 = gsl_sf_lngamma(c) + d1*log(1.0-x) - gsl_sf_lngamma(a+d2) - gsl_sf_lngamma(b+d2);
    if(ln_pre1 > GSL_LOG_DBL_MAX || ln_pre2 > GSL_LOG_DBL_MAX) {
      *result = 0.0; /* FIXME: should be Inf (?) */
      return GSL_EOVRFLW;
    }
    pre1 = exp(ln_pre1);
    pre2 = exp(ln_pre2);
    
    /* FIXME: F1= ...  F2= ... */
  }
  
  /* default case */
  {
    double prec;
    int status = hyperg_2F1_series(a, b, c, x, result, &prec);
    if(prec > locEPS) {
      return GSL_ELOSS;
    }
    else {
      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_hyperg_2F1_conj_impl(const double aR, const double aI, const double c,
                                const double x,
				double * result)
{
  double ax = fabs(x);
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rint(c)) < locEPS );

  if(ax >= 1.0) return GSL_EDOM;
  if(c_neg_integer) return GSL_EDOM;

  if(fabs(aR) < 5.0 && fabs(aI) < 5.0) {
    double prec;
    int status = hyperg_2F1_conj_series(aR, aI, c, x, result, &prec);
    return status;
  }
  else {
  }
}
