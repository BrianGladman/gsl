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
    } while(fabs(del/sum) > locEPS);
    *result = sum;
    *prec   = (GSL_MACH_EPS*delmax)/fabs(sum) + k*GSL_MACH_EPS;
    return GSL_SUCCESS;
  }
}


static
int
hyperg_2F1_luke(const double a, const double b, const double c,
                const double xin, 
                double * result, double * prec)
{
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

    if(*prec < locEPS || n > nmax) break;

#if defined(LUKE_INSTRUMENT)
    printf("%20.15g   %20.15g", An, Bn);
#endif

    if(fabs(An) > 1.e5 || fabs(Bn) > 1.e5) {
      An /= 1.e5;
      Bn /= 1.e5;
      Anm1 /= 1.e5;
      Bnm1 /= 1.e5;
      Anm2 /= 1.e5;
      Bnm2 /= 1.e5;
      Anm3 /= 1.e5;
      Bnm3 /= 1.e5;

#if defined(LUKE_INSTRUMENT)      
      printf("   -");
#endif
    }
    else if(fabs(An) < 1.e-5 || fabs(Bn) < 1.e-5) {
      An *= 1.e5;
      Bn *= 1.e5;
      Anm1 *= 1.e5;
      Bnm1 *= 1.e5;
      Anm2 *= 1.e5;
      Bnm2 *= 1.e5;
      Anm3 *= 1.e5;
      Bnm3 *= 1.e5;
      
#if defined(LUKE_INSTRUMENT) 
      printf("   +");
#endif
    }

#if defined(LUKE_INSTRUMENT)
    printf("  %20.15g  \n", An/Bn);
#endif

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

static void do_it(double a, double b, double c, double x)
{
  double y;
  double prec;
  int stat = gsl_sf_hyperg_2F1_impl(a, b, c, x, &y);
  printf("%5.3g  %5.3g  %5.3g   %10.6g", a, b, c, x);
  printf("   ");
  printf("%24.18g  %8.4g  %s", y, prec, gsl_strerror(stat));
  printf("\n");
}

void testy(void)
{
  do_it(1.0, 1.0, 1.0, -1.0);
  do_it(1.0, 1.0, 1.0, -0.8);
  do_it(1.0, 1.0, 1.0, -0.2);
  do_it(1.0, 1.0, 1.0, -1.e-10);
  do_it(1.0, 1.0, 1.0,  1.e-10);
  do_it(1.0, 1.0, 1.0,  0.2);
  do_it(1.0, 1.0, 1.0,  0.8);

  do_it(10.0, 1.0, 1.0, -1.0);
  do_it(10.0, 1.0, 1.0, -0.8);
  do_it(10.0, 1.0, 1.0, -0.2);
  do_it(10.0, 1.0, 1.0, -1.e-10);
  do_it(10.0, 1.0, 1.0,  1.e-10);
  do_it(10.0, 1.0, 1.0,  0.2);
  do_it(10.0, 1.0, 1.0,  0.8);
  
  do_it(50.0, 1.0, 1.0, -1.0);
  do_it(50.0, 1.0, 1.0, -0.8);
  do_it(50.0, 1.0, 1.0, -0.2);
  do_it(50.0, 1.0, 1.0, -1.e-10);
  do_it(50.0, 1.0, 1.0,  1.e-10);
  do_it(50.0, 1.0, 1.0,  0.2);
  do_it(50.0, 1.0, 1.0,  0.8);
  
  do_it(100.0, 1.0, 1.0, -1.0);
  do_it(100.0, 1.0, 1.0, -0.8);
  do_it(100.0, 1.0, 1.0, -0.2);
  do_it(100.0, 1.0, 1.0, -1.e-10);
  do_it(100.0, 1.0, 1.0,  1.e-10);
  do_it(100.0, 1.0, 1.0,  0.2);
  do_it(100.0, 1.0, 1.0,  0.8);

  do_it(1.0, 1.0, 50.0, -1.0);
  do_it(1.0, 1.0, 50.0, -0.8);
  do_it(1.0, 1.0, 50.0, -0.2);
  do_it(1.0, 1.0, 50.0, -1.e-10);
  do_it(1.0, 1.0, 50.0,  1.e-10);
  do_it(1.0, 1.0, 50.0,  0.2);
  do_it(1.0, 1.0, 50.0,  0.8);
  
  do_it(50.0, 50.0, 1.0, -1.0);
  do_it(50.0, 50.0, 1.0, -0.8);
  do_it(50.0, 50.0, 1.0, -0.2);
  do_it(50.0, 50.0, 1.0, -1.e-10);
  do_it(50.0, 50.0, 1.0,  1.e-10);
  do_it(50.0, 50.0, 1.0,  0.2);
  do_it(50.0, 50.0, 1.0,  0.8);
  

  exit(0);
}


/* a = aR + i aI, b = aR - i aI */
static
int
hyperg_2F1_conj_luke(const double aR, const double aI, const double c,
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
    double E  = -npam1_npbm1*(n-c-1) / (2*tnm3*npcm2*npcm1);

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


int
gsl_sf_hyperg_2F1_impl(double a, double b, const double c,
                       const double x,
                       double * result)
{
  int a_neg_integer;
  int b_neg_integer;
  int c_neg_integer;
  
  if(x < -1.0 || 1.0 <= x) return GSL_EDOM;

  a_neg_integer = ( a < 0.0  &&  fabs(a - rint(a)) < locEPS );
  b_neg_integer = ( b < 0.0  &&  fabs(b - rint(b)) < locEPS );
  c_neg_integer = ( c < 0.0  &&  fabs(c - rint(c)) < locEPS );

  if(c_neg_integer) {
    if(! (a_neg_integer && a > c + 0.1)) return GSL_EDOM;
    if(! (b_neg_integer && b > c + 0.1)) return GSL_EDOM;
  }  
  if(fabs(c-b) < locEPS || fabs(c-a) < locEPS) {
    return pow_omx(x, c-a-b, result);  /* (1-x)^(c-a-b) */
  }

  if(fabs(a) < 10.0 && fabs(b) < 10.0) {
    double prec;
    return hyperg_2F1_series(a, b, c, x, result, &prec);
  }
  else if(x < 0.0) {
    double prec;
    return hyperg_2F1_luke(a, b, c, x, result, &prec);
  }
  else {
    double xi  = x/(x-1.0);
    double pre;
    double prec;
    int stat_pre = pow_omx(x, -a, &pre);
    if(stat_pre == GSL_EUNDRFLW) {
      *result = 0.0;
      return stat_pre;
    }
    else if(stat_pre == GSL_EOVRFLW) {
      *result = 0.0; /* FIXME: should be Inf */
      return stat_pre;
    }
    else {
      double newF;
      int stat_hyp = hyperg_2F1_luke(a, c-b, c, xi, &newF, &prec);
      *result = pre * newF;
      return stat_hyp;
    }
  }
}


int gsl_sf_hyperg_2F1_impl_old(double a, double b, const double c,
                           const double x,
                           double * result
                           )
{
  const double ax = fabs(x);
  double d;
  int a_neg_integer;
  int b_neg_integer;
  int c_neg_integer;
  int d_integer;
  
  if(ax >= 1.0) return GSL_EDOM;

  if(fabs(a) > fabs(b)) {
    double tmp = a;
    a = b;
    b = tmp;
  }

  d  = c - a - b;
  a_neg_integer = ( a < 0.0  &&  fabs(a - rint(a)) < locEPS );
  b_neg_integer = ( b < 0.0  &&  fabs(b - rint(b)) < locEPS );
  c_neg_integer = ( c < 0.0  &&  fabs(c - rint(c)) < locEPS );
  d_integer	= ( fabs(d - rint(d)) < locEPS );

  if(c_neg_integer) {
    if(! (a_neg_integer && a > c + 0.1)) return GSL_EDOM;
    if(! (b_neg_integer && b > c + 0.1)) return GSL_EDOM;
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
    double prec_F1, prec_F2;
    int status_F1, status_F2;
    {
      /* these gamma functions appear in the denominator, so we
         catch their harmless domain errors and set the terms to zero
      */
      double ln_g1ca, ln_g1cb, ln_g2a, ln_g2b;
      int stat_1ca = gsl_sf_lngamma_impl(c-a, &ln_g1ca);
      int stat_1cb = gsl_sf_lngamma_impl(c-b, &ln_g1cb);
      int stat_2a  = gsl_sf_lngamma_impl(a, &ln_g2a);
      int stat_2b  = gsl_sf_lngamma_impl(b, &ln_g2b);
      int ok1 = (stat_1ca == GSL_SUCCESS && stat_1cb == GSL_SUCCESS);
      int ok2 = (stat_2a  == GSL_SUCCESS && stat_2b  == GSL_SUCCESS);
      
      double ln_gc, ln_gd, ln_gmd;
      gsl_sf_lngamma_impl( c, &ln_gc);
      gsl_sf_lngamma_impl( d, &ln_gd);
      gsl_sf_lngamma_impl(-d, &ln_gmd);
      
      if(ok1 && ok2) {
	double ln_pre1 = ln_gc + ln_gd  - ln_g1ca - ln_g1cb;
        double ln_pre2 = ln_gc + ln_gmd - ln_g2a  - ln_g2b + d*log(1.0-x);
	if(ln_pre1 > GSL_LOG_DBL_MAX || ln_pre2 > GSL_LOG_DBL_MAX) {
          *result = 0.0; /* FIXME: should be Inf */
          return GSL_EOVRFLW;
        }
        pre1 = exp(ln_pre1);
        pre2 = exp(ln_pre2);
      }
      else if(ok1 && !ok2) {
	double ln_pre1 = ln_gc + ln_gd  - ln_g1ca - ln_g1cb;
	if(ln_pre1 > GSL_LOG_DBL_MAX) {
	  *result = 0.0; /* FIXME: should be Inf */
          return GSL_EOVRFLW;
	}
	pre1 = exp(ln_pre1);
	pre2 = 0.0;
      }
      else if(!ok1 && ok2) {
	double ln_pre2 = ln_gc + ln_gmd - ln_g2a  - ln_g2b + d*log(1.0-x);
	if(ln_pre2 > GSL_LOG_DBL_MAX) {
	  *result = 0.0; /* FIXME: should be Inf */
          return GSL_EOVRFLW;
	}
        pre1 = 0.0;
	pre2 = exp(ln_pre2);
      }
      else {
        pre1 = 0.0;
	pre2 = 0.0;
	*result = 0.0;
	return GSL_EUNDRFLW;
      }
    }
    status_F1 = hyperg_2F1_series(  a,   b, 1.0-d, 1.0-x, &F1, &prec_F1);
    status_F2 = hyperg_2F1_series(c-a, c-b, 1.0+d, 1.0-x, &F2, &prec_F2);
    *result = pre1*F1 + pre2*F2;
    if(prec_F1 > locEPS || prec_F2 > locEPS)
      return GSL_ELOSS;
    else
      return GSL_SUCCESS;
  }
  
  if(x > 0.75 && d_integer) {
    int i;
    double ad = fabs(d);
    double ln_pre1, ln_pre2, pre1, pre2;
    double F1, F2;
    double d1, d2;
    double sum1 = 0.0, sum2 = 0.0;
    double lng_ad, lng_c;
    double lng_ad1, lng_ad2, lng_bd1, lng_bd2;
    int stat_ad, stat_c;
    int stat_ad1, stat_ad2, stat_bd1, stat_bd2;
    double ln_omx = log(1.0 - x);
    if(d >= 0.0) {
      d1 = d;
      d2 = 0.0;
    }
    else {
      d1 = 0.0;
      d2 = d;
    }
    stat_ad  = gsl_sf_lngamma_impl(ad,   &lng_ad);
    stat_c   = gsl_sf_lngamma_impl(c,    &lng_c);
    stat_ad1 = gsl_sf_lngamma_impl(a+d1, &lng_ad1);
    stat_ad2 = gsl_sf_lngamma_impl(a+d1, &lng_ad2);
    stat_bd1 = gsl_sf_lngamma_impl(a+d1, &lng_bd1);
    stat_bd2 = gsl_sf_lngamma_impl(a+d1, &lng_bd2);
    
    if(stat_ad != GSL_SUCCESS || stat_c != GSL_SUCCESS)
    ln_pre1 = lng_ad + lng_c + d2*ln_omx - lng_ad1 - lng_bd1;
    ln_pre2 =          lng_c + d1*ln_omx - lng_ad2 - lng_bd2;

    if(ln_pre1 > GSL_LOG_DBL_MAX || ln_pre2 > GSL_LOG_DBL_MAX) {
      *result = 0.0; /* FIXME: should be Inf (?) */
      return GSL_EOVRFLW;
    }
    pre1 = exp(ln_pre1);
    pre2 = exp(ln_pre2);
    
    {
      /* F1 sum */
      double term = 1.0;;
      for(i=0; i<ad; i++) {
        sum1 += term;
	term *= (a + d2 + i) * (b + d2 + i) / (1.0 + d2 + i) / (i+1) * (1.0-x);
      }
    }
    
    {
      /* F2 sum */
      double term;
      double Psi;
      double ln_dkfact;
      double p1 = -M_EULER;
      double p2, p3, p4;
      gsl_sf_lnfact_impl(ad, &ln_dkfact);
      term = exp(-ln_dkfact);
      gsl_sf_psi_impl(1 + ad, &p2);
      gsl_sf_psi_impl(a + d1, &p3);
      gsl_sf_psi_impl(b + d1, &p4);
      Psi = p1 + p2 - p3 - p4 - log(1.0-x);
      for(i=0; i<ad; i++) {
        sum2 += term * Psi;
	term *= (a + d1 + i) * (b + d1 + i) / (ad + i + 1) / (i+1) * (1.0-x);
	Psi  +=  1.0/(2.0 + i) + 1.0/(2.0 + ad + i)
	        -1.0/(a + d1 + 1 + i) - 1.0/(b + d1 + 1 + i);
      }
    }
    
    F1 = pre1 * sum1;
    F2 = pre2 * sum2;
    
    
  }
  

  {
    /* default case */
    double prec;
    int status = hyperg_2F1_series(a, b, c, x, result, &prec);
    if(prec > locEPS)
      return GSL_ELOSS;
    else
      return GSL_SUCCESS;
  }
}


int gsl_sf_hyperg_2F1_conj_impl(const double aR, const double aI, const double c,
                                const double x,
				double * result)
{
  const double ax = fabs(x);
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rint(c)) < locEPS );

  if(ax >= 1.0) return GSL_EDOM;
  if(c_neg_integer) return GSL_EDOM;

  if(ax < 0.2 && fabs(aR) < 20.0 && fabs(aI) < 20.0) {
    double prec;
    int status = hyperg_2F1_conj_series(aR, aI, c, x, result, &prec);
    return status;
  }
  else if(aR*aR+aI*aI < 50.0 /* fabs(aR) < 5.0 && fabs(aI) < 5.0 */) {
    double prec;
    int status = hyperg_2F1_conj_series(aR, aI, c, x, result, &prec);
    return status;
  }
  else {
    double prec;
    int status = hyperg_2F1_conj_luke(aR, aI, c, x, result, &prec);
    return status;
  }
}


int gsl_sf_hyperg_2F1_renorm_impl(const double a, const double b, const double c,
                                  const double x,
			          double * result
			          )
{
  const int a_neg_integer = ( a < 0.0  &&  fabs(a - rint(a)) < locEPS );
  const int b_neg_integer = ( b < 0.0  &&  fabs(b - rint(b)) < locEPS );
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rint(c)) < locEPS );
  
  if(c_neg_integer) {
    if((a_neg_integer && a > c) || (b_neg_integer && b > c)) {
      /* 2F1 terminates early */
      *result = 0.0;
      return GSL_SUCCESS;
    }
    else {
      /* 2F1 does not terminate early enough, so something survives */
      /* [Abramowitz+Stegun, 15.1.2] */
      double g1, g2, g3, g4, g5;
      int stat = 0;
      stat += gsl_sf_lngamma_impl(a-c+1, &g1);
      stat += gsl_sf_lngamma_impl(b-c+1, &g2);
      stat += gsl_sf_lngamma_impl(a, &g3);
      stat += gsl_sf_lngamma_impl(b, &g4);
      stat += gsl_sf_lngamma_impl(-c+2, &g5);
      if(stat != 0) {
        *result = 0.0;
        return GSL_EDOM;
      }
      else {
        double F;
        int stat_F = gsl_sf_hyperg_2F1_impl(a-c+1, b-c+1, -c+2, x, &F);
        double ln_pre = g1 + g2 - g3 - g4 - g5;
        /* FIXME: error handling */
        *result = exp(ln_pre) * gsl_sf_pow_int(x, -c+1) * F;
        return GSL_SUCCESS;
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
