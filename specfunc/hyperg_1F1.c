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
#define locMIN(a,b)     ((a) < (b) ? (a) : (b))
#define locEPS          (1000.0*GSL_MACH_EPS)



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

  int stat_b   = gsl_sf_lngamma_sgn_impl(b,   &lg_b,   &sgn_b);
  int stat_bma = gsl_sf_lngamma_sgn_impl(b-a, &lg_bma, &sgn_bma);
  
  if(stat_b == GSL_SUCCESS && stat_bma == GSL_SUCCESS) {
    double prec_F;
    double F;
    int stat_F = gsl_sf_hyperg_2F0_series_impl(a, 1.0+a-b, -1.0/x, -1, &F, &prec_F);
    if((stat_F == GSL_SUCCESS || stat_F == GSL_ELOSS) && F != 0) {
      double ln_F = log(fabs(F));
      double ln_pre = lg_b - lg_bma - a*log(-x);
      int stat_e = gsl_sf_exp_sgn_impl(ln_F + ln_pre, sgn_bma*sgn_b*F, result);
      if(stat_e == GSL_SUCCESS)
        return stat_F;
      else
        return stat_e;
    }
    else {
      *result = 0.0;
      return stat_F;
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

  int stat_b = gsl_sf_lngamma_sgn_impl(b, &lg_b, &sgn_b);
  int stat_a = gsl_sf_lngamma_sgn_impl(a, &lg_a, &sgn_a);

  if(stat_a == GSL_SUCCESS && stat_b == GSL_SUCCESS) {
    double prec_F;
    double F;
    int stat_F = gsl_sf_hyperg_2F0_series_impl(b-a, 1.0-a, 1.0/x, -1, &F, &prec_F);
    if((stat_F == GSL_SUCCESS || stat_F == GSL_ELOSS) && F != 0) {
      double ln_F = log(fabs(F));
      double ln_pre = lg_b - lg_a + x + (a-b)*log(x);
      int stat_e = gsl_sf_exp_sgn_impl(ln_F + ln_pre, sgn_a*sgn_b*F, result);
      if(stat_e == GSL_SUCCESS)
        return stat_F;
      else
        return stat_e;
    }
    else {
      *result = 0.0;
      return stat_F;
    }
  }
  else {
    *result = 0.0;
    return GSL_EDOM;
  }
}


/* Asymptotic result for x < 2b-4a, 2b-4a large.
 * [Abramowitz+Stegun, 13.5.21]
 */
static
int
hyperg_1F1_large2bm4a(const double a, const double b, const double x, double * result)
{
  double eta    = 2.0*b - 4.0*a;
  double cos2th = x/eta;
  double sin2th = 1.0 - cos2th;
  double th = acos(sqrt(cos2th));
  double pre_h  = 0.25*M_PI*M_PI*eta*eta*cos2th*sin2th;
  double ser;
  double lnpre;
  double lnr;
  double lg_b;
  gsl_sf_lngamma_impl(b, &lg_b);
  lnpre = lg_b + 0.5*x + 0.5*(1.0-b)*log(0.25*x*eta) - 0.25*log(pre_h);
  ser = sin(a*M_PI) + sin(0.25*eta*(2.0*th - sin(2.0*th)) + 0.25*M_PI);
  lnr = lnpre + log(fabs(ser));
  return gsl_sf_exp_sgn_impl(lnr, ser, result);
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


/* Series for 1F1(1,b,x)
 * b > 0
 */
static
int
hyperg_1F1_1_series(const double b, const double x, double * result)
{
  double term = 1.0;
  double sum  = 1.0;
  double n    = 1.0;
  while(fabs(term/sum) > 10.0*GSL_MACH_EPS) {
    term *= x/(b+n-1);
    sum  += term;
    n += 1.0;
  }
  *result = sum;
  return GSL_SUCCESS;
}


/* 1F1(1,b,x)
 * b >= 1, b integer
 */
static
int
hyperg_1F1_1_int(const int b, const double x, double * result)
{
  if(b < 1) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(b == 1) {
    return gsl_sf_exp_impl(x, result);
  }
  else if(b == 2) {
    return gsl_sf_exprel_impl(x, result);
  }
  else if(b == 3) {
    return gsl_sf_exprel_2_impl(x, result);
  }
  else {
    return gsl_sf_exprel_n_impl(b-1, x, result);
  }
}


/* 1F1(1,b,x)
 * b >=1, b real
 *
 * checked OK: [GJ] Thu Oct  1 16:46:35 MDT 1998
 */
static
int
hyperg_1F1_1(const double b, const double x, double * result)
{
  double ax = fabs(x);
  double ib = floor(b + 0.1);

  if(b < 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(b == 1.0) {
    return gsl_sf_exp_impl(x, result);
  }
  else if(b >= 1.4*ax) {
    return hyperg_1F1_1_series(b, x, result);
  }
  else if(fabs(b - ib) < locEPS && ib < INT_MAX) {
    return hyperg_1F1_1_int((int)ib, x, result);
  }
  else if(x > 0.0) {
    if(x > 100.0 && b < 0.75*x) {
      double prec;
      return hyperg_1F1_asymp_posx(1.0, b, x, result, &prec);
    }
    else if(b < 1.0e+05) {
      double bp = b + ceil(1.4*x-b) + 1.0;
      double M;
      hyperg_1F1_1_series(bp, x, &M);
      while(bp > b+0.1) {
        /* M(1,b-1) = x/(b-1) M(1,b) + 1 */
        bp -= 1.0;
        M   = 1.0 + x/bp * M;
      }
      *result = M;
      return GSL_SUCCESS;
    }
    else {
      return hyperg_1F1_large2bm4a(1.0, b, x, result);
    }
  }
  else {
    if(ax < 10.0 && b < 10.0) {
      return hyperg_1F1_1_series(b, x, result);
    }
    else if(ax >= 100.0 && locMAX(fabs(2.0-b),1.0) < 0.99*ax) {
      double prec;
      return hyperg_1F1_asymp_negx(1.0, b, x, result, &prec);
    }
    else {
      double prec;
      return hyperg_1F1_luke(1.0, b, x, result, &prec);
    }
  }
}


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
#if 0
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
#endif /* 0 */


/* 1F1'(a,b,x)/1F1(a,b,x)
 * Uses Gautschi's series transformation of the
 * continued fraction. This is apparently the best
 * method for getting this ratio in the stable region.
 * The convergence is monotone and supergeometric
 * when b > x.
 * Assumes a >= -1.
 */
static
int
hyperg_1F1_CF1_p_ser(const double a, const double b, const double x, double * result)
{
  if(a == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
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
}


/* 1F1(a+1,b,x)/1F1(a,b,x)
 *
 * I think this suffers from typical "anomalous convergence".
 * I could not find a region where it was truly useful.
 */
#if 0
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
#endif /* 0 */


/* 1F1(a,b+1,x)/1F1(a,b,x)
 *
 * This seemed to suffer from "anomalous convergence".
 * However, I have no theory for this recurrence.
 */
#if 0
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
#endif /* 0 */


/* 1F1(a,b,x)
 * |a| <= 1, b > 0
 */
static
int
hyperg_1F1_small_a_bgt0(const double a, const double b, const double x, double * result)
{
  double bma = b-a;
  double oma = 1.0-a;
  double ap1mb = 1.0+a-b;
  double abs_bma = fabs(bma);
  double abs_oma = fabs(oma);
  double abs_ap1mb = fabs(ap1mb);

  double ax = fabs(x);

  if(a == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(a == 1.0) {
    return hyperg_1F1_1(b, x, result);
  }
  else if(a == -1.0) {
    *result = 1.0 + a/b * x;
    return GSL_SUCCESS;
  }
  else if(b >= 1.4*ax) {
    double prec;
    return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
  }
  else if(x > 0.0) {
    if(x > 100.0 && abs_bma*abs_oma < 0.9*x) {
      double prec;
      return hyperg_1F1_asymp_posx(a, b, x, result, &prec);
    }
    else if(b < 1.0e+05) {
      double prec;
      double bp = b + ceil(1.4*x-b) + 1.0;
      double Mbp1;
      double Mb;
      double Mbm1;
      gsl_sf_hyperg_1F1_series_impl(a, bp+1.0, x, &Mbp1, &prec);
      gsl_sf_hyperg_1F1_series_impl(a, bp,     x, &Mb,   &prec);
      while(bp > b+0.1) {
        /* Do backward recursion. */
        Mbm1 = ((x+bp-1.0)*Mb - x*(bp-a)/bp*Mbp1)/(bp-1.0);
        bp -= 1.0;
	Mbp1 = Mb;
	Mb   = Mbm1;
      }
      *result = Mb;
      return GSL_SUCCESS;
    }
    else {
      return hyperg_1F1_large2bm4a(a, b, x, result);
    }
  }
  else {
    if(ax < 10.0 && b < 10.0) {
      double prec;
      return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
    }
    else if(ax >= 100.0 && locMAX(abs_ap1mb,1.0) < 0.99*ax) {
      double prec;
      return hyperg_1F1_asymp_negx(a, b, x, result, &prec);
    }
    else {
      double prec;
      return hyperg_1F1_luke(a, b, x, result, &prec);
    }
  }
}


/* 1F1(b+eps,b,x)
 * |eps|<=1, b > 0
 */
static
int
hyperg_1F1_beps_bgt0(const double eps, const double b, const double x, double * result)
{
  double Kummer_1F1;
  int stat_K = hyperg_1F1_small_a_bgt0(-eps, b, -x, &Kummer_1F1);
  if((stat_K == GSL_SUCCESS || stat_K == GSL_ELOSS) && Kummer_1F1 != 0.0) {
    double lnK = log(fabs(Kummer_1F1));
    int stat_e = gsl_sf_exp_sgn_impl(lnK + x, Kummer_1F1, result);
    if(stat_e == GSL_SUCCESS)
      return stat_K;
    else
      return stat_e;
  }
  else {
    *result = 0.0;
    return stat_K;
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
#if 0
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
#endif /* 0 */


/* Determine endpoint of diagonal recursion.
 *
 * given:  M(a,b)    and  M(a+1,b+2)
 * get:    M(a+1,b)  and  M(a+1,b+1)
 */
#if 0
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
#endif /* 0 */


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
    if(K == 0.0 || stat_K != GSL_SUCCESS) {
      *result = 0.0;
      return stat_K;
    }
    else {
      return gsl_sf_exp_sgn_impl(x + log(fabs(K)), K, result);
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
  else if(b > a && b >= 2*a + x) {
    /* Use the Gautschi CF series, then
     * recurse backward to a=0 for normalization.
     * This will work for either sign of x.
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
  else if(b > a && b < 2*a + x && b > x) {
    /* Use the Gautschi series representation of
     * the continued fraction. Then recurse forward
     * to the a=b line for normalization. This will
     * work for either sign of x, although we do need
     * to check for b > x, for when x is positive.
     */
    double rap;
    int stat_CF1 = hyperg_1F1_CF1_p_ser(a, b, x, &rap);
    double ra = 1.0 + x/a * rap;
    double ex;
    int stat_ex;

    if(stat_CF1 == GSL_SUCCESS || stat_CF1 == GSL_EMAXITER) {
      double Ma   = GSL_SQRT_DBL_MIN;
      double Map1 = ra * Ma;
      double Mnm1 = Ma;
      double Mn   = Map1;
      double Mnp1;
      int n;
      for(n=a+1; n<b; n++) {
    	Mnp1 = ((b-n)*Mnm1 + (2*n-b+x)*Mn)/n;
    	Mnm1 = Mn;
    	Mn   = Mnp1;
      }
      stat_ex = gsl_sf_exp_impl(x, &ex);  /* 1F1(b,b,x) */
      if(stat_ex == GSL_SUCCESS) {
    	*result = Ma/Mn * ex;
    	return stat_CF1;
      }
      else {
    	*result= 0.0;
    	return stat_ex;
      }
    }
    else {
      *result = 0.0;
      return stat_CF1;
    }
  }
  else if(x >= 0.0) {

    if(b < a) {
      /* The point b,b is below the b=2a+x line.
       * Forward recursion on a from b,b+1 is possible.
       * Note that a > b + 1 as well, since we already tried a = b + 1.
       */
      if(x + log(fabs(x/b)) < GSL_LOG_DBL_MAX-2.0) {
        double ex = exp(x);
        int n;
        double Mnm1 = ex;		  /* 1F1(b,b,x)   */
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
    else {
      /* b > a
       * b < 2a + x 
       * b <= x (otherwise we would have finished above)
       *
       * Gautschi anomalous convergence region. However, we can
       * recurse forward all the way from a=0,1 because we are
       * always underneath the b=2a+x line.
       */
      double Mnm1 = 1.0;    /* 1F1(0,b,x) */
      double Mn;	    /* 1F1(1,b,x)  */
      double Mnp1;
      int n;
      gsl_sf_exprel_n_impl(b-1, x, &Mn);
      for(n=1; n<a; n++) {
        Mnp1 = ((b-n)*Mnm1 + (2*n-b+x)*Mn)/n;
     	Mnm1 = Mn;
     	Mn   = Mnp1;
      }
      *result = Mn;
      return GSL_SUCCESS;
    }
  }
  else {
    /* x < 0
     * b < a (otherwise we would have tripped one of the above)
     */

    if(a <= 0.5*(b-x) || a >= -x) {
      /* Gautschi continued fraction is in the anomalous region,
       * so we must find another way. We recurse down in b,
       * from the a=b line.
       */
      double Manp1 = exp(x);
      double Man   = exp(x) * (1.0 + x/(a-1.0));
      double Manm1;
      int n;
      for(n=a-1; n>b; n--) {
        Manm1 = (-n*(1-n-x)*Man - x*(n-a)*Manp1)/(n*(n-1.0));
        Manp1 = Man;
        Man = Manm1;
      }
      *result = Man;
      return GSL_SUCCESS;
    }
    else {
      /* Pick a0 such that b ~= 2a0 + x, then
       * recurse down in b from a0,a0 to determine
       * the values near the line b=2a+x. Then recurse
       * forward on a from a0.
       */
      int a0 = ceil(0.5*(b-x));
      double Ma0b;    /* M(a0,b)   */
      double Ma0bp1;  /* M(a0,b+1) */
      double Ma0p1b;  /* M(a0+1,b) */
      double Mnm1;
      double Mn;
      double Mnp1;
      int n;
      {
    	double Ma0np1 = exp(x);
    	double Ma0n   = exp(x) * (1.0 + x/(a0-1.0));
    	double Ma0nm1;
    	for(n=a0-1; n>b; n--) {
    	  Ma0nm1 = (-n*(1-n-x)*Ma0n - x*(n-a0)*Ma0np1)/(n*(n-1.0));
    	  Ma0np1 = Ma0n;
    	  Ma0n = Ma0nm1;
    	}
	Ma0bp1 = Ma0np1;
    	Ma0b   = Ma0n;
	Ma0p1b = (b*(a0+x)*Ma0b + x*(a0-b)*Ma0bp1)/(a0*b);
      }

      Mnm1 = Ma0b;
      Mn   = Ma0p1b;
      for(n=a0+1; n<a; n++) {
    	Mnp1 = ((b-n)*Mnm1 + (2*n-b+x)*Mn)/n;
    	Mnm1 = Mn;
    	Mn   = Mnp1;
      }
      *result = Mn;
      return GSL_SUCCESS;
    }
  }
}


/* Evaluate a <= 0 cases directly. (Polynomial; Horner)
 * When the terms are all positive, this
 * must work. We will assume this here.
 */
static
int
hyperg_1F1_a_negint_poly(const int a, const double b, const double x, double * result)
{
  if(a == 0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else {
    int N = -a;
    double poly = 1.0;
    int k;
    for(k=N-1; k>=0; k--) {
      double t = (a+k)/(b+k) * (x/(k+1));
      double r = t + 1.0/poly;
      if(r > 0.9*DBL_MAX/poly) {
        *result = 0.0; /* FIXME: should be Inf */
	return GSL_EOVRFLW;
      }
      else {
        poly *= r;  /* P_n = 1 + t_n P_{n-1} */
      }
    }
    *result = poly;
    return GSL_SUCCESS;
  }
}


/* Assumes a <= -1,  b <= -1, and b <= a.
 */
static
int
hyperg_1F1_ab_negint(const int a, const int b, const double x, double * result)
{
  if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(x > 0.0) {
    return hyperg_1F1_a_negint_poly(a, b, x, result);
  }
  else {
    /* Apply a Kummer transformation to make x > 0 so
     * we can evaluate the polynomial safely. Of course,
     * this assumes b <= a, which must be true for
     * a<0 and b<0, since otherwise the thing is undefined.
     */
    double K;
    double ex;
    int stat_K = hyperg_1F1_a_negint_poly(b-a, b, -x, &K);
    int stat_e = gsl_sf_exp_impl(x, &ex);
    if(stat_K == GSL_SUCCESS && stat_e == GSL_SUCCESS) {
      *result = ex * K;
      return GSL_SUCCESS; 
    }
    else if(stat_K == GSL_EOVRFLW) {
      *result = 0.0;
      return stat_K;
    }
    else if(stat_e == GSL_EUNDRFLW) {
      *result = 0.0;
      return stat_e;
    }
    else {
      *result = 0.0;
      return GSL_EFAILED;
    }
  }
}


/* Handle case of generic positive a, b.
 */
static
int
hyperg_1F1_ab_pos(const double a, const double b, const double x, double * result)
{
  const double ax = fabs(x);

  if(   ( b < 10.0 && a < 10.0 && ax < 5.0 )
     || ( b > a*ax )
     || ( b > a && ax < 5.0 )
    ) {
    double prec;
    return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
  }
  else if(fabs(b-a) <= 1.0) {
    /* Directly handle b near a.
     */
    return hyperg_1F1_beps_bgt0(a-b, b, x, result);  /* a = b + eps */
  }
  else if(b > a && b >= 2*a + x) {
    /* Use the Gautschi CF series, then
     * recurse backward to a near 0 for normalization.
     * This will work for either sign of x.
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
      double Mn_true;
      int stat_Mt;
      double n;
      for(n=a; n>0.5; n -= 1.0) {
        Mnm1 = (n * Mnp1 - (2.0*n-b+x) * Mn) / (b-n);
        Mnp1 = Mn;
        Mn   = Mnm1;
      }
      stat_Mt = hyperg_1F1_small_a_bgt0(n, b, x, &Mn_true);
      if(stat_Mt == GSL_SUCCESS) {
        *result = (Ma/Mn) * Mn_true;
        return stat_CF1;
      }
      else {
        *result = 0.0;
	return stat_Mt;
      }
    }
    else {
      *result = 0.0;
      return stat_CF1;
    }
  }
  else if(b > a && b < 2*a + x && b > x) {
    /* Use the Gautschi series representation of
     * the continued fraction. Then recurse forward
     * to near the a=b line for normalization. This will
     * work for either sign of x, although we do need
     * to check for b > x, for when x is positive.
     */
    double Mn_true;
    int stat_Mt;
    double rap;
    int stat_CF1 = hyperg_1F1_CF1_p_ser(a, b, x, &rap);
    if(stat_CF1 == GSL_SUCCESS || stat_CF1 == GSL_EMAXITER) {
      double ra = 1.0 + x/a * rap;
      double Ma   = GSL_SQRT_DBL_MIN;
      double Mnm1 = Ma;
      double Mn   = ra * Mnm1;
      double Mnp1;
      double n;
      for(n=a+1.0; n<b-0.5; n += 1.0) {
        Mnp1 = ((b-n)*Mnm1 + (2*n-b+x)*Mn)/n;
        Mnm1 = Mn;
        Mn   = Mnp1;
      }
      stat_Mt = hyperg_1F1_beps_bgt0(n-b, b, x, &Mn_true);
      if(stat_Mt == GSL_SUCCESS || stat_Mt == GSL_ELOSS) {
        *result = Ma/Mn * Mn_true;
        return stat_Mt;
      }
      else {
        *result = 0.0;
        return stat_Mt;
      }
    }
    else {
      *result = 0.0;
      return stat_CF1;
    }
  }
  else if(x >= 0.0) {

    if(b < a) {
      /* Forward recursion on a from a=b+eps-1,b+eps.
       */
      double N   = floor(a-b);
      double eps = a - b - N;
      double M0, M1;
      int stat_0 = hyperg_1F1_beps_bgt0(eps-1.0, b, x, &M0);
      int stat_1 = hyperg_1F1_beps_bgt0(eps,     b, x, &M1);
      if(   (stat_0 == GSL_SUCCESS || stat_0 == GSL_ELOSS)
         && (stat_1 == GSL_SUCCESS || stat_1 == GSL_ELOSS)
	) {
        double Mam1 = M0;
        double Ma   = M1;
        double Map1;
        double ap;
        for(ap=b+eps; ap<a-0.1; ap += 1.0) {
          Map1 = ((b-ap)*Mam1 + (2.0*ap-b+x)*Ma)/ap;
	  Mam1 = Ma;
	  Ma   = Map1;
        }
        *result = Ma;
	if(stat_0 == GSL_ELOSS || stat_1 == GSL_ELOSS)
	  return GSL_ELOSS;
	else
          return GSL_SUCCESS;
      }
      else {
        *result = 0.0;
        return GSL_EFAILED;
      }
    }
    else {
      /* b > a
       * b < 2a + x 
       * b <= x
       *
       * Recurse forward on a from a=eps,eps+1.
       */
      double eps = a - floor(a);
      double Mnm1;
      double Mn;
      double Mnp1;
      int stat_0 = hyperg_1F1_small_a_bgt0(eps,     b, x, &Mnm1);
      int stat_1 = hyperg_1F1_small_a_bgt0(eps+1.0, b, x, &Mn);
      if(stat_0 == GSL_SUCCESS && stat_1 == GSL_SUCCESS) {
        double n;
        for(n=eps+1.0; n<a-0.1; n++) {
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
  }
  else {
    /* x < 0
     * b < a
     */

    if(a <= 0.5*(b-x) || a >= -x) {
      /* Recurse down in b, from near the a=b line, b=a+eps,a+eps-1.
       */
      double N   = floor(a - b);
      double eps = 1.0 + N - a + b;
      double Manp1;
      double Man;
      double Manm1;
      int stat_0 = hyperg_1F1_beps_bgt0(-eps,    a+eps,     x, &Manp1);
      int stat_1 = hyperg_1F1_beps_bgt0(1.0-eps, a+eps-1.0, x, &Man);
      if(   (stat_0 == GSL_SUCCESS || stat_0 == GSL_ELOSS)
         && (stat_1 == GSL_SUCCESS || stat_1 == GSL_ELOSS)
	) {
        double n;
        for(n=a+eps-1.0; n>b+0.1; n -= 1.0) {
          Manm1 = (-n*(1-n-x)*Man - x*(n-a)*Manp1)/(n*(n-1.0));
          Manp1 = Man;
          Man = Manm1;
        }
        *result = Man;
	if(stat_0 == GSL_ELOSS || stat_1 == GSL_ELOSS)
	  return GSL_ELOSS;
	else
          return GSL_SUCCESS;
      }
      else {
        *result = 0.0;
	return GSL_EFAILED;
      }
    }
    else {
      /* Pick a0 such that b ~= 2a0 + x, then
       * recurse down in b from a0,a0 to determine
       * the values near the line b=2a+x. Then recurse
       * forward on a from a0.
       */
      double epsa = a - floor(a);
      double a0   = floor(0.5*(b-x)) + epsa;
      double N    = floor(a0 - b);
      double epsb = 1.0 + N - a0 + b;
      double Ma0b;
      double Ma0bp1;
      double Ma0p1b;
      int stat_a0;
      double Mnm1;
      double Mn;
      double Mnp1;
      double n;
      {
    	double Ma0np1;
    	double Ma0n;
    	double Ma0nm1;
        int stat_0 = hyperg_1F1_beps_bgt0(-epsb,    a0+epsb,     x, &Ma0np1);
        int stat_1 = hyperg_1F1_beps_bgt0(1.0-epsb, a0+epsb-1.0, x, &Ma0n);
	if(   (stat_0 == GSL_SUCCESS || stat_0 == GSL_ELOSS)
	   && (stat_1 == GSL_SUCCESS || stat_1 == GSL_ELOSS)
	  ) {
    	  for(n=a0+epsb-1.0; n>b+0.1; n -= 1.0) {
    	    Ma0nm1 = (-n*(1-n-x)*Ma0n - x*(n-a0)*Ma0np1)/(n*(n-1.0));
            Ma0np1 = Ma0n;
            Ma0n = Ma0nm1;
          }
	  Ma0bp1 = Ma0np1;
          Ma0b   = Ma0n;
	  Ma0p1b = (b*(a0+x)*Ma0b+x*(a0-b)*Ma0bp1)/(a0*b); /* right-down hook */
	  if(stat_0 == GSL_ELOSS || stat_1 == GSL_ELOSS)
	    stat_a0 = GSL_ELOSS;
	  else
	    stat_a0 = GSL_SUCCESS;
	}
	else {
	  *result = 0.0;
	  return GSL_EFAILED;
	}
      }

      Mnm1 = Ma0b;
      Mn   = Ma0p1b;
      for(n=a0+1.0; n<a-0.1; n += 1.0) {
    	Mnp1 = ((b-n)*Mnm1 + (2*n-b+x)*Mn)/n;
    	Mnm1 = Mn;
    	Mn   = Mnp1;
      }
      *result = Mn;
      return stat_a0;
    }
  }



#if 0
  else {
    /* 2a + x >= b >= a - 1
     */
    double ra;
    int stat_CF1;
    int stat_Mt;
    double Mn_true;
    double Ma;
    double Mnm1;
    double Mn;
    double Mnp1;
    double n;

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
       * M(a+1,b,x)/M(a,b,x) = M(b-a-1,b,-x)/M(b-a,b,-x)
       *
       * and the second ratio is in the Gautschi stable
       * region and can be evaluated. [The Kummer transform
       * is a reflection about b=2a].
       */
      double rap_Kummer;
      double ra_Kummer;
      stat_CF1 = hyperg_1F1_CF1_p_ser(b-a-1.0, b, -x, &rap_Kummer);
      ra_Kummer = 1.0 + (-x/(b-a-1.0)) * rap_Kummer;
      ra = 1.0/ra_Kummer;
    }

    /* Recurse forward to near a=b to determine normalization.
     * Since b < 2a + x, this is stable.
     */
    Ma   = GSL_SQRT_DBL_MIN;
    Mnm1 = Ma;
    Mn   = ra * Mnm1;

    for(n=a+1.0; n<b-0.5; n += 1.0) {
      Mnp1 = ((b-n)*Mnm1 + (2*n-b+x)*Mn)/n;
      Mnm1 = Mn;
      Mn   = Mnp1;
    }
    stat_Mt = hyperg_1F1_beps_bgt0(n-b, b, x, &Mn_true);
    if(stat_Mt == GSL_SUCCESS || stat_Mt == GSL_ELOSS) {
      *result = Ma/Mn * Mn_true;
      return stat_Mt;
    }
    else {
      *result = 0.0;
      return stat_Mt;
    }
  }
#endif
  
}


/* Assumes b != integer
 */
static
int
hyperg_1F1_ab_neg(const double a, const double b, const double x,
                  double * result)
{
  const double abs_a = fabs(a);
  const double abs_b = fabs(b);
  const double size_a = locMAX(abs_a, 1.0);
  const double size_b = locMAX(abs_b, 1.0);

  if(   x > 0.0
     && size_b > size_a
     && size_a*log(M_E*x/size_b) < GSL_LOG_MACH_EPS+7.0
    ) {
    /* Series terms are positive definite up until
     * there is a sign change. But by then the
     * terms are small due to the last condition.
     */
    double prec;
    return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
  }
  else if(x > 0.0) {
    /* [Abramowitz+Stegun, 13.1.3]
     *
     * M(1+a-b,2-b,x) = Gamma(a)/Gamma(b) x^(b-1) *
     *                  { Gamma(2-b)/Gamma(1+a-b) M(a,b,x) - (1-b) U(a,b,x) }
     */
    double bp = 2.0 - b;
    double ap = a - b + 1.0;
    double lnpre;
    double lg_ap, lg_bp;
    double sg_ap;
    double lnc1;
    double lg_2mbp, lg_1papmbp;
    double sg_2mbp, sg_1papmbp;
    double M, U;
    double term_M;
    double inner, lninner;
    int stat_F;
    int stat_U;
    int stat_e;

    gsl_sf_lngamma_sgn_impl(ap, &lg_ap, &sg_ap);
    gsl_sf_lngamma_impl(bp, &lg_bp);
    lnpre = lg_ap - lg_bp + (bp-1.0)*log(x);

    gsl_sf_lngamma_sgn_impl(2.0-bp,    &lg_2mbp,    &sg_2mbp);
    gsl_sf_lngamma_sgn_impl(1.0+ap-bp, &lg_1papmbp, &sg_1papmbp);
    lnc1 = lg_2mbp - lg_1papmbp;

    stat_F = gsl_sf_hyperg_1F1_impl(ap, bp, x, &M);
    stat_U = gsl_sf_hyperg_U_impl(ap, bp, x, &U);

    stat_e = gsl_sf_exp_sgn_impl(lnc1+log(fabs(M)), sg_2mbp*sg_1papmbp*M, &term_M);

    inner = term_M - (1.0-bp) * U;
    lninner = log(fabs(inner));
printf("--:  %22.18g   %22.18g\n", term_M, (1.0-bp) * U);

    return gsl_sf_exp_sgn_impl(lninner + lnpre, sg_ap*inner, result);
  }
  else {
    *result = 0.0;
    return GSL_EUNIMPL;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

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
  else if(x > 100.0  && locMAX(1.0,fabs(b-a))*locMAX(1.0,fabs(1-a)) < 0.5 * x) {
    /* x -> +Inf asymptotic
     */
    double prec;
    return hyperg_1F1_asymp_posx(a, b, x, result, &prec);
  }
  else if(x < -100.0 && locMAX(1.0,fabs(a))*locMAX(1.0,fabs(1+a-b)) < 0.5 * fabs(x)) {
    /* x -> -Inf asymptotic
     */
    double prec;
    return hyperg_1F1_asymp_negx(a, b, x, result, &prec);
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
  const double rinta = floor(a + 0.5);
  const double rintb = floor(b + 0.5);
  const double rintbma = floor(bma + 0.5);
  const int a_integer   = ( fabs(a - rinta) < locEPS );
  const int b_integer   = ( fabs(b - rintb) < locEPS );
  const int bma_integer = ( fabs(bma - rintbma) < locEPS );
  const int b_neg_integer   = ( b < -0.1 &&  b_integer );
  const int a_neg_integer   = ( a < -0.1 &&  a_integer );
  const int bma_neg_integer = ( bma < -0.1 &&  bma_integer );

  if(x == 0.0) {
    /* Testing for this before testing a and b
     * is somewhat arbitrary. The result is that
     * we can have 1F1(0,0,0) = 1. Whatever.
     */
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(b == 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(a == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(a == b) {
    /* case: a==b; exp(x)
     * It's good to test exact equality now.
     * We also test approximate equality later.
     */
    return gsl_sf_exp_impl(x, result);
  }
  else if(fabs(b) < 10.0*locEPS) {
    /* Note that neither a nor b is zero, since
     * we eliminated that with the above tests.
     */
    if(fabs(a) < 10.0*locEPS) {
      /* a and b near zero: 1 + a/b (exp(x)-1)
       */
      double exm1;
      int stat_e = gsl_sf_expm1_impl(x, &exm1);
      if(stat_e == GSL_SUCCESS) {
        double sa = ( a > 0.0 ? 1.0 : -1.0 );
        double sb = ( b > 0.0 ? 1.0 : -1.0 );
        double se = ( exm1 > 0.0 ? 1.0 : -1.0 );
        double lnr = log(fabs(a)) - log(fabs(b)) + log(fabs(exm1));
        double hx;
        int stat_hx = gsl_sf_exp_sgn_impl(lnr, sa * sb * se, &hx);
        if(stat_hx == GSL_SUCCESS) {
          *result = (hx == DBL_MAX ? hx : 1.0 + hx);  /* excessive paranoia ? */
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
    else {
      /* b near zero and a not near zero
       */
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
  }
  else if(   (a_integer && a > INT_MIN && a < INT_MAX)
          && (b_integer && b > INT_MIN && b < INT_MAX)
    ) {
    /* Check for reduction to the integer case.
     * Make the arbitrary "near an integer" test.
     */
    return gsl_sf_hyperg_1F1_int_impl((int)rinta, (int)rintb, x, result);
  }
  else if(a_neg_integer && ((b <= a && x > 0.0) || (b > 0.0 && x < 0.0)) && a > INT_MIN) {
    /* The polynomial evaluation is safe since all
     * the terms are positive definite.
     */
    return hyperg_1F1_a_negint_poly((int)rinta, b, x, result);
  }
  else if(bma_neg_integer && ((a < 0.0 && x < 0.0) || (b > 0.0 && x > 0.0)) && bma > INT_MIN) {
    /* Kummer transformed version of safe polynomial.
     */
    double K;
    int stat_K = hyperg_1F1_a_negint_poly((int)rintbma, b, -x, &K);
    if(K != 0.0 && stat_K == GSL_SUCCESS) {
      return gsl_sf_exp_sgn_impl(x + log(fabs(K)), K, result);
    }
    else {
      *result = 0.0;
      return stat_K;
    }
  }
  else if(a_neg_integer && a > INT_MIN && x > 0.0) {
    /* Combine [Abramowitz+Stegun, 13.6.9 + 13.6.27]
     * M(-n,1+y,x) = (-1)^n / (y+1)_n U(-n,1+y,x)
     */
    int n = -(int)rinta;
    double sgn = ( GSL_IS_ODD(n) ? -1.0 : 1.0 );
    double lnpoch;
    double sgpoch;
    double U;
    int stat_p = gsl_sf_lnpoch_sgn_impl(b, n, &lnpoch, &sgpoch);
    int stat_U = gsl_sf_hyperg_U_impl(-n, b, x, &U);
    if(stat_p == GSL_SUCCESS) {
      if(U != 0.0 && (stat_U == GSL_SUCCESS || stat_U == GSL_ELOSS)) {
        double lnr = -lnpoch + log(fabs(U));
        int stat_e = gsl_sf_exp_sgn_impl(lnr, sgn * sgpoch * U, result);
	if(stat_e == GSL_SUCCESS)
	  return stat_U;
	else
	  return stat_e;
      }
      else {
        *result = 0.0;
        return stat_U;
      }
    }
    else {
      *result = 0.0;
      return stat_p;
    }
  }
  else if(b_neg_integer) {
    /* b is a neg integer, but a is not even an integer,
     * so the thing is undefined.
     */
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(   (fabs(x) < 5.0 && fabs(a) < 10.0 && fabs(b) < 10.0)
          || (b > 0.8*locMAX(fabs(a),1.0)*fabs(x))
    ) {
    /* Arguments small enough to evaluate series directly
     * or series is dominated and safe.
     */
    double prec;
    return gsl_sf_hyperg_1F1_series_impl(a, b, x, result, &prec);
  }
  else if(   (fabs(x) < 5.0 && fabs(bma) < 10.0 && fabs(b) < 10.0)
          || (b > 0.8*locMAX(fabs(bma),1.0)*fabs(x))
    ) {
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
  else if(   x < -30.0
          && locMAX(fabs(a),1.0)*locMAX(fabs(1.0+a-b),1.0) < 0.99*fabs(x)
          && !b_neg_integer
          && !bma_neg_integer
    ) {
    /* Large negative x asymptotic.
     */
    double prec;
    return hyperg_1F1_asymp_negx(a, b, x, result, &prec);
  }
  else if(   x > 100.0
          && locMAX(fabs(bma),1.0)*locMAX(fabs(1.0-a),1.0) < 0.99*fabs(x)
          && !b_neg_integer
          && !a_neg_integer
    ) {
    /* Large positive x asymptotic.
     */
    double prec;
    return hyperg_1F1_asymp_posx(a, b, x, result, &prec);
  }
  else if(b > locEPS && fabs(bma) < GSL_SQRT_MACH_EPS && fabs(b) > fabs(x) ) {
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
    double eps = bma;
    double exab;
    int stat_e = gsl_sf_exp_impl(a*x/b, &exab);
    if(stat_e == GSL_SUCCESS) {
      double v2 = a/(2.0*b*b*(b+1.0));
      double v3 = a*(b-2.0*a)/(3.0*b*b*b*(b+1.0)*(b+2.0));
      double v  = v2 + v3 * x;
      *result = exab * (1.0 + eps*x*x*v);
      return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      return stat_e;
    }
  }
  else if(-1.0 <= a && a <= 1.0 && b > 0.0) {
    /* Handle small a explicitly. This helps clarify
     * the thinking for the recursize cases below.
     * It is also a minor premature optimization.
     */
    return hyperg_1F1_small_a_bgt0(a, b, x, result);
  }
  else if(a < 0.0 && b > 0.0) {
    /* Use Kummer to reduce it to the generic positive case.
     * Note that b > a, strictly, since we already trapped b = a.
     */
    double Kummer_1F1;
    int stat_K = hyperg_1F1_ab_pos(b-a, b, -x, &Kummer_1F1);
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
  else if(a > 0.0 && b < 0.0) {
    /* Use Kummer to reduce it to the generic negative case.
     */
    double Kummer_1F1;
    int stat_K = hyperg_1F1_ab_neg(b-a, b, -x, &Kummer_1F1);
    if(stat_K == GSL_SUCCESS) {
      if(Kummer_1F1 == 0.0) {
        *result = 0.0;
	return stat_K;
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
  else if(a < 0.0 && b < 0.0) {
    /* generic negative case
     */
    return hyperg_1F1_ab_neg(a, b, x, result);
  }
  else {
    /* a > 0 and b > 0
     * generic positive case
     */
    return hyperg_1F1_ab_pos(a, b, x, result);
  }
  
  
  
#if 0  
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
#endif

}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

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


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

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
