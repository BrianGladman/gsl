/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "hyperg.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_hyperg.h"

#define locEPS       (1000.0*GSL_MACH_EPS)
#define locMAX(a,b)  ((a) > (b) ? (a) : (b))
#define locMIN(a,b)  ((a) < (b) ? (a) : (b))


/* Evaluate u_{N+1}/u_N by Steed's continued fraction method.
 *
 * u_N := Gamma[a+N]/Gamma[a] U(a + N, b, x)
 */
static
int
hyperg_U_CF1(const double a, const double b, const int N, const double x, double * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 20000;
  int n = 1;
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = -(a + N);
  double b1 =  (b - 2.0*a - x - 2.0*(N+1));
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
    an = -(a + N + n - b)*(a + N + n - 1.0);
    bn =  (b - 2.0*a - x - 2.0*(N+n));
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


/* Large x asymptotic for  x^a U(a,b,x)
 * Based on SLATEC D9CHU() [W. Fullerton]
 *
 * Uses a rational approximation due to Luke.
 * See [Luke, Algorithms for the Computation of Special Functions, p. 252]
 *     [Luke, Utilitas Math. (1977)]
 *
 * z^a U(a,b,z) ~ 2F0(a,1+a-b,-1/z)
 *
 * This assumes that a is not a negative integer and
 * that 1+a-b is not a negative integer. If one of them
 * is, then the 2F0 actually terminates, the above
 * relation is an equality, and the sum should be
 * evaluated directly [see below].
 */
static
int
d9chu(const double a, const double b, const double x, double * result)
{
  const double EPS   = 8.0 * GSL_MACH_EPS;  /* EPS = 4.0D0*D1MACH(4)   */
  const double SQEPS = GSL_MACH_EPS;        /* SQEPS = SQRT(D1MACH(4)) */
  const int maxiter = 500;
  double aa[4], bb[4];
  int i;

  double bp = 1.0 + a - b;
  double ab = a*bp;
  double ct2 = 2.0 * (x - ab);
  double sab = a + bp;
  
  double ct3 = sab + 1.0 + ab;
  double anbn = ct3 + sab + 3.0;
  double ct1 = 1.0 + 2.0*x/anbn;

  bb[0] = 1.0;
  aa[0] = 1.0;

  bb[1] = 1.0 + 2.0*x/ct3;
  aa[1] = 1.0 + ct2/ct3;
  
  bb[2] = 1.0 + 6.0*ct1*x/ct3;
  aa[2] = 1.0 + 6.0*ab/anbn + 3.0*ct1*ct2/ct3;

  for(i=4; i<maxiter; i++) {
    int j;
    double c2;
    double d1z;
    double g1, g2, g3;
    double x2i1 = 2*i - 3;
    ct1   = x2i1/(x2i1-2.0);
    anbn += x2i1 + sab;
    ct2   = (x2i1 - 1.0)/anbn;
    c2    = x2i1*ct2 - 1.0;
    d1z   = 2.0*x2i1*x/anbn;
    
    ct3 = sab*ct2;
    g1  = d1z + ct1*(c2+ct3);
    g2  = d1z - c2;
    g3  = ct1*(1.0 - ct3 - 2.0*ct2);
    
    bb[3] = g1*bb[2] + g2*bb[1] + g3*bb[0];
    aa[3] = g1*aa[2] + g2*aa[1] + g3*aa[0];
    
    if(fabs(aa[3]*bb[0]-aa[0]*bb[3]) < EPS*fabs(bb[3]*bb[0])) break;
    
    for(j=0; j<3; j++) {
      aa[j] = aa[j+1];
      bb[j] = bb[j+1];
    }
  }
  
  *result = aa[3]/bb[3];
  
  if(i == maxiter) {
    return GSL_EMAXITER;
  }
  else if(*result < SQEPS || *result > 1.0/SQEPS) {
    return GSL_ELOSS;
  }
  else {
    return GSL_SUCCESS;
  }
}


/* Evaluate asymptotic for z^a U(a,b,z) ~ 2F0(a,1+a-b,-1/z)
 * We check for termination of the 2F0 as a special case.
 * Assumes x > 0.
 * Also assumes a,b are not too large compared to x.
 */
static
int
hyperg_zaU_asymp(const double a, const double b, const double x, double *result)
{
  const double ap = a;
  const double bp = 1.0 + a - b;
  const int ap_neg_int = ( ap < 0.0 && fabs(ap - rint(ap)) < locEPS );
  const int bp_neg_int = ( bp < 0.0 && fabs(bp - rint(bp)) < locEPS );

  if(ap_neg_int || bp_neg_int) {
    /* Evaluate 2F0 polynomial.
     */
    double mxi = -1.0/x;
    double nmax = -(int)(locMIN(ap,bp) - 0.1);
    double tn  = 1.0;
    double sum = 1.0;
    double n   = 1.0;
    while(n <= nmax) {
      double apn = (ap+n-1.0);
      double bpn = (bp+n-1.0);
      tn  *= ((apn/n)*mxi)*bpn;
      sum += tn;
      n += 1.0;
    }
    *result = sum;
    return GSL_SUCCESS;
  }
  else {
    return d9chu(a,b,x,result);
  }
}


/* Evaluate finite sum which appears below.
 */
static
int
hyperg_U_finite_sum(int N, double a, double b, double x, double xeps,
                    double * result)
{
  int i;
  double sum;

  if(N <= 0) {
    double t = 1.0;
    double poch;
    int stat_poch;

    sum = 1.0;
    for(i=1; i<= -N; i++) {
      double xi1 = i - 1;
      t   *= (a+xi1)*x/((b+xi1)*(xi1+1.0));
      sum += t;
    }

    stat_poch = gsl_sf_poch_impl(1.0+a-b, -a, &poch);

    if(stat_poch == GSL_SUCCESS || stat_poch == GSL_EUNDRFLW) {
      sum *= poch;
      *result = sum;
      return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      return stat_poch;
    }
  }
  else {
    int M = N - 2;
    if(M < 0) {
      *result = 0.0;
      return GSL_SUCCESS;
    }
    else {
      double gbm1;
      double gamr;
      int stat_gbm1;
      int stat_gamr;
      double t = 1.0;

      sum = 1.0;
      for(i=1; i<=M; i++) {
        t   *= (a-b+i)*x/((1.0-b+i)*i);
        sum += t;
      }

      stat_gbm1 = gsl_sf_gamma_impl(b-1.0, &gbm1);
      stat_gamr = gsl_sf_gammainv_impl(a,  &gamr);

      if(stat_gbm1 == GSL_SUCCESS) {
        sum = sum * gbm1 * gamr * pow(x,1.0-N) * xeps;
	*result = sum;
	return GSL_SUCCESS;
      }
      else {
        *result = 0.0;
        return stat_gbm1;
      }
    }
  }
}



/* Based on SLATEC DCHU() [W. Fullerton]
 * Assumes x > 0.
 * This is just a series summation method, and
 * it is not good for large a.
 *
 * I patched up the window for 1+a-b near zero. [GJ]
 */
static
int
hyperg_U_series(const double a, const double b, const double x, double * result)
{
  const double EPS      = 2.0 * GSL_MACH_EPS;  /* EPS = D1MACH(3) */
  const double SQRT_EPS = M_SQRT2 * GSL_SQRT_MACH_EPS;

/*
  if(locMAX(fabs(a),1.0)*locMAX(fabs(1.0+a-b),1.0) < 0.99 * fabs(x)) {
    double asymp;
    int stat_asymp = hyperg_zaU_asymp(a, b, x, &asymp);
    *result = asymp * pow(x, -a);
    return stat_asymp;
  }  
  else 
  */
  if(fabs(1.0 + a - b) < SQRT_EPS) {
    /* ALGORITHM IS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X
     */
    /* We can however do the following:
     * U(a,b,x) = U(a,a+1,x) when 1+a-b=0
     * and U(a,a+1,x) = x^(-a).
     */
    double lnr = -a * log(x);
    return gsl_sf_exp_impl(lnr, result);
  }
  else {
    double aintb = ( b < 0.0 ? ceil(b-0.5) : floor(b+0.5) );
    double beps  = b - aintb;
    int N = aintb;
    
    double lnx  = log(x);
    double xeps = exp(-beps*lnx);

    /* Evaluate finite sum.
     */
    double sum;
    int stat_sum = hyperg_U_finite_sum(N, a, b, x, xeps, &sum);


    /* Evaluate infinite sum. */

    int istrt = ( N < 1 ? 1-N : 0 );
    double xi = istrt;

    double gamr;
    int stat_gamr = gsl_sf_gammainv_impl(1.0+a-b, &gamr);
    double powx   = gsl_sf_pow_int(x, istrt);
    double sarg   = beps*M_PI;
    double sfact  = ( sarg != 0.0 ? sarg/sin(sarg) : 1.0 );
    double factor = sfact * ( GSL_IS_ODD(N) ? -1.0 : 1.0 ) * gamr * powx;
   
    double pochai;
    double gamri1;
    double gamrni;
    int stat_pochai = gsl_sf_poch_impl(a, xi, &pochai);
    int stat_gamri1 = gsl_sf_gammainv_impl(xi + 1.0, &gamri1);
    int stat_gamrni = gsl_sf_gammainv_impl(aintb + xi, &gamrni);

    double pochaxibeps;
    double gamrxi1beps;
    int stat_pochaxibeps = gsl_sf_poch_impl(a, xi-beps, &pochaxibeps);
    int stat_gamrxi1beps = gsl_sf_gammainv_impl(xi + 1.0 - beps, &gamrxi1beps);
    double b0 = factor * pochaxibeps * gamrni * gamrxi1beps;

    if(fabs(xeps-1.0) < 0.5) {
      /*
       C  X**(-BEPS) IS CLOSE TO 1.0D0, SO WE MUST BE
       C  CAREFUL IN EVALUATING THE DIFFERENCES.
       */
      int i;
      double pch1ai;
      double pch1i;
      double poch1bxibeps;
      int stat_pch1ai = gsl_sf_pochrel_impl(a + xi, -beps, &pch1ai);
      int stat_pch1i  = gsl_sf_pochrel_impl(xi + 1.0 - beps, beps, &pch1i);
      int stat_poch1bxibeps = gsl_sf_pochrel_impl(b+xi, -beps, &poch1bxibeps);
      double c0 = factor * pochai * gamrni * gamri1
                  * (-poch1bxibeps + pch1ai - pch1i + beps*pch1ai*pch1i);

      /*
       C  XEPS1 = (1.0 - X**(-BEPS))/BEPS = (X**(-BEPS) - 1.0)/(-BEPS)
       */
      double dexprl;
      int stat_dexprl = gsl_sf_exprel_impl(-beps*lnx, &dexprl);
      double xeps1 = lnx * dexprl;

      double dchu = sum + c0 + xeps1*b0;
      double xn = N;
      
      for(i=1; i<1000; i++) {
        double xi  = istrt + i;
        double xi1 = istrt + i - 1;
	double tmp = (a-1.0)*(xn+2.0*xi-1.0) + xi*(xi-beps);
	double t;
        b0 = (a+xi1-beps)*b0*x/((xn+xi1)*(xi-beps));
        c0 = (a+xi1)*c0*x/((b+xi1)*xi) - tmp * b0 / (xi*(b+xi1)*(a+xi1-beps));
        t = c0 + xeps1*b0;
	dchu += t;
        if (fabs(t) < EPS*fabs(dchu)) break;
      }
     
      *result = dchu;
      if(i == 1000) {
        return GSL_EMAXITER;
      }
      else {
        return GSL_SUCCESS;
      }
    }
    else {
      /*
       C  X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE
       C  STRAIGHTFORWARD FORMULATION IS STABLE.
       */
      int i;
      double dchu;
      double dgamrbxi;
      int stat_dgamrbxi = gsl_sf_gammainv_impl(b+xi, &dgamrbxi);
      double a0 = factor * pochai * dgamrbxi * gamri1 / beps;
      b0 = xeps * b0 / beps;

      dchu = sum + a0 - b0;
      
      for(i=1; i<1000; i++) {
        double xi = istrt + i;
        double xi1 = istrt + i - 1;
	double t;
        a0 = (a+xi1)*a0*x/((b+xi1)*xi);
        b0 = (a+xi1-beps)*b0*x/((aintb+xi1)*(xi-beps));
        t = a0 - b0;
        dchu += t;
        if(fabs(t) < EPS*fabs(dchu)) break;
      }
      
      *result = dchu;
      if(i == 1000) {
        return GSL_EMAXITER;
      }
      else {
        return GSL_SUCCESS;
      }
    }
  }
}


static
int
hyperg_U_small_a(const double a, const double b, const double x, double * result)
{
  if(   (fabs(b) > 400.0 && fabs(x) < 0.95 * fabs(b))
     || (fabs(b) > 50.0  && fabs(x) < 0.70 * fabs(b))
     ) {
    return gsl_sf_hyperg_U_large_b_impl(a, b, x, result);
  }
  else {
    return hyperg_U_series(a, b, x, result);
  }
}


int
gsl_sf_hyperg_U_impl(const double a, const double b, const double x, double * result)
{
  double a = 0.5;
  double b = 3.0;
  double x = 1.0;
  double f;
  hyperg_U_CF1(a, b, 10, x, &f);
  printf("%22.18g\n", f);
  exit(0);
  
  if(x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }

  if(locMAX(fabs(a),1.0)*locMAX(fabs(1.0+a-b),1.0) < fabs(x)) {
    double asymp;
    int stat_asymp = hyperg_zaU_asymp(a, b, x, &asymp);
    *result = asymp * pow(x, -a);
    return stat_asymp;
  }


  return hyperg_U_small_a(a, b, x, result);

#if 0
  if(b <= 0.0) {
    /* Use the reflection formula
     * U(a,b,x) = x^(1-b) U(1+a-b,2-b,x)
     */
    double ap = 1.0 + a - b;
    double bp = 2.0 - b;
    double powx;
    double U;
    int stat_powx = gsl_sf_exp_impl((1.0-b)*log(x), &powx);
    int stat_U = hyperg_U(ap, bp, x, &U);
    if(stat_powx == GSL_EUNDRFLW || stat_U == GSL_EUNDRFLW) {
      *result = 0.0;
      return GSL_EUNDRFLW;
    }
    else if(stat_powx == GSL_SUCCESS && stat_U == GSL_SUCCESS) {
      *result = powx * U;
      return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      return GSL_EFAILED;
    }
  }
#endif
}


int
gsl_sf_hyperg_U_e(const double a, const double b, const double x, double * result)
{
  int status = gsl_sf_hyperg_U_impl(a, b, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("  gsl_sf_hyperg_U_e", status);
  }
  return status;
}


double
gsl_sf_hyperg_U(const double a, const double b, const double x)
{
  double y;
  int status = gsl_sf_hyperg_U_impl(a, b, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_hyperg_U", status);
  }
  return y;
}
