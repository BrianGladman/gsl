/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "hyperg.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_bessel.h"
#include "gsl_sf_laguerre.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_hyperg.h"

#define locEPS       (1000.0*GSL_MACH_EPS)


/* Log[U(a,2a,x)]
 * [Abramowitz+stegun, 13.6.21]
 * Assumes x > 0, a > 1/2.
 */
static
int
hyperg_lnU_beq2a(const double a, const double x, double * result)
{
  const double nu = a - 0.5;
  const double lnpre = 0.5*(x - M_LNPI) - nu*log(x);
  double lnK;
  gsl_sf_bessel_lnKnu_impl(nu, 0.5*x, &lnK);
  *result = lnpre + lnK;
  return GSL_SUCCESS;
}


/* Evaluate u_{N+1}/u_N by Steed's continued fraction method.
 *
 * u_N := Gamma[a+N]/Gamma[a] U(a + N, b, x)
 *
 * u_{N+1}/u_N = (a+N) U(a+N+1,b,x)/U(a+N,b,x)
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
  const double rintap = floor(ap + 0.5);
  const double rintbp = floor(bp + 0.5);
  const int ap_neg_int = ( ap < 0.0 && fabs(ap - rintap) < locEPS );
  const int bp_neg_int = ( bp < 0.0 && fabs(bp - rintbp) < locEPS );

  if(ap_neg_int || bp_neg_int) {
    /* Evaluate 2F0 polynomial.
     */
    double mxi = -1.0/x;
    double nmax = -(int)(GSL_MIN(ap,bp) - 0.1);
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

    int stat_gamall = GSL_ERROR_SELECT_5(stat_sum, stat_gamr, stat_pochai, stat_gamri1, stat_gamrni);

    double pochaxibeps;
    double gamrxi1beps;
    int stat_pochaxibeps = gsl_sf_poch_impl(a, xi-beps, &pochaxibeps);
    int stat_gamrxi1beps = gsl_sf_gammainv_impl(xi + 1.0 - beps, &gamrxi1beps);

    int stat_all = GSL_ERROR_SELECT_3(stat_gamall, stat_pochaxibeps, stat_gamrxi1beps);

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

      stat_all = GSL_ERROR_SELECT_5(stat_all, stat_dexprl, stat_poch1bxibeps, stat_pch1i, stat_pch1ai);

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
        return stat_all;
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

      stat_all = GSL_ERROR_SELECT_2(stat_all, stat_dgamrbxi);

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
        return stat_all;
      }
    }
  }
}


/* Assumes b > 0 and x > 0.
 */
static
int
hyperg_U_small_ab(const double a, const double b, const double x, double * result)
{
  if(a == -1.0) {
    /* U(-1,c+1,x) = Laguerre[c,0,x] = -b + x
     */
    *result = -b + x;
    return GSL_SUCCESS;
  }
  else if(a == 0.0) {
    /* U(0,c+1,x) = Laguerre[c,0,x] = 1
     */
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(GSL_MAX_DBL(fabs(a),1.0)*GSL_MAX_DBL(fabs(1.0+a-b),1.0) < 0.99*fabs(x)) {
    double asymp;
    int stat_asymp = hyperg_zaU_asymp(a, b, x, &asymp);
    *result = asymp * pow(x, -a);
    return stat_asymp;
  }
  else {
    return hyperg_U_series(a, b, x, result);
  }
}


/* Assumes b > 0 and x > 0.
 */
static
int
hyperg_U_small_a_bgt0(const double a, const double b, const double x,
                      double * result,
		      double * ln_multiplier
		      )
{
  if(a == 0.0) {
    *result = 1.0;
    *ln_multiplier = 0.0;
    return GSL_SUCCESS;
  }
  else if(   (b > 5000.0 && x < 0.90 * fabs(b))
          || (b >  500.0 && x < 0.50 * fabs(b))
    ) {
    int stat = gsl_sf_hyperg_U_large_b_impl(a, b, x, result, ln_multiplier);
    if(stat == GSL_EOVRFLW)
      return GSL_SUCCESS;
    else
      return stat;
  }
  else if(b > 15.0) {
    /* Recurse up from b near 1.
     */
    double eps = b - floor(b);
    double b0  = 1.0 + eps;
    double Ubm1;
    double Ub;
    double Ubp1;
    double bp;
    int stat_0 = hyperg_U_small_ab(a, b0,     x, &Ubm1);
    int stat_1 = hyperg_U_small_ab(a, b0+1.0, x, &Ub);
    if(   (stat_0 == GSL_SUCCESS || stat_0 == GSL_ELOSS)
       && (stat_1 == GSL_SUCCESS || stat_1 == GSL_ELOSS)
      ) {
      for(bp = b0+1.0; bp<b-0.1; bp += 1.0) {
        Ubp1 = ((1.0+a-bp)*Ubm1 + (bp+x-1.0)*Ub)/x;
        Ubm1 = Ub;
        Ub   = Ubp1;
      }
      *result = Ub;
      *ln_multiplier = 0.0;
      if(stat_0 == GSL_ELOSS || stat_1 == GSL_ELOSS)
        return GSL_ELOSS;
      else
        return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      *ln_multiplier = 0.0;
      return GSL_EFAILED;
    }
  }
  else {
    *ln_multiplier = 0.0;
    return hyperg_U_small_ab(a, b, x, result);
  }
}


/* We use this to keep track of large
 * dynamic ranges in the recursions.
 * This can be important because sometimes
 * we want to calculate a very large and
 * a very small number and the answer is
 * the product, of order 1. This happens,
 * for instance, when we apply a Kummer
 * transform to make b positive and
 * both x and b are large.
 */
#define RESCALE_2(u0,u1,factor,count)      \
do {                                       \
  double au0 = fabs(u0);                   \
  if(au0 > factor) {                       \
    u0 /= factor;                          \
    u1 /= factor;                          \
    count++;                               \
  }                                        \
  else if(au0 < 1.0/factor) {              \
    u0 *= factor;                          \
    u1 *= factor;                          \
    count--;                               \
  }                                        \
} while (0)


/* Specialization to b >= 1, for integer parameters.
 * Assumes x > 0.
 */
static
int
hyperg_U_int_bge1(const int a, const int b, const double x,
                  double * result,
	          double * ln_multiplier)
{
  if(a == 0) {
    *result = 1.0;
    *ln_multiplier = 0.0;
    return GSL_SUCCESS;
  }
  else if(a == -1) {
    *result = -b + x;
    *ln_multiplier = 0.0;
    return GSL_SUCCESS;
  }
  else if(b == a + 1) {
    /* U(a,a+1,x) = x^(-a)
     */
    return gsl_sf_exp_impl(-a*log(x), result);
  }
  else if(fabs(a)*GSL_MAX_DBL(fabs(1+a-b),1.0) < 0.99*fabs(x)) {
    double asymp;
    int stat_asymp = hyperg_zaU_asymp(a, b, x, &asymp);
    *result = asymp;
    *ln_multiplier = -a*log(x);  /* log(pow(x, -a)) */
    return stat_asymp;
  }
  else if(abs(a) < 5 && b < 5 && x < 2.0) {
    *ln_multiplier = 0.0;
    return hyperg_U_series(a, b, x, result);
  }
  else if(a < 0) {
    /* Recurse backward from a = -1,0.
     */
    const double scale_factor = GSL_SQRT_DBL_MAX;
    int scale_count = 0;
    double Uap1 = 1.0;     /* U(0,b,x)  */
    double Ua   = -b + x;  /* U(-1,b,x) */
    double Uam1;
    int ap;

    for(ap=-1; ap>a; ap--) {
      Uam1 = ap*(b-ap-1.0)*Uap1 + (x+2.0*ap-b)*Ua;
      Uap1 = Ua;
      Ua   = Uam1;
      RESCALE_2(Ua,Uap1,scale_factor,scale_count);
    }
    *result = Ua;
    *ln_multiplier = ( scale_count != 0 ? scale_count*log(scale_factor) : 0.0 );
    return GSL_SUCCESS;
  }
  else if(b >= 2.0*a + x) {
    /* Recurse forward from a = 0,1.
     */
    const double scale_factor = GSL_SQRT_DBL_MAX;
    int scale_count = 0;
    double Uam1 = 1.0;  /* U(0,b,x) */
    double Ua;
    double Uap1;
    int ap;
    double lm;
    int stat_1 = hyperg_U_small_a_bgt0(1.0, b, x, &Ua, &lm);  /* U(1,b,x) */

    if(stat_1 == GSL_SUCCESS || stat_1 == GSL_ELOSS) {

      Uam1 *= exp(-lm);

      for(ap=1; ap<a; ap++) {
        Uap1 = -(Uam1 + (b-2.0*ap-x)*Ua)/(ap*(1.0+ap-b));
        Uam1 = Ua;
        Ua   = Uap1;
        RESCALE_2(Ua,Uam1,scale_factor,scale_count);
      }

      *result = Ua;
      *ln_multiplier = lm + scale_count * log(scale_factor);
      return stat_1;
    }
    else {
      *result = 0.0;
      *ln_multiplier = 0.0;
      return stat_1;
    }
  }
  else {
    if(b <= x) {
      /* Recurse backward either to the b=a+1 line
       * or to a=0, whichever we hit.
       */
      const double scale_factor = GSL_SQRT_DBL_MAX;
      int scale_count = 0;
      int stat_CF1;
      double ru;
      int a_target;
      double lnU_target;

      if(b < a + 1) {
        a_target = b-1;
	lnU_target = -a_target*log(x);
      }
      else {
        a_target = 0;
	lnU_target = 0.0;
      }

      stat_CF1 = hyperg_U_CF1(a, b, 0, x, &ru);

      if(stat_CF1 == GSL_SUCCESS || stat_CF1 == GSL_EMAXITER) {
        double Ua   = 1.0;
        double Uap1 = ru/a * Ua;
        double Uam1;
        int ap;
        for(ap=a; ap>a_target; ap--) {
          Uam1 = -((b-2.0*ap-x)*Ua + ap*(1.0+ap-b)*Uap1);
          Uap1 = Ua;
          Ua   = Uam1;
	  RESCALE_2(Ua,Uap1,scale_factor,scale_count);
        }

        if(Ua == 0.0) {
	  *result = 0.0;
	  *ln_multiplier = 0.0;
	  return GSL_EZERODIV;
	}
	else {
	  double lnscl = -scale_count*log(scale_factor);
	  double lnpre = lnU_target + lnscl;
	  int stat_e = gsl_sf_exp_mult_impl(lnpre, 1.0/Ua, result);
	  if(stat_e == GSL_SUCCESS) {
	    *ln_multiplier = 0.0;
	  }
	  else {
	    *result = 1.0/Ua;
	    *ln_multiplier = lnpre;
	  }
	  return stat_CF1;
	}
      }
      else {
        *result = 0.0;
	*ln_multiplier = 0.0;
	return stat_CF1;
      }
    }
    else {
      /* Recurse backward to near the b=2a+x line, then
       * determine normalization by either direct evaluation
       * or by a forward recursion. The direct evaluation
       * is needed when x is small (which is precisely
       * when it is easy to do).
       */
      const double scale_factor = GSL_SQRT_DBL_MAX;
      int scale_count_for = 0;
      int scale_count_bck = 0;
      int a0 = 1;
      int a1 = a0 + ceil(0.5*(b-x) - a0);
      double Ua1_bck;
      double Ua1_for;
      int stat_for;
      int stat_bck;
      double lm_for;

      {
        /* Recurse back to determine U(a1,b), sans normalization.
         */
        double ru;
        int stat_CF1 = hyperg_U_CF1(a, b, 0, x, &ru);
        double Ua   = 1.0;
        double Uap1 = ru/a * Ua;
        double Uam1;
        int ap;
        for(ap=a; ap>a1; ap--) {
          Uam1 = -((b-2.0*ap-x)*Ua + ap*(1.0+ap-b)*Uap1);
          Uap1 = Ua;
          Ua   = Uam1;
	  RESCALE_2(Ua,Uap1,scale_factor,scale_count_bck);
        }
        Ua1_bck  = Ua;
        stat_bck = stat_CF1;
      }

      if(b == 2*a1 && a1 > 1) {
        /* This can happen when x is small, which is
	 * precisely when we need to be careful with
	 * this evaluation.
	 */
	hyperg_lnU_beq2a((double)a1, x, &lm_for);
	Ua1_for = 1.0;
        stat_for = GSL_SUCCESS;
      }
      else if(b == 2*a1 - 1 && a1 > 1) {
        /* Similar to the above. Happens when x is small.
	 * Use
	 *   U(a,2a-1) = (x U(a,2a) - U(a-1,2(a-1))) / (2a - 2)
	 */
	double lnU00, lnU12;
	double U00, U12;
	hyperg_lnU_beq2a(a1-1.0, x, &lnU00);
	hyperg_lnU_beq2a(a1,	 x, &lnU12);
	lm_for = GSL_MAX(lnU00, lnU12);
	U00 = exp(lnU00 - lm_for);
	U12 = exp(lnU12 - lm_for);
	Ua1_for = (x * U12 - U00) /(2.0*a1 - 2.0);
	stat_for = GSL_SUCCESS;
      }
      else {
        /* Recurse forward to determine U(a1,b) with
         * absolute normalization.
         */
        double Uam1 = 1.0;  /* U(a0-1,b,x) = U(0,b,x) */
        double Ua;
        double Uap1;
        int ap;
        stat_for = hyperg_U_small_a_bgt0(a0, b, x, &Ua, &lm_for); /* U(1,b,x) */

        Uam1 *= exp(-lm_for);

        for(ap=a0; ap<a1; ap++) {
          Uap1 = -(Uam1 + (b-2.0*ap-x)*Ua)/(ap*(1.0+ap-b));
          Uam1 = Ua;
          Ua   = Uap1;
	  RESCALE_2(Ua,Uam1,scale_factor,scale_count_for);
        }
        Ua1_for = Ua;
      }

      /* Now do the matching to produce the final result.
       */
      if(Ua1_bck == 0.0) {
        *result = 0.0;
        *ln_multiplier = 0.0;
        return GSL_EZERODIV;
      }
      else if(Ua1_for == 0.0) {
        /* Should never happen. */
        *result = 0.0;
	*ln_multiplier = 0.0;
	return GSL_EUNDRFLW;
      }
      else {
        double lns = (scale_count_for - scale_count_bck)*log(scale_factor);
        double lnr = lm_for + log(fabs(Ua1_for)) - log(fabs(Ua1_bck)) + lns;
	double sgn = GSL_SIGN(Ua1_for) * GSL_SIGN(Ua1_bck);
        int stat_e = gsl_sf_exp_sgn_impl(lnr, sgn, result);
        if(stat_e == GSL_SUCCESS) {
          *ln_multiplier = 0.0;
        }
        else {
          *result = sgn;
          *ln_multiplier = lnr;
        }
	if(stat_bck == GSL_EFAILED || stat_for == GSL_EFAILED)
          return GSL_EFAILED;
        else if(stat_bck == GSL_ELOSS || stat_for == GSL_ELOSS)
          return GSL_ELOSS;
        else 
          return GSL_SUCCESS;
      }
    }
  }
}


/* Handle b >= 1 for generic a,b values.
 */
static
int
hyperg_U_bge1(const double a, const double b, const double x,
              double * result,
	      double * ln_multiplier)
{
  const double rinta = floor(a+0.5);
  const int a_neg_integer = (a < 0.0 && fabs(a - rinta) < locEPS);

  *ln_multiplier = 0.0;

  if(a == 0.0) {
    *ln_multiplier = 0.0;
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(a_neg_integer && fabs(rinta) < INT_MAX) {
    /* U(-n,b,x) = (-1)^n n! Laguerre[n,b-1,x]
     */
    const int n = -(int)rinta;
    const double sgn = (GSL_IS_ODD(n) ? -1.0 : 1.0);
    double lnfact;
    double L;
    const int stat_L = gsl_sf_laguerre_n_impl(n, b-1.0, x, &L);
    gsl_sf_lnfact_impl(n, &lnfact);
    if(L != 0.0 || (stat_L == GSL_SUCCESS || stat_L == GSL_ELOSS)) {
      const int stat_e = gsl_sf_exp_mult_impl(lnfact, sgn*L, result);
      if(stat_e == GSL_SUCCESS) {
        *ln_multiplier = 0.0;
      }
      else {
        *ln_multiplier = lnfact;
	*result = sgn*L;
      }
      return stat_L;
    }
    else {
      *ln_multiplier = 0.0;
      *result = 0.0;
      return stat_L;
    }
  }
  else if(GSL_MAX_DBL(fabs(a),1.0)*GSL_MAX_DBL(fabs(1.0+a-b),1.0) < 0.99*fabs(x)) {
    double asymp;
    int stat_asymp = hyperg_zaU_asymp(a, b, x, &asymp);
    *result = asymp;
    *ln_multiplier = -a*log(x);  /* log(pow(x, -a)) */
    return stat_asymp;
  }
  else if(fabs(a) <= 1.0) {
    return hyperg_U_small_a_bgt0(a, b, x, result, ln_multiplier);
  }
  else if(fabs(a) < 5.0 && b < 5.0 && x < 2.0) {
    *ln_multiplier = 0.0;
    return hyperg_U_series(a, b, x, result);
  }
  else if(a < 0.0) {
    /* Recurse backward on a and then upward on b.
     */
    const double scale_factor = GSL_SQRT_DBL_MAX;
    const double a0 = a - floor(a) - 1.0;
    const double b0 = b - floor(b) + 1.0;
    int scale_count = 0;
    double lm_0, lm_1, lm_max;
    double Uap1;
    double Ua;
    double Uam1;
    double ap;
    int stat_0 = hyperg_U_small_a_bgt0(a0+1.0, b0, x, &Uap1, &lm_0);
    int stat_1 = hyperg_U_small_a_bgt0(a0,     b0, x, &Ua,   &lm_1);
    lm_max = GSL_MAX(lm_0, lm_1);
    Uap1 *= exp(lm_0-lm_max);
    Ua   *= exp(lm_1-lm_max);

    if(   (stat_0 == GSL_SUCCESS || stat_0 == GSL_ELOSS)
       && (stat_1 == GSL_SUCCESS || stat_1 == GSL_ELOSS)
      ) {

      /* Downward recursion on a.
       */
      for(ap=a0; ap>a+0.1; ap -= 1.0) {
        Uam1 = ap*(b0-ap-1.0)*Uap1 + (x+2.0*ap-b0)*Ua;
	Uap1 = Ua;
	Ua   = Uam1;
	RESCALE_2(Ua,Uap1,scale_factor,scale_count);
      }

      if(b < 2.0) {
        /* b == b0, so no recursion necessary
	 */
	*result = Ua;
	*ln_multiplier = lm_max + scale_count * log(scale_factor);
      }
      else {
        /* Upward recursion on b.
         */
        double Ubm1 = Ua;                                 /* U(a,b0)   */
        double Ub   = (a*(b0-a-1.0)*Uap1 + (a+x)*Ua)/x;   /* U(a,b0+1) */
        double Ubp1;
        double bp;
        for(bp=b0+1.0; bp<b-0.1; bp += 1.0) {
          Ubp1 = ((1.0+a-bp)*Ubm1 + (bp+x-1.0)*Ub)/x;
          Ubm1 = Ub;
          Ub   = Ubp1;
	  RESCALE_2(Ub,Ubm1,scale_factor,scale_count);
        }
        *result = Ub;
        *ln_multiplier = lm_max + scale_count * log(scale_factor);
      }
      if(stat_0 == GSL_ELOSS || stat_1 == GSL_ELOSS)
        return GSL_ELOSS;
      else
        return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      *ln_multiplier = 0.0;
      return GSL_EFAILED;
    }
  }
  else if(b >= 2*a + x) {
    /* Recurse forward from a near zero.
     * Note that we cannot cross the singularity at
     * the line b=a+1, because the only way we could
     * be in that little wedge is if a < 1. But we
     * have already dealt with the small a case.
     */
    const double scale_factor = GSL_SQRT_DBL_MAX;
    int scale_count = 0;
    double lm_0, lm_1, lm_max;
    double a0 = a - floor(a);
    double Uam1;
    double Ua;
    double Uap1;
    double ap;
    int stat_0 = hyperg_U_small_a_bgt0(a0-1.0, b, x, &Uam1, &lm_0);
    int stat_1 = hyperg_U_small_a_bgt0(a0,     b, x, &Ua,   &lm_1);
    lm_max = GSL_MAX(lm_0, lm_1);
    Uam1 *= exp(lm_0-lm_max);
    Ua   *= exp(lm_1-lm_max);
    
    if(   (stat_0 == GSL_SUCCESS || stat_0 == GSL_ELOSS)
       && (stat_1 == GSL_SUCCESS || stat_1 == GSL_ELOSS)
      ) {
      for(ap=a0; ap<a-0.1; ap += 1.0) {
        Uap1 = -(Uam1 + (b-2.0*ap-x)*Ua)/(ap*(1.0+ap-b));
        Uam1 = Ua;
        Ua   = Uap1;
        RESCALE_2(Ua,Uam1,scale_factor,scale_count);
      }
      *result = Ua;
      *ln_multiplier = lm_max + scale_count * log(scale_factor);
      if(stat_0 == GSL_ELOSS || stat_1 == GSL_ELOSS)
        return GSL_ELOSS;
      else
        return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      *ln_multiplier = 0.0;
      return GSL_EFAILED;
    }
  }
  else {
    if(b <= x) {
      /* Recurse backward to a near zero.
       */
      const double scale_factor = GSL_SQRT_DBL_MAX;
      int scale_count = 0;
      double lm_0;
      double a0 = a - floor(a);
      double Uap1;
      double Ua;
      double Uam1;
      double U0;
      int stat_U0;
      double ap;
      double ru;
      double r;
      int stat_CF1 = hyperg_U_CF1(a, b, 0, x, &ru);
      r = ru/a;
      Ua   = GSL_SQRT_DBL_MIN;
      Uap1 = r * Ua;
      for(ap=a; ap>a0+0.1; ap -= 1.0) {
        Uam1 = -((b-2.0*ap-x)*Ua + ap*(1.0+ap-b)*Uap1);
        Uap1 = Ua;
        Ua   = Uam1;
	RESCALE_2(Ua,Uap1,scale_factor,scale_count);
      }
      
      stat_U0 = hyperg_U_small_a_bgt0(a0, b, x, &U0, &lm_0);

      if(stat_U0 == GSL_SUCCESS) {
        *result = GSL_SQRT_DBL_MIN*(U0/Ua);
	*ln_multiplier = lm_0 - scale_count * log(scale_factor);
	return stat_CF1;
      }
      else if(stat_U0 == GSL_ELOSS) {
        *result = GSL_SQRT_DBL_MIN*(U0/Ua);
	*ln_multiplier = lm_0 - scale_count * log(scale_factor);
	return GSL_ELOSS;
      }
      else {
        *result = 0.0;
	*ln_multiplier = 0.0;
	return stat_U0;
      }
    }
    else {
      /* Recurse backward to near the b=2a+x line, then
       * forward from a near zero to get the normalization.
       */
      const double scale_factor = GSL_SQRT_DBL_MAX;
      int scale_count_for = 0;
      int scale_count_bck = 0;
      double eps = a - floor(a);
      double a0 = ( eps == 0.0 ? 1.0 : eps );
      double a1 = a0 + ceil(0.5*(b-x) - a0);
      double Ua1_bck;
      double Ua1_for;
      int stat_for;
      int stat_bck;
      double lm_for;

      {
        /* Recurse back to determine U(a1,b), sans normalization.
         */
        double Uap1;
        double Ua;
        double Uam1;
        double ap;
        double ru;
        double r;
        int stat_CF1 = hyperg_U_CF1(a, b, 0, x, &ru);
        r = ru/a;
        Ua   = GSL_SQRT_DBL_MIN;
        Uap1 = r * Ua;
        for(ap=a; ap>a1+0.1; ap -= 1.0) {
          Uam1 = -((b-2.0*ap-x)*Ua + ap*(1.0+ap-b)*Uap1);
          Uap1 = Ua;
          Ua   = Uam1;
	  RESCALE_2(Ua,Uap1,scale_factor,scale_count_bck);
        }
        Ua1_bck = Ua;
        stat_bck = stat_CF1;
      }
      {
        /* Recurse forward to determine U(a1,b) with
         * absolute normalization.
         */
	double lm_0, lm_1;
        double Uam1;
        double Ua;
        double Uap1;
        double ap;
        int stat_0 = hyperg_U_small_a_bgt0(a0-1.0, b, x, &Uam1, &lm_0);
        int stat_1 = hyperg_U_small_a_bgt0(a0,     b, x, &Ua,   &lm_1);
	lm_for = GSL_MAX(lm_0, lm_1);
	Uam1 *= exp(lm_0 - lm_for);
	Ua   *= exp(lm_1 - lm_for);

        if(   (stat_0 == GSL_SUCCESS || stat_0 == GSL_ELOSS)
           && (stat_1 == GSL_SUCCESS || stat_1 == GSL_ELOSS)
          ) {
          for(ap=a0; ap<a1-0.1; ap += 1.0) {
            Uap1 = -(Uam1 + (b-2.0*ap-x)*Ua)/(ap*(1.0+ap-b));
            Uam1 = Ua;
            Ua   = Uap1;
	    RESCALE_2(Ua,Uam1,scale_factor,scale_count_for);
          }
          Ua1_for = Ua;
          if(stat_0 == GSL_ELOSS || stat_1 == GSL_ELOSS)
            stat_for = GSL_ELOSS;
          else 
            stat_for = GSL_SUCCESS;
        }
        else {
          *result = 0.0;
	  *ln_multiplier = 0.0;
	  return GSL_EFAILED;
        }
      }

      *result = GSL_SQRT_DBL_MIN*Ua1_for/Ua1_bck;
      *ln_multiplier = lm_for + (scale_count_for - scale_count_bck)*log(scale_factor);

      if(stat_bck == GSL_EFAILED || stat_for == GSL_EFAILED)
        return GSL_EFAILED;
      else if(stat_bck == GSL_ELOSS || stat_for == GSL_ELOSS)
        return GSL_ELOSS;
      else 
        return GSL_SUCCESS;
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/


int
gsl_sf_hyperg_U_int_impl(const int a, const int b, const double x, double * result)
{
  if(x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    double ln_multiplier;
    double ln_pre;
    double U;
    int stat_U;
    if(b >= 1) {
      stat_U = hyperg_U_int_bge1(a, b, x, &U, &ln_multiplier);
      ln_pre = 0.0;
    }
    else {
      /* Use the reflection formula
       * U(a,b,x) = x^(1-b) U(1+a-b,2-b,x)
       */
      int ap = 1 + a - b;
      int bp = 2 - b;
      stat_U = hyperg_U_int_bge1(ap, bp, x, &U, &ln_multiplier);
      ln_pre = (1.0-b)*log(x);
    }

    if(U != 0.0 && (stat_U == GSL_SUCCESS || stat_U == GSL_ELOSS)) {
      double ln_r = ln_pre + ln_multiplier;
      int stat_e = gsl_sf_exp_mult_impl(ln_r, U, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_U);
    }
    else {
      *result = 0.0;
      return stat_U;
    }
  }
}


int
gsl_sf_hyperg_U_impl(const double a, const double b, const double x, double * result)
{
  const double rinta = floor(a + 0.5);
  const double rintb = floor(b + 0.5);
  const int a_integer = ( fabs(a - rinta) < locEPS );
  const int b_integer = ( fabs(b - rintb) < locEPS );

  if(x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(a == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(a_integer && b_integer) {
    return gsl_sf_hyperg_U_int_impl(rinta, rintb, x, result);
  }
  else {
    double ln_multiplier;
    double ln_pre;
    double U;
    int stat_U;
    if(b >= 1.0) {
      /* Use b >= 1 function.
       */
      stat_U = hyperg_U_bge1(a, b, x, &U, &ln_multiplier);
      ln_pre = 0.0;
    }
    else {
      /* Use the reflection formula
       * U(a,b,x) = x^(1-b) U(1+a-b,2-b,x)
       */
      double ap = 1.0 + a - b;
      double bp = 2.0 - b;
      stat_U = hyperg_U_bge1(ap, bp, x, &U, &ln_multiplier);
      ln_pre = (1.0-b)*log(x);
    }

    if(U != 0.0 && (stat_U == GSL_SUCCESS || stat_U == GSL_ELOSS)) {
      double ln_r = ln_pre + ln_multiplier;
      int stat_e = gsl_sf_exp_mult_impl(ln_r, U, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_U);
    }
    else {
      *result = 0.0;
      return stat_U;
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/


int
gsl_sf_hyperg_U_int_e(const int a, const int b, const double x, double * result)
{
  int status = gsl_sf_hyperg_U_int_impl(a, b, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("  gsl_sf_hyperg_U_int_e", status);
  }
  return status;
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



/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/


double
gsl_sf_hyperg_U_int(const int a, const int b, const double x)
{
  double y;
  int status = gsl_sf_hyperg_U_int_impl(a, b, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_hyperg_U_int", status);
  }
  return y;
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
