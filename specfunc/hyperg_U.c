/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_hyperg.h"

#define locMAX(a,b)  ((a) > (b) ? (a) : (b))


/* Large x asymptotic for  x^a U(a,b,x)
 * Based on SLATEC D9CHU() [W. Fullerton]
 *
 * Uses a rational approximation due to Luke.
 */
static
int
d9chu(const double a, const double b, const double x, double * result)
{
  const double EPS   = 8.0 * GSL_MACH_EPS;  /* EPS = 4.0D0*D1MACH(4)   */
  const double SQEPS = GSL_MACH_EPS;        /* SQEPS = SQRT(D1MACH(4)) */
  const int maxiter = 500;
  double aa[5], bb[5];
  int i;

  double bp = 1.0 - a + b;
  double ab = a*bp;
  double ct2 = 2.0 * (x - ab);
  double sab = a + bp;
  
  double ct3 = sab + 1.0 + ab;
  double anbn = ct3 + sab + 3.0;
  double ct1 = 1.0 + 2.0*x/anbn;

  bb[1] = 1.0;
  aa[1] = 1.0;

  bb[2] = 1.0 + 2.0*x/ct3;
  aa[2] = 1.0 + ct2/ct3;
  
  bb[3] = 1.0 + 6.0*ct1*x/ct3;
  aa[3] = 1.0 + 6.0*ab/anbn + 3.0*ct1*ct2/ct3;

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
    
    bb[4] = g1*bb[3] + g2*bb[2] + g3*bb[1];
    aa[4] = g1*aa[3] + g2*aa[2] + g3*aa[1];
    
    if(fabs(aa[4]*bb[1]-aa[1]*bb[4]) < EPS*fabs(bb[4]*bb[1])) break;
    
    for(j=1; j<=3; j++) {
      aa[j] = aa[j+1];
      bb[j] = bb[j+1];
    }
  }
  
  *result = aa[4]/bb[4];
  
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
 *
 * If fixed up the window for 1+a-b near zero. [GJ]
 */
int
gsl_sf_hyperg_U_impl(const double a, const double b, const double x, double * result)
{
  const double EPS      = 2.0 * GSL_MACH_EPS;  /* EPS = D1MACH(3) */
  const double SQRT_EPS = M_SQRT2 * GSL_SQRT_MACH_EPS;

  if(x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(locMAX(fabs(a),1.0)*locMAX(fabs(1.0+a-b),1.0) < 0.99 * fabs(x)) {
    double d9;
    int stat_d9 = d9chu(a, b, x, &d9);
    *result = d9 * pow(x, -a);
    return stat_d9;
  }  
  else if(fabs(1.0 + a - b) < SQRT_EPS) {
    /* ALGORITHM IS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X
     */
    /* We can however do the following:
     * U(a,b,x) = U(a,a+1,x) when 1+a-b=0
     * and U(a,a+1,x) = x^(-a).
     */
    double lnx = log(x);
    double lnr = -a * lnx;
    return gsl_sf_exp_impl(lnr, result);
  }
  else {
    double aintb = ( b < 0.0 ? (int)(b-0.5) : (int)(b+0.5) );
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

    double powx;
    double gamr;
    int stat_gamr = gsl_sf_gammainv_impl(1.0+a-b, &gamr);
    int stat_powx = gsl_sf_pow_int_impl(x, istrt);
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
      int stat_pch1i  = gsl_sf_pochrel_impl(xi + 1.0, -beps, &pch1i);
      int stat_poch1bxibeps = gsl_sf_pochrel_impl(b+xi, -beps, &poch1bxibeps);
      double c0 = factor * pochai * gamrni * gamri1
                  * (-poch1bxibeps + pch1ai - pch1i + beps*pch1ai*pch1i);

      /*
       C  XEPS1 = (1.0 - X**(-BEPS))/BEPS = (X**(-BEPS) - 1.0)/(-BEPS)
       */
      double dexprl;
      int stat_dexprl = gsl_sf_expm1_impl(-beps*lnx, &dexprl);
      double xeps1 = lnx * dexprl / (-beps*lnx);

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

