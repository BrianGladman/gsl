/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_hyperg.h"



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


#if 0
/* Based on SLATEC DCHU() [W. Fullerton]
 * -  C This routine is not valid when 1+A-B is close to zero if X is small.
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
    /* ALGORITHM IS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X */
    *result = 0.0;
    return GSL_EFAILED;
  }
  else {
    double aintb = ( b < 0.0 ? (int)(b-0.5) : (int)(b+0.5) );
    double beps  = b - aintb;
    int N = aintb;
    
    double lnx  = log(x);
    double xeps = exp(-beps*lnx);
    
    double sum;
   
    /* Evaluate finite sum.
     */
    if(N <= 0) {
      double t = 1.0;
      int i;
      sum = 1.0;
      for(i=1; i<= -N; i++) {
	double xi1 = i - 1;
	t   *= (a+xi1)*x/((b+xi1)*(xi1+1.0));
	sum += t;
      }
      sum = sum * DPOCH(1.0+a-b, -a);
    }
    else {
      int M = N - 2;
      if(M < 0) {
        sum = 0.0;
      }
      else {
        double t = 1.0;
	sum = 1.0;
	for(i=1; i<=M; i++) {
	  t   *= (a-b+i)*x/((1.0-b+i)*i);
	  sum += t;
	}
	sum = sum * DGAMMA(b-1.0) * DGAMR(a) * pow(x,1.0-N) * xeps;
      }
    }
      

    /* Evaluate infinite sum.
     */
    istrt = 0;
    if(N < 1) istrt = 1 - N;
    xi = istrt;
    
    factor = (-1.0)**N * DGAMR(1.0D0+A-B) * X**ISTRT;
    if(beps != 0.0) factor *= beps*M_PI/sin(beps*M_PI);
    
    POCHAI = DPOCH (A, XI);
    GAMRI1 = DGAMR (XI+1.0D0);
    GAMRNI = DGAMR (AINTB+XI);
    B0 = FACTOR * DPOCH(A,XI-BEPS) * GAMRNI * DGAMR(XI+1.0-BEPS);


    if(fabs(xeps-1.0) < 0.5) {
/*    
    C
C X**(-BEPS) IS CLOSE TO 1.0D0, SO WE MUST BE CAREFUL IN EVALUATING THE
C DIFFERENCES.
C
*/
      PCH1AI = DPOCH1 (A+XI, -BEPS);
      PCH1I = DPOCH1 (XI+1.0D0-BEPS, BEPS);
      C0 = FACTOR * POCHAI * GAMRNI * GAMRI1 * (
       -DPOCH1(B+XI,-BEPS) + PCH1AI - PCH1I + BEPS*PCH1AI*PCH1I);

/*C
C XEPS1 = (1.0 - X**(-BEPS))/BEPS = (X**(-BEPS) - 1.0)/(-BEPS)
*/
      XEPS1 = ALNX*DEXPRL(-BEPS*ALNX);

      DCHU = SUM + C0 + XEPS1*B0;
      XN = N;
      
      for(i=1; i<1000; i++) {
     
        XI = ISTRT + I;
        XI1 = ISTRT + I - 1;
        B0 = (A+XI1-BEPS)*B0*X/((XN+XI1)*(XI-BEPS));
        C0 = (A+XI1)*C0*X/((B+XI1)*XI)
     	 - ((A-1.0D0)*(XN+2.D0*XI-1.0D0) + XI*(XI-BEPS)) * B0
     	 / (XI*(B+XI1)*(A+XI1-BEPS));
        T = C0 + XEPS1*B0;
        DCHU = DCHU + T;
        IF (ABS(T).LT.EPS*ABS(DCHU)) GO TO 130
     }
     /*
      CALL XERMSG ('SLATEC', 'DCHU',
     +   'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES', 3, 2)
     */
    }
    else {
    /*
    C
C X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE STRAIGHTFORWARD
C FORMULATION IS STABLE.
C
*/
      A0 = FACTOR * POCHAI * DGAMR(B+XI) * GAMRI1 / BEPS;
      B0 = XTOEPS * B0 / BEPS;

      DCHU = SUM + A0 - B0;
      
      for(i=1; i<1000; i++) {
      
        XI = ISTRT + I;
        XI1 = ISTRT + I - 1;
        A0 = (A+XI1)*A0*X/((B+XI1)*XI);
        B0 = (A+XI1-BEPS)*B0*X/((AINTB+XI1)*(XI-BEPS));
        T = A0 - B0;
        DCHU = DCHU + T;
        IF (ABS(T).LT.EPS*ABS(DCHU)) GO TO 130
	
      }
      /*
      CALL XERMSG ('SLATEC', 'DCHU',
     +   'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES', 3, 2)  
      */
    }

    /* 130  RETURN */
  }
}
#endif
