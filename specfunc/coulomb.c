/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_laguerre.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_coulomb.h"


/* normalization for hydrogenic wave functions */
static double R_norm(int n, int l, double Z)
{
  int i;
  double A = 2.*Z/n;
  double term1 = A*A*A /(2.*n);
  double term2 = 1.;
  for(i=n+l; i>=n-l; i--) { term2 *= (double)i; }
  return sqrt(term1/term2);
}

double gsl_sf_hydrogenicR_1(double Z, double r)
{
  double A = 2.*Z;
  double norm = A*sqrt(0.5*A);
  double ea = exp(-Z*r);
  return norm*ea;
}

double gsl_sf_hydrogenicR_2(int l, double Z, double r)
{
  double A = Z;
  double norm = R_norm(2, l, Z);
  double ea = exp(-Z*r);
  if(l == 0) {
    return norm * ea * laguerre_1(2.*l+1, A*r);
  }
  else if(l == 1) {
    double pp = A*r;
    return norm * ea * pp;
  }
  else {
    char buff[100];
    sprintf(buff,"hydrogenicR_2: l= %d", l);
    GSL_ERROR_RETURN(buff, GSL_EDOM, 0.);
  }
}

double gsl_sf_hydrogenicR(int n, int l, double Z, double r)
{
  if(n < 1 || l > n-1 || Z <= 0.) {
    return 0.;
  }
  else {
    double A = 2.*Z/n;
    double norm = R_norm(n, l, Z);
    double rho = A*r;
    double ea = exp(-0.5*rho);
    double pp = pow(rho, l);
    return norm * ea * pp * laguerre_cp(n-l-1, 2.*l+1., rho);
  }
}


/* the L=0 normalization constant */
static double C0sq(double eta)
{
  double twopieta = 2.*M_PI*eta;

  if(fabs(eta) < GSL_MACH_EPS) {
    return 1.;
  }
  else if(twopieta > GSL_LOG_DBL_MAX) {
    return 0.;
  }
  else {
    return twopieta/expm1(twopieta);
  }
}


/* the full definition of C_L(eta) for any valid L and eta
   [Abramowitz and Stegun 14.1.7]
   This depends on the complex gamma function. For large
   arguments the phase of the complex gamma function is not
   very accurately determined. However the modulus is, and that
   is all that we need to calculate C_L.
 */
static double CLeta(double L, double eta)
{
  double ln1; /* log of numerator Gamma function */
  double ln2; /* log of denominator Gamma function */

  if(fabs(eta) < GSL_MACH_EPS) {
    ln1 = gsl_sf_lngamma(L+1.);
  }
  else {
    double p1;  /* phase of numerator Gamma -- not used */
    gsl_sf_complex_lngamma(L+1., eta, &ln1, &p1);
  }
  ln2 = gsl_sf_lngamma(2.*L+2.);
  
  return exp(L*M_LN2 - 0.5*eta*M_PI + ln1 - ln2);
}


double gsl_sf_coulomb_CL(double lam, double eta)
{
  if(lam <= -1.) {
    char buff[100];
    sprintf(buff,"coulomb_CL: lam= %g <= -1", lam);
    GSL_ERROR_RETURN(buff, GSL_EDOM, 0.);
  }
  if(fabs(lam) < 10.*GSL_MACH_EPS) {
    /* saves a calculation of complex_lngamma(),
       otherwise not necessary
     */
    return sqrt(C0sq(eta));
  }
  else {
    return CLeta(lam, eta);
  }
}


int gsl_sf_coulomb_CL_list(double lam_min, int count, double eta, double * cl)
{
  int ell;
  cl[0] = gsl_sf_coulomb_CL(lam_min, eta);

  for(ell=1; ell<count; ell++) {
    double L = lam_min + ell;
    cl[ell] = cl[ell-1] * sqrt(L*L + eta*eta)/(L*(2.*L+1.));
  }
  
  return GSL_SUCCESS;
}



/* junk for coulfg() */
#define CFG_ABORT	 1.e5 /* 2.e4 */
#define CFG_TM30	 1.e-30
#define CFG_RT2DPI	 0.79788456080286535587989211986876373 /* sqrt(2/pi) */

#define CFG_ACCUR	 GSL_MACH_EPS
#define CFG_ACC	         (10.*CFG_ACCUR)
#define CFG_ACC4	 (CFG_ACC*100.*100.)


/* coulfg() mode control */
#define Mode_F   3
#define Mode_FG  2
#define Mode_FGp 1


/* zero argument calculation of coulomb wave functions */
static void coulfg_zero_x(double eta, double xlmin, double xlmax,
			  double *fc, double *gc, double *fcp, double *gcp,
			  int mode)
{
  double delta_lam = xlmax - xlmin + CFG_ACC;
  int  i_delta_lam = (int)delta_lam;

  int i;
  for(i=0; i<=i_delta_lam; i++) { fc[i] = 0.; }
  if(mode==Mode_FGp || mode==Mode_FG) {
    if(fabs(xlmin) > GSL_MACH_EPS || fabs(xlmax) > GSL_MACH_EPS) {
      GSL_ERROR("coulfg_zero_x: x=0.0: G,Gprime undefined for L>0",
	      	GSL_EDOM
	      	);
    }
    else {
      gc[0] = 1./sqrt(C0sq(eta));
    }
  }
  if(mode==Mode_FGp) {
    if(fabs(xlmin) > GSL_MACH_EPS || fabs(xlmax) > GSL_MACH_EPS) {
      fcp[0] = 0.;
    }
    else {
      fcp[0] = sqrt(C0sq(eta));
    }
    gcp[0] = 0.;
  }
}


/* small argument calculation of coulomb wave functions 
   based on expansion in terms of spherical Bessel functions
   [Abramowitz and Stegun 14.4.5]
 */
static void coulfg_small_args(double x, double eta, double xlmin, double xlmax,
			      double *fc, double *gc, double *fcp, double *gcp,
			      int mode)
{
  int i;
  double delta_lam = xlmax - xlmin + CFG_ACC;
  int  i_delta_lam = (int)delta_lam;

  /* everything is small so we use
     F_L ~ C_L(eta) x^(L+1)
     G_L ~ 1/(2L+1) 1/C_L(eta) 1/x^L
     Fp_L ~ (L+1) C_L(eta) x^L = (L+1)/x F_L
     Gp_L ~ -L/(2L+1) 1/C_L(eta) 1/x^(L+1) = -L/x G_L
   */
  double * cl = (double *)malloc((i_delta_lam+1)*sizeof(double));
  if(cl == 0){
    GSL_ERROR("coulfg_small_args: out of memory", GSL_ENOMEM);
    return;
  }
  coulomb_CL_list(xlmin, i_delta_lam+1, eta, cl);
  for(i=0; i<=i_delta_lam; i++) { 
    fc[i] = cl[i] * pow(x,i+1);
    if(mode==Mode_FGp || mode==Mode_FG) {
      gc[i] = 1./(2.*i+1.) /cl[i] /pow(x,i);
    }
    if(mode==Mode_FGp) {
      fcp[i] = fc[i] * (i+1.)/x;
      gcp[i] = gc[i] * (-i)/x;
    }
  }
  free(cl);
}




/* ------------------------------------------------------------ 

 The following hack job is converted from some fortran code.
   There is a small function jwkb() and the main function
   coulfg(). See the comments with coulfg().
   
  ------------------------------------------------------------ */


/* overflow exponent, can be != 0 if scaling required,
   which occurs if the wkb method is invoked and
   it generates large exponents
 */
static int over_exp_ = 0;
int coul_wave_overflow_exp(void) { return over_exp_; }


#define ALOGE  0.4342945
static void jwkb(double x, double eta, double xl,
		 double *fjwkb, double *gjwkb, int *iexp)
{
  /*
    COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS FOR XL.GE. 0
    AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
  */
  double gh2  = x*(2.*eta - x);
  double xll1 = Max(xl*(xl + 1.), 0.);

  if(gh2 + xll1 > 0.0) {   /* if( x < turning point) */
    double hll = xll1 + 6./35.;
    double hl  = sqrt(hll);
    double sl  = eta/hl + hl/x;
    double rl2 = 1.0 + eta*eta/hll;
    double gh  = sqrt(gh2 + hll)/x;
    double phi = x*gh - 0.5*(hl*log((gh+sl)*(gh+sl)/rl2) - log(gh));
    double phi10;

    if(eta != 0.0) phi -= eta*atan2(x*gh,x - eta);
    phi10 = -phi*ALOGE;
    *iexp = (int)phi10;

    /* scale the results if the exponent would be an overflow */
    if(*iexp > 300)  {
      char buff[100];
      sprintf(buff,"coul: overflow in G: scale all F,G by 10^%d", *iexp);
      GSL_ERROR(buff, GSL_EFAILED);
      over_exp_ = *iexp;
      *gjwkb = pow(10., phi10 - *iexp);
    }
    if(*iexp <= 300) {
      *gjwkb = exp(-phi);
      *iexp  = 0;
    }
    *fjwkb = 0.5/(gh * *gjwkb);
  }
  else {
    GSL_ERROR("jwkb: called above turning point: INTERNAL ERROR",
	      GSL_EFAILED
	      );
  }
}
#undef ALOGE


/*
 SUBROUTINE COULFG(X,ETA,XLMIN,XLMAX, FC,GC,FCP,GCP, MODE1,KFN,IFAIL)

C  REVISED IJT WITH L-T ALGORITHMN FOR CONTINUED FRACTIONS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								       C
C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD	       C
C								       C
C  A. R. BARNETT	   MANCHESTER  MARCH   1981		       C
C								       C
C  ORIGINAL PROGRAM 'RCWFN'	 IN    CPC  8 (1974) 377-395	       C
C		  + 'RCWFF'	 IN    CPC 11 (1976) 141-142	       C
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314	       C
C  THIS VERSION WRITTEN UP	 IN    CPC 27 (1982) 147-166	       C
C								       C
C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
C   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    C
C   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   C
C   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     C
C   THE DIRAC EQUATION ,ALSO SPHERICAL & CYLINDRICAL BESSEL EQUATIONS  C
C								       C
C  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    C
C  STARTING ARRAY ELEMENT IS M1 = MAX0(IDINT(XLMIN+ACCUR),0) + 1       C
C      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES		       C
C								       C
C  IF 'MODE' = 1  GET F,G,F',G'	  FOR INTEGER-SPACED LAMBDA VALUES     C
C	     = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    C
C	     = 3      F		      CALL TO AT LEAST LENGTH (1)      C
C  IF 'KFN'  = 0 REAL	     COULOMB FUNCTIONS ARE RETURNED	       C
C	     = 1 SPHERICAL   BESSEL	 "	"     "		       C
C	     = 2 CYLINDRICAL BESSEL	 "	"     "		       C
C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT			       C
C								       C
C  PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    C
C   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     C
C   COULFG IS CODED FOR REAL*8 ON IBM OR EQUIVALENT  ACCUR = 10**-16   C
C   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   C
C   FOR MANTISSAS OF 56 & 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) C
C   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*/

/* I turned off the KFN input, since this will only be used
   for calculating Coulomb wave functions explicitly. Of course,
   there is no reason this cannot be called with eta=0, in which
   case Bessel functions is what you get...
 */

/* This troglodyte code makes some sense if you read the paper:
   Barnett, Comp. Phys. Comm., 21, 297 (1981).
   But not much.
 */

/* Behavioural change: the original fortran version would load
   the computed values into the arrays starting at an index given
   by the lowest requested l value (integer part). This is dumb.
   I changed it so that it always loads them starting at index 0.
   You can easily get the other effect by passing an offset
   pointer if you like.
 */

static void coulfg(double x, double eta, double xlmin, double xlmax,
		   double *fc, double *gc, double *fcp, double *gcp,
		   int mode)
{
  double acch  = sqrt(CFG_ACC);

  double delta_lam = xlmax - xlmin + CFG_ACC;
  int  i_delta_lam = (int)delta_lam;
  double   upper_l = xlmin + i_delta_lam;

  int iexp = 1;
  double fjwkb = 0.0;
  double gjwkb = 0.0;

  if(xlmin <= -1. || xlmax < xlmin) {
    char buff[100];
    sprintf(buff, "coulfg: problem with lambda inputs: %g %g", xlmax, xlmin);
    GSL_ERROR(buff, GSL_EDOM);
    return;
  }

  if(x < 0.) {
    int i;
    char buff[100];
    sprintf(buff,"coulfg: x= %g   < 0", x);
    GSL_ERROR(buff, GSL_EDOM);
    for(i=0; i<=i_delta_lam; i++) {
      fc[i] = 0.;
      if(mode==Mode_FGp || mode==Mode_FG) {
	gc[i] = 0.;
      }
      if(mode==Mode_FGp) {
	fcp[i] = 0.;
	gcp[i] = 0.;
      }
    }
  }
  else if(x == 0.) {
    coulfg_zero_x(eta, xlmin, xlmax, fc, gc, fcp, gcp, mode);
  }
  else if(x < 0.001 && fabs(eta*x) < 0.01) {
    coulfg_small_args(x, eta, xlmin, xlmax, fc, gc, fcp, gcp, mode);
  }
  else if(x < acch) {
    char buff[100];
    sprintf(buff,"coulfg: x= %g  eta= %g  x.eta= %g  not yet implemented",
	    x, eta, x*eta);
    GSL_ERROR(buff, GSL_EDOM);
  }
  else {
    /* check if x is below the turning point */
    int xlturn = (x*(x - 2.*eta) < xlmin*(xlmin + 1.) ? 1 : 0);

    double e2mm1 = eta*eta + xlmin*(xlmin + 1.);

    /* starting index */
    /* int M1  = Max((int)(xlmin + CFG_ACC),0) */ /* + 1 */
    int M1 = 0;

    double x_inv = 1./x;
    double fcl = 1.0;
    double gcl;
    double pk  = upper_l + 1.0;
    double px  = pk  + CFG_ABORT;
    double df;

    int L;
    double xl;
    double P;
    double Q;

    double F = eta/pk + pk*x_inv;
    double D;
    double C;

    double alpha;
    double beta;
    double gam;
    double fcm;
    double W;
    double gpl;

    if(fabs(F) < CFG_TM30) F = CFG_TM30;
    D = 0.0;
    C = F;

    if(fabs(delta_lam - floor(delta_lam+0.5)) > 100.*CFG_ACC) {
      char buff[100];
      sprintf(buff, "coulfg: xlmax-xlmin= %25.18g  not an integer",
	      delta_lam);
      GSL_ERROR_MESSAGE(buff, GSL_EDOM);
    }

    /* Compute first continued fraction, evaluating F_prime()/F()
       at the upper lambda value.
     */
    do {
      double pk1 = pk + 1.0;
      double ek  = eta / pk;
      double rk2 = 1. + ek*ek;
      double tk  = (pk + pk1)*(x_inv + ek/pk1);
      D   =  tk - rk2 * D;
      C   =  tk - rk2 / C;
      if(fabs(C) < CFG_TM30) C = CFG_TM30;
      if(fabs(D) < CFG_TM30) D = CFG_TM30;
      D = 1.0/D;
      df = D * C;
      F  = F * df; /* the result ???...  I guess so */
      if(D < 0.) fcl = -fcl;
      pk = pk1;
      if( pk > px ) {
	char buff[100];
	sprintf(buff, "coulfg: first continued fraction not converging");
	GSL_ERROR(buff, GSL_EFAILED);
	return;
      }
    }
    while(fabs(df-1.) > CFG_ACC);

    if(i_delta_lam != 0) {
      /* downward recurrence from the upper l value to the minimum
	 requested l point
	 array gc[] (if present) stores rl
       */
      int lp;
      double fpl;
      fcl *= CFG_TM30;
      fpl = fcl*F;
      if(mode == Mode_FGp) fcp[M1 + i_delta_lam] = fpl;
      fc[M1 + i_delta_lam] = fcl;
      xl  = upper_l;
      for(lp=1; lp<=i_delta_lam; lp++) {
	double el = eta/xl;
	double rl = sqrt(1. + el*el);
	double sl =  el  + xl*x_inv;
	double fc_lm1 = (fcl *sl + fpl)/rl;
	fpl   =  fc_lm1*sl - fcl *rl;
	fcl   =  fc_lm1;
	L     =  M1 + i_delta_lam  - lp;
	fc[L] =  fcl;
	if(mode == Mode_FGp) fcp[L]  = fpl;
	if(mode != Mode_F  /* && eta_ne_zero */) gc[L+1] = rl;
	--xl;
      }
      
      if(fcl == 0.) fcl = CFG_ACC;
      F	= fpl/fcl;
    }

    /*
      NOW WE HAVE REACHED LAMBDA = xlmin = XLM
      EVALUATE CF2 = P + I.Q	 AGAIN USING STEED'S ALGORITHM
     */
    
    if(xlturn) jwkb(x, eta, Max(xlmin,0.), &fjwkb, &gjwkb, &iexp);
    
    /* if(iexp != 1) fprintf(stderr,"iexp= %d\n", iexp); */

    if(iexp > 1 || gjwkb > 1.0/(acch*100.)) {
      /* ARRIVE HERE IF G(xlmin) > 10**6 OR iexp > 70 and xlturn = true */
      W	 = fjwkb;
      gam = gjwkb*W;
      P	  = F;
      Q	  = 1.;
    }
    else {
      double ta =  2.0*CFG_ABORT;
      double wi =  2.*eta;
      double ar = -e2mm1;
      double ai =  eta;
      double br =  2.*(x - eta);
      double bi =  2.;
      double dr =  br/(br*br + bi*bi);
      double di = -bi/(br*br + bi*bi);

      double dp = -x_inv*(ar*di + ai*dr);
      double dq =  x_inv*(ar*dr - ai*di);
      
      double A;
      double B;

      xlturn = 0;
      pk =  0.0;
      P	 =  0.0;
      Q	 =  1. - eta*x_inv;
      
      do {
	P += dp;
	Q += dq;
	pk += 2.;
	ar += pk;
	ai += wi;
	bi += 2.;
	D  = ar*dr - ai*di + br;
	di = ai*dr + ar*di + bi;
	C  = 1./(D*D + di*di);
	dr =  C*D;
	di = -C*di;
	A  = br*dr - bi*di - 1.;
	B  = bi*dr + br*di;
	C  = dp*A  - dq*B;
	dq = dp*B  + dq*A;
	dp = C;
	if(pk > ta) {
	  char buff[100];
	  sprintf(buff,"coulfg: second continued fraction not converging");
	  GSL_ERROR(buff, GSL_EFAILED);
	  return;
	}
      }
      while(fabs(dp)+fabs(dq) >= (fabs(P)+fabs(Q))*CFG_ACC);
      
      /*
	SOLVE FOR fcm = F AT LAMBDA = xlmin,THEN FIND NORM FACTOR W=W/fcm
       */
      gam = (F - P)/Q;
      if(Q <= CFG_ACC4*fabs(P)) {
	GSL_ERROR("coulfg: final Q < abs(P)*CFG_ACC*1e4", GSL_EFAILED);
      }
      W	= 1.0/sqrt((F - P)*gam + Q);
    }
    
    /*
      NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
     */
    alpha = 0.0;
    beta  = 1.0;
    fcm   = ( fcl > 0. ? W : -W) * beta; /* fcm   = DSIGN(W,fcl)*beta; */
    fc[M1] = fcm;

    if(mode == Mode_FGp || mode == Mode_FG) {
      if(! xlturn)   gcl =  fcm*gam;
      if(  xlturn)   gcl =  gjwkb*beta;
      gc[M1]  = gcl;
      gpl =  gcl*(P - Q/gam) - alpha*gcl;
    }
    if(mode == Mode_FGp) {
      gcp[M1] = gpl;
      fcp[M1] = fcm*(F - alpha);
    }
    
    if(i_delta_lam == 0 ) return;
    
    /*
      UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS rl
      RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
      XL   = xlmin HERE  AND rl = 1. , el = 0. FOR BESSELS
     */
    W = beta*W/fabs(fcl);

    for(L=M1; L<=M1+i_delta_lam-1; L++) {

      double gcl1;
      
      xl += 1.0;

      if(mode == Mode_FGp || mode == Mode_FG) {
	/* if(eta_ne_zero) */ double el = eta/xl;
	/* if(eta_ne_zero) */ double rl = gc[L+1];
	double sl = el + xl*x_inv;
	gcl1  = ((sl - alpha)*gcl - gpl)/rl;
	gpl   = rl*gcl - (sl + alpha)*gcl1;
	gcl   = gcl1;
	gc[L+1] = gcl1;
      }
      if(mode == Mode_FGp) {
	gcp[L+1] = gpl;
	fcp[L+1] = W*(fcp[L+1] - alpha*fc[L+1]);
      }
      fc[L+1] = W * fc[L+1];
    }
  }
}



void gsl_sf_coulomb_wave_F(double x, double eta,
		     double lam_min, double lam_max,
		     double * fc)
{
  coulfg(x, eta, lam_min, lam_max,
	 fc, (double *) 0, (double *) 0, (double *) 0,
	 Mode_F
	 );
}

void gsl_sf_coulomb_wave_FG(double x, double eta,
		     double lam_min, double lam_max,
		     double * fc, double * gc)
{
  coulfg(x, eta, lam_min, lam_max,
	 fc, gc, (double *) 0, (double *) 0,
	 Mode_FG
	 );
}

void gsl_sf_coulomb_wave_FGp(double x, double eta,
		      double lam_min, double lam_max,
		      double * fc, double * gc,
		      double * fcp, double * gcp)

{
  coulfg(x, eta, lam_min, lam_max,
	 fc, gc, fcp, gcp,
	 Mode_FGp
	 );
}


int gsl_sf_coulomb_wave_sphF(double x, double eta,
		      double lam_min, double lam_max,
		      double * fc)
{
  int mode;
  double delta_lam = lam_max - lam_min + CFG_ACC;
  int i_delta_lam  = (int)delta_lam;
  int i;

  if(fabs(x) <= 100.*GSL_MACH_EPS) {
    /* if x==0 then it is nonzero and finite only for l=0 */
    double ell = lam_min;
    mode = 0;
    for(i=0; i<=i_delta_lam; i++) {
      if(fabs(ell) < 1000.*GSL_MACH_EPS) {
	fc[i] = sqrt(C0sq(eta));
      }
      else {
	fc[i] = 0.;
      }
      ell += 1.;
    }
  }
  else {
    mode = 1;
    coulfg(x, eta, lam_min, lam_max,
	   fc, (double *) 0, (double *) 0, (double *) 0,
	   Mode_F
	   );
    for(i=0; i<=i_delta_lam; i++) { fc[i] = fc[i] / x; }
  }
  return mode;
}
