/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_log.h"
#include "gsl_sf_psi.h"
#include "gsl_sf_gamma.h"


int
gsl_sf_lnpoch_impl(const double a, const double x, double * result)
{
  if(a <= 0.0 || a+x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(fabs(x) < 0.01*a) {
    double gs_ax, gs_a;
    double ln_1pxa, ln_1pxa_mxa;
    gsl_sf_gammastar_impl(a+x, &gs_ax);
    gsl_sf_gammastar_impl(a,   &gs_a);
    gsl_sf_log_1plusx_impl(x/a, &ln_1pxa);         /* log(1 + x/a)       */
    gsl_sf_log_1plusx_mx_impl(x/a, &ln_1pxa_mxa);  /* log(1 + x/a) - x/a */
    *result = x*log(a)+ (x-0.5)*ln_1pxa + a*ln_1pxa_mxa + log(gs_ax/gs_a);
    return GSL_SUCCESS;
  }
  else {
    double lg_apn, lg_a;
    int stat_apn = gsl_sf_lngamma_impl(a+x, &lg_apn);
    int stat_a   = gsl_sf_lngamma_impl(a,   &lg_a);
    if(stat_apn == GSL_SUCCESS && stat_a == GSL_SUCCESS) {
      *result = lg_apn - lg_a;
      return GSL_SUCCESS;
    }
    else if(stat_apn == GSL_EDOM || stat_a == GSL_EDOM){
      *result = 0.0;
      return GSL_EDOM;
    }
    else {
      *result = 0.0;
      return GSL_FAILURE;
    }
  }
}


int
gsl_sf_lnpoch_sgn_impl(const double a, const double x,
                       double * result, double * sgn)
{
  double lg_apn, lg_a;
  double s_apn, s_a;
  int stat_apn = gsl_sf_lngamma_sgn_impl(a+x, &lg_apn, &s_apn);
  int stat_a   = gsl_sf_lngamma_sgn_impl(a,   &lg_a,   &s_a);
  if(stat_apn == GSL_SUCCESS && stat_a == GSL_SUCCESS) {
    *result = lg_apn - lg_a;
    *sgn = s_a * s_apn;
    return GSL_SUCCESS;
  }
  else if(stat_apn == GSL_EDOM || stat_a == GSL_EDOM){
    *result = 0.0;
    *sgn = 0.0;
    return GSL_EDOM;
  }
  else {
    *result = 0.0;
    *sgn = 0.0;
    return GSL_FAILURE;
  }
}


int
gsl_sf_poch_impl(const double a, const double x, double * result)
{
  double lnpoch;
  double sgn;
  int stat_lnpoch = gsl_sf_lnpoch_sgn_impl(a, x, &lnpoch, &sgn);
  if(stat_lnpoch != GSL_SUCCESS) {
    *result = 0.0;
    return stat_lnpoch;
  }
  else {
    double abspoch;
    int stat_exp = gsl_sf_exp_impl(lnpoch, &abspoch);
    *result = sgn * abspoch;
    return stat_exp;
  }
}



/* Based on SLATEC DPOCH1(). */

const static double bern[21] = {
   0.0   /* no element 0 */,  
  +0.833333333333333333333333333333333e-01,
  -0.138888888888888888888888888888888e-02,
  +0.330687830687830687830687830687830e-04,
  -0.826719576719576719576719576719576e-06,
  +0.208767569878680989792100903212014e-07,
  -0.528419013868749318484768220217955e-09,
  +0.133825365306846788328269809751291e-10,
  -0.338968029632258286683019539124944e-12,
  +0.858606205627784456413590545042562e-14,
  -0.217486869855806187304151642386591e-15,
  +0.550900282836022951520265260890225e-17,
  -0.139544646858125233407076862640635e-18,
  +0.353470703962946747169322997780379e-20,
  -0.895351742703754685040261131811274e-22,
  +0.226795245233768306031095073886816e-23,
  -0.574472439520264523834847971943400e-24,
  +0.145517247561486490186626486727132e-26,
  -0.368599494066531017818178247990866e-28,
  +0.933673425709504467203255515278562e-30,
  -0.236502241570062993455963519636983e-31
};


/*
C When ABS(X) is so small that substantial cancellation will occur if
C the straightforward formula is used, we use an expansion due
C to Fields and discussed by Y. L. Luke, The Special Functions and Their
C Approximations, Vol. 1, Academic Press, 1969, page 34.
C
C The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
C        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) .
C In order to maintain significance in POCH1, we write for positive a
C        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q)
C                       = 1.0 + Q*EXPREL(Q) .
C Likewise the polynomial is written
C        POLY = 1.0 + X*POLY1(A,X) .
C Thus,
C        POCH1(A,X) = (POCH(A,X) - 1) / X
C                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X)
C
*/
int
gsl_sf_pochrel_impl(const double a, const double x, double * result)
{
  /*
   SQTBIG = 1.0D0/SQRT(24.0D0*D1MACH(1))
   ALNEPS = LOG(D1MACH(3))
   */
  const double SQTBIG = 1.0/(2.0*M_SQRT2*M_SQRT3*GSL_SQRT_DBL_MIN);
  const double ALNEPS = GSL_LOG_MACH_EPS - M_LN2;

  if(x == 0.0) {
    return gsl_sf_psi_impl(a, result);
  }
  else {
    double absx = fabs(x);
    double absa = fabs(a);
    
    if(absx > 0.1*absa || absx*log(GSL_MAX(absa,2.0)) > 0.1) {
      double lnpoch, sgn;
      int stat_poch = gsl_sf_lnpoch_sgn_impl(a, x, &lnpoch, &sgn);
      if(lnpoch > GSL_LOG_DBL_MAX) {
        *result = 0.0;
	return GSL_EOVRFLW;
      }
      else {
        *result = (sgn*exp(lnpoch) - 1.0)/x;
	return stat_poch;
      }
    }
    else {
      double q;
      double var, alnvar;
      double poly1;
      double dpoch1;
      const double bp   = (  (a < -0.5) ? 1.0-a-x : a );
      const int    incr = ( (bp < 10.0) ? 11.0-bp : 0 );
      const double b    = bp + incr;
      double dexprl;
      int stat_dexprl;
      
      int i;

      var    = b + 0.5*(x-1.0);
      alnvar = log(var);
      q = x*alnvar;

      poly1 = 0.0;
      
      if(var < SQTBIG) {
	int nterms = (int)(-0.5*ALNEPS/alnvar + 1.0);
        const double var2 = (1.0/var)/var;
	const double rho  = 0.5 * (x + 1.0);
	double term = var2;
        double gbern[24];
	int k, j;

	gbern[1] = 1.0;
	gbern[2] = -rho/12.0;
	poly1 = gbern[2] * term;

	if(nterms > 20) {
	  /* NTERMS IS TOO BIG, MAYBE D1MACH(3) IS BAD */
	  nterms = 20;
	  /* FIXME: something to do here...? */
	}

	for(k=2; k<=nterms; k++) {
	  double gbk = 0.0;
	  for(j=1; j<=k; j++) {
	    gbk += bern[k-j+1]*gbern[j];
	  }
	  gbern[k+1] = -rho*gbk/k;
	  
	  term  *= (2*k-2-x)*(2*k-1-x)*var2;
	  poly1 += gbern[k+1]*term;
	}
      }

      stat_dexprl = gsl_sf_expm1_impl(q, &dexprl);
      if(stat_dexprl != GSL_SUCCESS) {
        *result = 0.0;
	return stat_dexprl;
      }
      dexprl = dexprl/q;
      poly1 *= (x - 1.0);
      dpoch1 = dexprl * (alnvar + q * poly1) + poly1;

      /* for(ii=1; ii<=incr; ii++) { */
      for(i=incr-1; i >= 0; i--) {
        /*
         C WE HAVE DPOCH1(B,X), BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
         C TO OBTAIN DPOCH1(BP,X).
         */
	/* i = incr - ii; */
	double binv = 1.0/(bp+i);
	dpoch1 = (dpoch1 - binv) / (1.0 + x*binv);
      }

      if(bp == a) {
        *result = dpoch1;
	return GSL_SUCCESS;
      }
      else {
        /*
         C WE HAVE DPOCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A
         C REFLECTION FORMULA TO OBTAIN DPOCH1(A,X).
         */
        double sinpxx = sin(M_PI*x)/x;
        double sinpx2 = sin(0.5*M_PI*x);
        double trig   = sinpxx/tan(M_PI*b) - 2.0*sinpx2*(sinpx2/x);
        *result = dpoch1 * (1.0 + x*trig) + trig;
	return GSL_SUCCESS;
      }

    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_lnpoch_e(const double a, const double x, double * result)
{
  int status = gsl_sf_lnpoch_impl(a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_lnpoch_e", status);
  }
  return status;
}


int
gsl_sf_lnpoch_sgn_e(const double a, const double x, double * result, double * sgn)
{
  int status = gsl_sf_lnpoch_sgn_impl(a, x, result, sgn);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_lnpoch_sgn_e", status);
  }
  return status;
}


int
gsl_sf_poch_e(const double a, const double x, double * result)
{
  int status = gsl_sf_poch_impl(a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_poch", status);
  }
  return status;
}


int gsl_sf_pochrel_e(const double a, const double x, double * result)
{
  int status = gsl_sf_pochrel_impl(a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_pochrel_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*/

double gsl_sf_lnpoch(const double a, const double x)
{
  double y;
  int status = gsl_sf_lnpoch_impl(a, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_lnpoch", status);
  }
  return y;
}


double gsl_sf_poch(const double a, const double x)
{
  double y;
  int status = gsl_sf_poch_impl(a, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_poch", status);
  }
  return y;
}

double gsl_sf_pochrel(const double a, const double x)
{
  double y;
  int status = gsl_sf_pochrel_impl(a, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_pochrel", status);
  }
  return y;
}
