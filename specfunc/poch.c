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


/* ((a)_x - 1)/x in the "small x" region where
 * cancellation must be controlled.
 *
 * Based on SLATEC DPOCH1().
 */
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
static
int
pochrel_smallx(const double a, const double x, double * result)
{
  /*
   SQTBIG = 1.0D0/SQRT(24.0D0*D1MACH(1))
   ALNEPS = LOG(D1MACH(3))
   */
  const double SQTBIG = 1.0/(2.0*M_SQRT2*M_SQRT3*GSL_SQRT_DBL_MIN);
  const double ALNEPS = GSL_LOG_DBL_EPSILON - M_LN2;

  if(x == 0.0) {
    return gsl_sf_psi_impl(a, result);
  }
  else {
    const double bp   = (  (a < -0.5) ? 1.0-a-x : a );
    const int	 incr = ( (bp < 10.0) ? 11.0-bp : 0 );
    const double b    = bp + incr;
    double dpoch1;
    double dexprl;
    int stat_dexprl;
    int i;

    double var    = b + 0.5*(x-1.0);
    double alnvar = log(var);
    double q = x*alnvar;

    double poly1 = 0.0;

    if(var < SQTBIG) {
      const int nterms = (int)(-0.5*ALNEPS/alnvar + 1.0);
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
        /* nterms = 20; */
	*result = 0.0;
	return GSL_ESANITY;
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

    for(i=incr-1; i >= 0; i--) {
      /*
       C WE HAVE DPOCH1(B,X), BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
       C TO OBTAIN DPOCH1(BP,X).
       */
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


/* Handle the case of a>0 and a+x>0.
 * Assumes neither is equal to zero.
 */
static
int
lnpoch_pos(const double a, const double x, double * result)
{
  double absx = fabs(x);

  if(absx > 0.1*a || absx*log(GSL_MAX(a,2.0)) > 0.1) {
    double g1, g2;
    gsl_sf_gammainv_impl(a,   &g1);
    gsl_sf_gammainv_impl(a+x, &g2);
    *result = -log(g2/g1);
    return GSL_SUCCESS;
  }
#if 0
  else if(fabs(x) < 0.1*a) {
    /* Be careful about the implied subtraction.
     * Note that both a+x and and a must both be
     * large here since they are not both small
     * (the above test failed) and x is not
     * relatively large. So we calculate
     * using Stirling for Log[Gamma(z)],
     *
     *   Log[Gamma(a+x)/Gamma(a)] = x(Log[a]-1) + (x+a-1/2)Log[1+x/a]
     *                              + (1/(1+eps)   - 1) / (12 a)
     *                              - (1/(1+eps)^3 - 1) / (360 a^3)
     *                              + (1/(1+eps)^5 - 1) / (1260 a^5)
     *                              - (1/(1+eps)^7 - 1) / (1680 a^7)
     *                              + ...
     */
    const double eps = x/a;
    const double den = 1.0 + eps;
    const double d3 = den*den*den;
    const double d5 = d3*den*den;
    const double d7 = d5*den*den;
    const double c1 = -eps/den;
    const double c3 = -eps*(3.0+eps*(3.0+eps))/d3;
    const double c5 = -eps*(5.0+eps*(10.0+eps*(10.0+eps*(5.0+eps))))/d5;
    const double c7 = -eps*(7.0+eps*(21.0+eps*(35.0+eps*(35.0+eps*(21.0+eps*(7.0+eps))))))/d7;
    const double a4 = a*a*a*a;
    const double a6 = a4*a*a;
    const double ser = (c1 + c3/(30.0*a*a) + c5/(105.0*a4) +c7/(140.0*a6))/(12.0*a);

    double term1 = x * (log(a) - 1.0);
    double term2;
    double ln_1peps;
    gsl_sf_log_1plusx_impl(eps, &ln_1peps);  /* log(1 + x/a) */
    term2 = (x + a - 0.5) * ln_1peps;

    *result = term1 + term2 + ser;
    return GSL_SUCCESS;
  }
#endif
  else {
    double poch_rel;
    int stat_p = pochrel_smallx(a, x, &poch_rel);
    double eps = x*poch_rel;
    int stat_e = gsl_sf_log_1plusx_impl(eps, result);
    return GSL_ERROR_SELECT_2(stat_e, stat_p);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

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
  else {
    return lnpoch_pos(a, x, result);
  }
}


int
gsl_sf_lnpoch_sgn_impl(const double a, const double x,
                       double * result, double * sgn)
{
  if(a == 0.0 || a+x == 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(a > 0.0 && a+x > 0.0) {
    return lnpoch_pos(a, a+x, result);
  }
  else if(a < 0.0 && a+x < 0.0) {
    /* Reduce to positive case using reflection.
     */
    double sin_1 = sin(M_PI * (1.0 - a));
    double sin_2 = sin(M_PI * (1.0 - a - x));
    if(sin_1 == 0.0 || sin_2 == 0.0) {
      *result = 0.0;
      *sgn = 0.0;
      return GSL_EDOM;
    }
    else {
      double lnp_pos;
      int stat_pp = lnpoch_pos(1.0-a, -x, &lnp_pos);
      *result = log(fabs(sin_1/sin_2)) - lnp_pos;
      *sgn = GSL_SIGN(sin_1*sin_2);
      return stat_pp;
    }
  }
  else {
    /* Evaluate gamma ratio directly.
     */
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
}


int
gsl_sf_poch_impl(const double a, const double x, double * result)
{
  if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else {
    double lnpoch;
    double sgn;
    int stat_lnpoch = gsl_sf_lnpoch_sgn_impl(a, x, &lnpoch, &sgn);
    int stat_exp    = gsl_sf_exp_impl(lnpoch, result);
    *result *= sgn;
    return GSL_ERROR_SELECT_2(stat_exp, stat_lnpoch);
  }
}


int
gsl_sf_pochrel_impl(const double a, const double x, double * result)
{
  const double absx = fabs(x);
  const double absa = fabs(a);

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
    return pochrel_smallx(a, x, result);
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
