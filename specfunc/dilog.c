/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_integration.h>
#include "gsl_sf_dilog.h"


/* data for summation loop for dilog()
 * CC(N)=(N(N+1)(N+2))**2  [with index shift for conversion to C]
 */
static double CC[30] = {
  36.,         576.,       36.e+2,      144.e+2,     441.e+2,
  112896.,     254016.,    5184.e+2,    9801.e+2,    17424.e+2,
  2944656.,    4769856.,   74529.e+2,   112896.e+2,  166464.e+2,
  23970816.,   33802596.,  467856.e+2,  636804.e+2,  853776.e+2,
  112911876.,  147476736., 19044.e+4,   24336.e+4,   3080025.e+2,
  386358336.,  480661776., 5934096.e+2, 7273809.e+2, 8856576.e+2
};

/* summation loop used by most cases in dilog() */
static double do_sum(double y, double by, double dl)
{
  int i;
  double a, test;
  double b = 4.*y*y/by;
  for(i=0; i<30; i++) {
    b = b*y;
    a = b/CC[i];
    test = dl;
    dl += a;
    if(dl == test) return dl;
  }
  return dl;
}


/* based on:
   ALGORITHM 490 COLLECTED ALGORITHMS FROM ACM.
   ALGORITHM APPEARED IN COMM. ACM, VOL. 18, NO. 4N, P. 200.
   DOUBLE PRECISION FUNCTION DILOG(X)
   used: ZERO OF DILOG ON THE POSITIVE REAL AXIS, X0=12.59517...

   I checked this with Mathematica Re[PolyLog[2,x]]; it works. [GJ]
*/
double gsl_sf_dilog(double x)
{
  if(x > 12.6) {
    double log_x = log(x);
    double Y  = 1./x;
    double BY = -1. - Y*(4.+Y);
    double dl = 3.28986813369645287 -
      0.5*log_x*log_x + (Y*(4.+5.75*Y)+3.*(1.+Y)*(1.-Y)*log(1.-Y))/BY;
    if(dl + 4.*Y == dl)
      return dl;
    else
      return do_sum(Y, BY, dl);
  }
  else if(x >= 12.59) {
    /* DILOG COMPUTED FROM TAYLOR SERIES ABOUT ZERO OF DILOG, X0. */
    double X0 = 12.5951703698450184;
    double Y  = x/X0 - 1.;
    double Z  = 1./11.5951703698450184;
    double W  = Y*Z;
    double C1 = (3.*X0-2.)/6.;
    double C2 = ((11.*X0-15.)*X0+6.)/24.;
    double C3 = (((50.*X0-104.)*X0+84.)*X0-24.)/120.;
    double C4 = ((((274.*X0-770.)*X0+940.)*X0-540.)*X0+120.)/720.;
    return Y*(1.-Y*(0.5-Y*(1./3.-Y*(0.25-Y*(.2-Y/6.)))))*log(Z)
      - W*X0*Y*(0.5-W*(C1-W*(C2-W*(C3-W*C4))));
  }
  else if(x >= 2.) {
    /* same as first case... */
    double log_x = log(x);
    double Y  = 1./x;
    double BY = -1. - Y*(4.+Y);
    double dl = 3.28986813369645287 -
      0.5*log_x*log_x + (Y*(4.+5.75*Y)+3.*(1.+Y)*(1.-Y)*log(1.-Y))/BY;
    if(dl + 4.*Y == dl)
      return dl;
    else
      return do_sum(Y, BY, dl);
  }
  else if(x > 1.) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(7) WITH
     * X=1/X + EQ(6), AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
     */
    double Y = 1. - 1./x;
    double DX = log(x);
    double BY = 1. + Y*(4.+Y);
    double dl = 1.64493406684822643 + DX*(0.5*DX-log(x-1.)) 
      + (Y*(4.+5.75*Y)-3.*(1.+Y)*DX/x)/BY;
    return do_sum(Y, BY, dl);
  }
  else if(x == 1.) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(2). */
    return 1.64493406684822643;
  }
  else if(x > 0.5) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(7),
     * AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
     */
    double Y = 1. - x;
    double DX = log(x);
    double BY = -1. - Y*(4.+Y);
    double dl = 1.64493406684822643 
      - DX*log(Y) + (Y*(4.+5.75*Y)+3.*(1.+Y)*DX*x)/BY;
    return do_sum(Y, BY, dl);
  }
  else if(x > 0.01) {
    /* DILOG COMPUTED FROM DESCRIPTION OF THIS ALGORITHM, EQ(4) */
    double Y = x;
    double BY = 1. + Y*(4.+Y);
    double dl = (Y*(4.+5.75*Y)+3.*(1.+Y)*(1.-Y)*log(1.-Y))/BY;
    return do_sum(Y, BY, dl);
  }
  else if(x < -1.) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.245, EQ(12) WITH
     * X=-X, AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
     */
    double Y = 1./(1.-x);
    double DX = log(-x);
    double DY = log(Y);
    double BY = 1. + Y*(4.+Y);
    double dl = -1.64493406684822643 +
      0.5*DY*(DY+2.*DX) + (Y*(4.+5.75*Y) + 3.*(1.+Y)*(1.-Y)*(DX+DY))/BY;
    if(dl + 4.*Y == dl)
      return dl;
    else
      return do_sum(Y, BY, dl);
  }
  else if(x < -0.01) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(8),
     * AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
     */
    double Y = x/(x-1.);
    double DX = log(1.-x);
    double BY = -1. - Y*(4.+Y);
    double dl = (Y*(4.+5.75*Y)-3.*(1.+Y)*(1.-Y)*DX)/BY - 0.5*DX*DX;
    return do_sum(Y, BY, dl);
  }
  else {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(1). */
    return x*(1.+
	      x*(0.25+x*(1./9.+
			 x*(0.0625+x*(4.e-2+
				      x*(1./36.+x*(1./49.+x/64.)))))));
  }
}


/* dilog_c() for small modulus */
static inline double dilog_c_smallr(double r, double cos_theta)
{
  return cos_theta * r + 0.25 * (2.*cos_theta*cos_theta - 1.) * r*r;
}


/* integrand in definition of dilog_c() */
static double dummy_cos_theta;
static double dilog_c_integrand(double s)
{
  return -0.5 * log(1. - 2.*dummy_cos_theta*s + s*s) / s;
}


/* calculate the defining integral for the (incomplete) dilog_c() */
static double dilog_c_integral(double a, double b, double cos_theta, double eps)
{
  set_integ_info("called from dilog_c_integral()");
  dummy_cos_theta = cos_theta;
  return simpson(dilog_c_integrand, a, b, eps);
}


double gsl_sf_dilog_c(double r, double cos_theta)
{
  if(r < 0.) {
    GSL_MESSAGE("dilog_c: r < 0");
    return 0.;
  }
  else if(fabs(cos_theta) > 1.) {
    GSL_MESSAGE("dilog_c: fabs(cos_theta) > 1");
    return 0.;
  }
  else if(r == 0.) {
    return 0.;
  }
  else if(fabs(r-1.) < 1.e-10) {
    double arg = M_PI - acos(cos_theta);
    return 0.25 * arg*arg - M_PI*M_PI/12.;
  }
  else {
    /* restrict to (0,1) using the standard functional relation */
    double t = ( r > 1. ? 1./r : r);
    double tmp_result;
    
    /* combination of series expansion and integral (if necessary);
     * the series is used to handle very small t, which could upset
     * the integrator 
     */
    double lo_cut = 1.e-6;
    if(t <= lo_cut) {
      tmp_result  = dilog_c_smallr(t, cos_theta);
    }
    else {
      tmp_result  = dilog_c_smallr(lo_cut, cos_theta);
      tmp_result += dilog_c_integral(lo_cut, t, cos_theta, lo_cut);
    }
    
    if(r > 1.) {
      double lgr = log(r);
      double arg = M_PI - acos(cos_theta);
      return -tmp_result - 0.5*lgr*lgr + 0.5*arg*arg - M_PI*M_PI/6.;
    }
    else {
      return tmp_result;
    }
  }
}
