/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_errno.h>
#include <gsl_math.h>
#include "gsl_sf_clausen.h"
#include "gsl_sf_expint.h"
#include "gsl_sf_trig.h"
#include "gsl_sf_log.h"
#include "gsl_sf_dilog.h"

#define locMin(a,b)  ((a) < (b) ? (a) : (b))
#define locMax(a,b)  ((a) > (b) ? (a) : (b))


/* data for summation loop for dilog()
 * CC(N)=(N(N+1)(N+2))**2  [with index shift for conversion to C]
 */
static double CC[30] = {
  36.0,         576.0,       36.0e+2,      144.0e+2,     441.0e+2,
  112896.0,     254016.0,    5184.0e+2,    9801.0e+2,    17424.0e+2,
  2944656.0,    4769856.0,   74529.0e+2,   112896.0e+2,  166464.0e+2,
  23970816.0,   33802596.0,  467856.0e+2,  636804.0e+2,  853776.0e+2,
  112911876.0,  147476736.0, 19044.0e+4,   24336.0e+4,   3080025.0e+2,
  386358336.0,  480661776.0, 5934096.0e+2, 7273809.0e+2, 8856576.0e+2
};

/* summation loop used by most cases in dilog() */
static double do_sum(const double y, const double by, double dl)
{
  int i;
  double a;
  double b = 4.0*y*y/by;
  double pre_dl;
  for(i=0; i<30; i++) {
    b = b*y;
    a = b/CC[i];
    pre_dl = dl;
    dl += a;
    if(dl == pre_dl) break;
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
int
gsl_sf_dilog_impl(const double x, double * result)
{
  if(x > 12.6) {
    double log_x = log(x);
    double Y  = 1.0/x;
    double BY = -1.0 - Y*(4.+Y);
    double dl = 3.28986813369645287 -
      0.5*log_x*log_x + (Y*(4.+5.75*Y)+3.*(1.+Y)*(1.-Y)*log(1.-Y))/BY;
    if(dl + 4.0*Y == dl) {
      *result = dl;
      return GSL_SUCCESS;
    }
    else {
      *result = do_sum(Y, BY, dl);
      return GSL_SUCCESS;
    }
  }
  else if(x >= 12.59) {
    /* DILOG COMPUTED FROM TAYLOR SERIES ABOUT ZERO OF DILOG, X0. */
    double X0 = 12.5951703698450184;
    double Y  = x/X0 - 1.;
    double Z  = 1.0/11.5951703698450184;
    double W  = Y*Z;
    double C1 = (3.*X0-2.)/6.;
    double C2 = ((11.*X0-15.)*X0+6.)/24.;
    double C3 = (((50.*X0-104.)*X0+84.)*X0-24.)/120.;
    double C4 = ((((274.*X0-770.)*X0+940.)*X0-540.)*X0+120.)/720.;
    *result = Y*(1.-Y*(0.5-Y*(1./3.-Y*(0.25-Y*(.2-Y/6.)))))*log(Z)
      - W*X0*Y*(0.5-W*(C1-W*(C2-W*(C3-W*C4))));
    return GSL_SUCCESS;
  }
  else if(x >= 2.0) {
    /* same as first case... */
    double log_x = log(x);
    double Y  = 1.0/x;
    double BY = -1.0 - Y*(4.0+Y);
    double dl = 3.28986813369645287 -
      0.5*log_x*log_x + (Y*(4.+5.75*Y)+3.*(1.+Y)*(1.-Y)*log(1.-Y))/BY;
    if(dl + 4.0*Y == dl) {
      *result = dl;
      return GSL_SUCCESS;
    }
    else {
      *result = do_sum(Y, BY, dl);
      return GSL_SUCCESS;
    }
  }
  else if(x > 1.0) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(7) WITH
     * X=1/X + EQ(6), AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
     */
    double Y = 1.0 - 1.0/x;
    double DX = log(x);
    double BY = 1.0 + Y*(4.+Y);
    double dl = 1.64493406684822643 + DX*(0.5*DX-log(x-1.)) 
      + (Y*(4.0+5.75*Y)-3.*(1.+Y)*DX/x)/BY;
    *result = do_sum(Y, BY, dl);
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(2). */
    *result = 1.64493406684822643;
    return GSL_SUCCESS;
  }
  else if(x > 0.5) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(7),
     * AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
     */
    double Y = 1.0 - x;
    double DX = log(x);
    double BY = -1.0 - Y*(4.0+Y);
    double dl = 1.64493406684822643 
      - DX*log(Y) + (Y*(4.0+5.75*Y)+3.0*(1.0+Y)*DX*x)/BY;
    *result = do_sum(Y, BY, dl);
    return GSL_SUCCESS;
  }
  else if(x > 0.01) {
    /* DILOG COMPUTED FROM DESCRIPTION OF THIS ALGORITHM, EQ(4) */
    double Y = x;
    double BY = 1.0 + Y*(4.0+Y);
    double dl = (Y*(4.0+5.75*Y)+3.0*(1.0+Y)*(1.0-Y)*log(1.0-Y))/BY;
    *result = do_sum(Y, BY, dl);
    return GSL_SUCCESS;
  }
  else if(x < -1.0) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.245, EQ(12) WITH
     * X=-X, AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
     */
    double Y = 1.0/(1.0-x);
    double DX = log(-x);
    double DY = log(Y);
    double BY = 1.0 + Y*(4.0+Y);
    double dl = -1.64493406684822643 +
      0.5*DY*(DY+2.0*DX) + (Y*(4.0+5.75*Y) + 3.0*(1.0+Y)*(1.0-Y)*(DX+DY))/BY;
    if(dl + 4.0*Y == dl) {
      *result = dl;
      return GSL_SUCCESS;
    }
    else {
      *result = do_sum(Y, BY, dl);
      return GSL_SUCCESS;
    }
  }
  else if(x < -0.01) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(8),
     * AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
     */
    double Y = x/(x-1.0);
    double DX = log(1.0-x);
    double BY = -1.0 - Y*(4.0+Y);
    double dl = (Y*(4.0+5.75*Y)-3.0*(1.0+Y)*(1.0-Y)*DX)/BY - 0.5*DX*DX;
    *result = do_sum(Y, BY, dl);
    return GSL_SUCCESS;
  }
  else {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(1). */
    *result = x*(1.0+
	         x*(0.25+x*(1.0/9.0+
			    x*(0.0625+x*(4.0e-2+
				         x*(1.0/36.0+x*(1.0/49.0+x/64.0)))))));
    return GSL_SUCCESS; 
  }
}


/* Evaluate the series representation for Li2(z):
 *
 *   Li2(z) = Sum[ |z|^k / k^2 Exp[i k arg(z)], {k,1,Infinity}]
 *   |z|    = r
 *   arg(z) = theta
 *   
 * Assumes 0 < r < 1. 
 */
static
int
dilogc_series_1(double r, double cos_theta, double sin_theta,
                double * real_result, double * imag_result)
{
  double alpha = 1.0 - cos_theta;
  double beta  = sin_theta;
  double ck = cos_theta;
  double sk = sin_theta;
  double rk = r;
  double real_sum = r*ck;
  double imag_sum = r*sk;
  int kmax = 50 + (int)(20.0/(-log(r)));
  int k;
  for(k=2; k<kmax; k++) {
    double ck_tmp = ck;
    ck = ck - (alpha*ck + beta*sk);
    sk = sk - (alpha*sk - beta*ck_tmp);
    rk *= r;
    real_sum += rk/((double)k*k) * ck;
    imag_sum += rk/((double)k*k) * sk;
  }
  
  *real_result = real_sum;
  *imag_result = imag_sum;
  return GSL_SUCCESS;
}


/* Evaluate a series for Li_2(z) when |z| is near 1.
 * This is uniformly good away from z=1.
 *
 *   Li_2(z) = Sum[ a^n/n! H_n(theta), {n, 0, Infinity}]
 *
 * where
 *   H_n(theta) = Li_2(Exp(I theta))
 *   a = ln(r)
 *
 *  H_0(t) = Gl_2(t) + i Cl_2(t)
 *  H_1(t) = 1/2 ln(2(1-c)) + I atan2(-s, 1-c)
 *  H_2(t) = -1/2 + I/2 s/(1-c)
 *  H_3(t) = -1/2 /(1-c)
 *  H_4(t) = -I/2 s/(1-c)^2
 *  H_5(t) = 1/2 (2 + c)/(1-c)^2
 *  H_6(t) = I/2 s/(1-c)^5 (8(1-c) - s^2 (3 + c))
 *
 *  assumes: 0 <= theta <= 2Pi
 */
static
int
dilogc_series_2(double r, double theta, double cos_theta, double sin_theta,
                double * real_result, double * imag_result)
{
  double a = log(r);
  double omc = 1.0 - cos_theta;
  double H_re[7];
  double H_im[7];
  double an, nfact;
  double sum_re, sum_im;
  int n;

  H_re[0] = M_PI*M_PI/6.0 + 0.25*(theta*theta - 2.0*M_PI*fabs(theta));
  gsl_sf_clausen_impl(theta, &(H_im[0]));
  
  H_re[1] = 0.5*log(2.0*omc);
  H_im[1] = atan2(-sin_theta, omc);
  
  H_re[2] = -0.5;
  H_im[2] = 0.5 * sin_theta/omc;
  
  H_re[3] = -0.5/omc;
  H_im[3] = 0.0;
  
  H_re[4] = 0.0;
  H_im[4] = -0.5*sin_theta/(omc*omc);
  
  H_re[5] = 0.5 * (2.0 + cos_theta)/(omc*omc);
  H_im[5] = 0.0;
  
  H_re[6] = 0.0;
  H_im[6] = 0.5 * sin_theta/(omc*omc*omc*omc*omc)
            * (8*omc - sin_theta*sin_theta*(3 + cos_theta));
 
  sum_re = H_re[0];
  sum_im = H_im[0];
  an = 1.0;
  nfact = 1.0;
  for(n=1; n<=6; n++) {
    double t;
    an *= -a;
    nfact *= n;
    t = an/nfact;
    sum_re += t * H_re[n];
    sum_im += t * H_im[n];
  }
  
  *real_result = sum_re;
  *imag_result = sum_im;
  return GSL_SUCCESS;
}


/* complex dilogarithm in the unit disk
 * assumes:  r < 1  and  0 <= theta <= 2Pi
 */
static
int
complex_dilog_unitdisk(double r, double theta,
                       double * real_dl, double * imag_dl)
{
  const double zeta2 = M_PI*M_PI/6.0;
  int stat_dilog;
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  double x = r * cos_theta;
  double y = r * sin_theta;
  double x_tmp, y_tmp, r_tmp;
  double result_re_tmp, result_im_tmp;

  /* Reflect away from z = 1 if
   * we are too close.
   */
  if(x > 0.5) {
    x_tmp = 1.0 - x;
    y_tmp = -y;
    r_tmp = sqrt(x_tmp*x_tmp + y_tmp*y_tmp);
  }
  else {
    x_tmp = x;
    y_tmp = y;
    r_tmp = r;
  }

  /* Calculate dilog of the transformed variable.
   */
  if(r_tmp < 0.98) {
    double cos_theta_tmp = x_tmp / r_tmp;
    double sin_theta_tmp = y_tmp / r_tmp;
    stat_dilog = dilogc_series_1(r_tmp, cos_theta_tmp, sin_theta_tmp,
				 &result_re_tmp, &result_im_tmp
				 );
  }
  else {
    double cos_theta_tmp = x_tmp / r_tmp;
    double sin_theta_tmp = y_tmp / r_tmp;
    double theta_tmp = atan2(y_tmp, x_tmp);
    stat_dilog = dilogc_series_2(r_tmp, theta_tmp, cos_theta_tmp, sin_theta_tmp,
                                 &result_re_tmp, &result_im_tmp
				 );
  }

  /* Unwind reflection if necessary.
   *
   * Li2(z) = -Li2(1-z) + zeta(2) - ln(z) ln(1-z)
   */
  if(x > 0.5) {
    double lnz    =  log(r);                 /*  log(|z|)   */
    double lnomz  =  log(r_tmp);             /*  log(|1-z|) */
    double argz   =  theta;                  /*  arg(z)     */
    double argomz =  atan2(y_tmp, x_tmp);    /*  arg(1-z)   */
    *real_dl = -result_re_tmp + zeta2 - lnz*lnomz + argz*argomz;
    *imag_dl = -result_im_tmp - argz*lnomz - argomz*lnz;
  }
  else {
    *real_dl = result_re_tmp;
    *imag_dl = result_im_tmp;
  }

  return stat_dilog;
}


int
gsl_sf_complex_dilog_impl(const double r, double theta,
                          double * real_dl, double * imag_dl)
{
  if(theta < 0.0 || theta > 2.0*M_PI) {
    gsl_sf_angle_restrict_pos_impl(&theta);
  }

  if(r == 0.0) {
    *real_dl = 0.0;
    *imag_dl = 0.0;
    return GSL_SUCCESS;
  }

  /* Trap cases of real-valued argument.
   */
  if(theta == 0.0) {
    *imag_dl = ( r > 1.0 ? -M_PI * log(r) : 0.0 );
    return gsl_sf_dilog_impl(r, real_dl);
  }
  if(theta == M_PI) {
    *imag_dl = 0.0;
    return gsl_sf_dilog_impl(-r, real_dl);
  }

  /* Trap unit circle case.
   */
  if(r == 1.0) {
    *real_dl = M_PI*M_PI/6.0 + 0.25*(theta*theta - 2.0*M_PI*fabs(theta));
    return gsl_sf_clausen_impl(theta, imag_dl);
  }

  /* Generic case.
   */
  {
    int stat_dilog;
    double r_tmp, theta_tmp;
    double result_re_tmp, result_im_tmp;

    /* Reduce argument to unit disk.
     */
    if(r > 1.0) {
      r_tmp     = 1.0 / r;
      theta_tmp = 2.0*M_PI - theta;
    }
    else {
      r_tmp     = r;
      theta_tmp = theta;
    }

    /* Calculate in the unit disk.
     */
    stat_dilog = complex_dilog_unitdisk(r_tmp, theta_tmp,
                                        &result_re_tmp, &result_im_tmp
					);

    /* Unwind the inversion if necessary. We calculate
     * the imaginary part explicitly if using the inversion
     * because there is no simple relationship between
     * arg(1-z) and arg(1 - 1/z), which is the "omega"
     * term in [Lewin A.2.5 (1)].
     */
    if(r > 1.0) {
      const double zeta2 = M_PI*M_PI/6.0;
      double x = r * cos(theta);
      double y = r * sin(theta);
      double omega = atan2(y, 1.0-x);
      double lnr = log(r);
      double Cl_a, Cl_b, Cl_c;
      double pmt = M_PI - theta;
      gsl_sf_clausen_impl(2.0*omega, &Cl_a);
      gsl_sf_clausen_impl(2.0*theta, &Cl_b);
      gsl_sf_clausen_impl(2.0*(omega+theta), &Cl_c);
      *real_dl = -result_re_tmp - 0.5*lnr*lnr + 0.5*pmt*pmt - zeta2;
      *imag_dl = omega*log(r) + 0.5*(Cl_a + Cl_b - Cl_c);
    }
    else {
      *real_dl = result_re_tmp;
      *imag_dl = result_im_tmp;
    }
    
    return stat_dilog;
  }
}


int gsl_sf_dilog_e(const double x, double * result)
{
  int status = gsl_sf_dilog_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_dilog_e", status);
  }
  return status;
}

int gsl_sf_complex_dilog_e(const double x, const double y, double * result_re, double * result_im)
{
  int status = gsl_sf_complex_dilog_impl(x, y, result_re, result_im);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_dilog_e", status);
  }
  return status;
}


double gsl_sf_dilog(const double x)
{
  double y;
  int status = gsl_sf_dilog_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_dilog", status);
  }
  return y;
}
