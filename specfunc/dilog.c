/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_errno.h>
#include <gsl_math.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_clausen.h"
#include "gsl_sf_expint.h"
#include "gsl_sf_trig.h"
#include "gsl_sf_log.h"
#include "gsl_sf_dilog.h"

#define locMin(a,b)  ((a) < (b) ? (a) : (b))
#define locMax(a,b)  ((a) > (b) ? (a) : (b))


/* Chebyshev fit for f(t), -1 < t < 1
 * f(t) := g(y(t))
 * g(y) := Li2(e^y) - y^2 / 2
 * ln(100) < y < ln(1000)
 *     100 < x < 1000
 */
static double li2_100_data[16] = {
 -61.028079324811232797820044037,
 -13.250456556257821054795103038,
 -0.6639108587048147840321307640,
  0.0002200540327603079836383037,
 -0.0000313868258373686202867481,
  3.6260369837091164247821887091e-06,
 -3.5496362345786979782455944857e-07,
  3.0620871708238888177065740622e-08,
 -2.4246508584548732667041886525e-09,
  1.8455504629240591893145007108e-10,
 -1.4152531131794409561489951395e-11,
  1.1297627664081362099337605665e-12,
 -9.4449717256599517834598188335e-14,
  8.1522854869580797300554951414e-15,
 -7.1420495403769695224676811122e-16,
  6.2912377518369338425394488383e-17
};
static struct gsl_sf_cheb_series li2_100_cs = {
  li2_100_data,
  15,
  -1, 1,
  (double *)0,
  (double *)0
};


/* Chebyshev fit for f(t), -1 < t < 1
 * f(t) := g(y(t))
 * g(y) := Li2(e^y) - y^2 / 2
 * ln(12.8) < y < ln(100)
 *     12.8 < x < 100
 */
static double li2_128_data[20] = {
 -20.143648665645142102724844394,
 -7.3205904240532221814737357600,
 -0.5366068386310610280100336649,
  0.0014474301595135692246987985,
 -0.0001955700828883060185739608,
  0.0000224631238841696967842219,
 -2.3705284953339821946312727933e-06,
  2.4678700660662671110541595250e-07,
 -2.6602697657502708483759115435e-08,
  3.0086702014507668782188160992e-09,
 -3.5329565849464834813281331964e-10,
  4.2506826109785483748269984540e-11,
 -5.2021743313231527081025406483e-12,
  6.4580207435845740872418744638e-13,
 -8.1198790560514584576622308967e-14,
  1.0324722357928494025287106650e-14,
 -1.3256357964592402022034221203e-15,
  1.7163997468657032520536099126e-16,
 -2.2381415252860542554097578635e-17,
  2.8879257066187275357654284735e-18
};
static struct gsl_sf_cheb_series li2_128_cs = {
  li2_128_data,
  19,
  -1, 1,
  (double *)0,
  (double *)0
};


/* data for summation loop for dilog()
 * (n(n+1)(n+2))^2
 */
const static double n_np1_np2_sq[30] = {
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
  double a;
  double b = 4.0*y*y/by;
  double previous_dl;
  int i;
  for(i=0; i<30; i++) {
    b *= y;
    a  = b/n_np1_np2_sq[i];
    previous_dl = dl;
    dl += a;
    if(fabs(a) < GSL_MACH_EPS) break;
  }
  return dl;
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



/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_dilog_impl(const double x, double * result)
{
  if(x > 1000.0) {
    /* Li2(x) = Li2(A) + (ln(A)^2 - ln(x)^2)/2 + S(A) - S(x)
     *
     * S(s) = 1/s + 1/2^2 1/s^2 + 1/3^2 1/s^3 + ...
     *
     * A = 10.0
     */
    const double ix = 1.0/x;
    const double lx = log(x);
    const double ln_A  = 2.3025850929940456840;  /* log(10) */
    const double li2_A = 0.5363012873578627366;  /* Li2(10) */
    const double S_A   = 0.1026177910993911311;
    const double S_x   = ix * (1.0 + ix*(0.25 + ix*(1.0/9.0 + ix*(1.0/16.0 + ix/25.0))));
    *result = li2_A + 0.5*ln_A*ln_A + S_A - S_x - 0.5*lx*lx;
    return GSL_SUCCESS;
  }
  else if(x > 100.0) {
    const double a = 4.60517018598809136803598;  /* log(100)  */
    const double b = 6.90775527898213705205397;  /* log(1000) */
    const double y = log(x);
    const double t = 2.0*(y-a)/(b-a) - 1.0;
    const double c = gsl_sf_cheb_eval(&li2_100_cs, t);
    *result = c + 0.5*y*y;
    return GSL_SUCCESS;
  }
  else if(x > 12.8) {
    const double a = 2.54944517092557148190263; /* log(12.8) */
    const double b = 4.60517018598809136803598; /* log(100)  */
    const double y = log(x);
    const double t = 2.0*(y-a)/(b-a) - 1.0;
    const double c = gsl_sf_cheb_eval(&li2_128_cs, t);
    *result = c + 0.5*y*y;
    return GSL_SUCCESS;
  }
  else if(x > 12.4) {
    const double x0 = 12.5951703698450184;   /* positive real zero of dilog(x) */
    const double dx = x-x0;
    const double c1 = -0.19456574163185890;
    const double c2 =  0.004300177565288112;
    const double c3 = -0.0001291882631106331;
    const double c4 =  3.44864872694839e-06;
    const double c5 =  5.66899694544e-10;
    const double c6 = -1.26641834906114e-08;
    const double c7 =  1.63966793864394e-09;
    const double c8 = -1.64221074630073e-10;
    const double c9 =  1.49644905020987e-11;
    *result = dx*(c1 + dx*(c2 + dx*(c3 + dx*(c4 + dx*(c5 + dx*(c6 + dx*(c7 + dx*(c8 + dx*c9))))))));
    return GSL_SUCCESS;
  }
  else if(x >= 2.0) {
    const double log_x = log(x);
    const double inv_x = 1.0/x;
    const double den   = -(1.0 + 4.0*inv_x + inv_x*inv_x);
    const double term1 = 0.25*inv_x * (16.0 + 23.0*inv_x);
    const double term2 = (1.0+inv_x)*(1.0-inv_x) * log(1.0-inv_x);
    const double c0 = 3.28986813369645287;
    const double dl = c0 - 0.5*log_x*log_x + (term1+3.0*term2)/den;
    *result = do_sum(inv_x, den, dl);
    return GSL_SUCCESS;
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
  else if(x > -0.01) {
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
  else if(x > -1.0) {
    /* DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(1). */
    *result = x*(1.0+
	         x*(0.25+x*(1.0/9.0+
			    x*(0.0625+x*(4.0e-2+
				         x*(1.0/36.0+x*(1.0/49.0+x/64.0)))))));
    return GSL_SUCCESS;
  }
  else {
    /* x < -1.0 */
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


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_dilog_e(const double x, double * result)
{
  int status = gsl_sf_dilog_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_dilog_e", status);
  }
  return status;
}


int
gsl_sf_complex_dilog_e(const double x, const double y, double * result_re, double * result_im)
{
  int status = gsl_sf_complex_dilog_impl(x, y, result_re, result_im);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_dilog_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*/

double
gsl_sf_dilog(const double x)
{
  double y;
  int status = gsl_sf_dilog_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_dilog", status);
  }
  return y;
}
