/* specfunc/dilog.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_clausen.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_dilog.h>


/* Evaluate series for real dilog(x)
 * Sum[ x^k / k^2, {k,1,Infinity}]
 *
 * Converges rapidly for |x| < 1/2.
 */
static
int
dilog_series(const double x, gsl_sf_result * result)
{
  const int kmax = 1000;
  double sum  = x;
  double term = x;
  int k;
  for(k=2; k<kmax; k++) {
    double rk = (k-1.0)/k;
    term *= x;
    term *= rk*rk;
    sum += term;
    if(fabs(term/sum) < GSL_DBL_EPSILON) break;
  }

  result->val  = sum;
  result->err  = 2.0 * fabs(term);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  if(k == kmax)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


/* Assumes x >= 0.0
 */
static
int
dilog_xge0(const double x, gsl_sf_result * result)
{
  if(x > 2.0) {
    const double log_x = log(x);
    gsl_sf_result ser;
    int stat_ser = dilog_series(1.0/x, &ser);
    double t1 = M_PI*M_PI/3.0;
    double t2 = ser.val;
    double t3 = 0.5*log_x*log_x;
    result->val  = t1 - t2 - t3;
    result->err  = GSL_DBL_EPSILON * fabs(log_x) + ser.err;
    result->err += GSL_DBL_EPSILON * (fabs(t1) + fabs(t2) + fabs(t3));
    result->val += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_ser;
  }
  else if(x > 1.01) {
    const double log_x    = log(x);
    const double log_term = log_x * (log(1.0-1.0/x) + 0.5*log_x);
    gsl_sf_result ser;
    int stat_ser = dilog_series(1.0 - 1.0/x, &ser);
    double t1 = M_PI*M_PI/6.0;
    double t2 = ser.val;
    double t3 = log_term;
    result->val  = t1 + t2 - t3;
    result->err  = GSL_DBL_EPSILON * fabs(log_x) + ser.err;
    result->err += GSL_DBL_EPSILON * (fabs(t1) + fabs(t2) + fabs(t3));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_ser;
  }
  else if(x > 1.0) {
    /* series around x = 1.0 */
    const double eps = x - 1.0;
    const double lne = log(eps);
    const double c0 = M_PI*M_PI/6.0;
    const double c1 =   1.0 - lne;
    const double c2 = -(1.0 - 2.0*lne)/4.0;
    const double c3 =  (1.0 - 3.0*lne)/9.0;
    const double c4 = -(1.0 - 4.0*lne)/16.0;
    const double c5 =  (1.0 - 5.0*lne)/25.0;
    const double c6 = -(1.0 - 6.0*lne)/36.0;
    const double c7 =  (1.0 - 7.0*lne)/49.0;
    const double c8 = -(1.0 - 8.0*lne)/64.0;
    result->val = c0+eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*(c7+eps*c8)))))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    result->val = M_PI*M_PI/6.0;
    result->err = 2.0 * GSL_DBL_EPSILON * M_PI*M_PI/6.0;
    return GSL_SUCCESS;
  }
  else if(x > 0.5) {
    const double log_x = log(x);
    gsl_sf_result ser;
    int stat_ser = dilog_series(1.0-x, &ser);
    double t1 = M_PI*M_PI/6.0;
    double t2 = ser.val;
    double t3 = log_x*log(1.0-x);
    result->val  = t1 - t2 - t3;
    result->err  = GSL_DBL_EPSILON * fabs(log_x) + ser.err;
    result->err += GSL_DBL_EPSILON * (fabs(t1) + fabs(t2) + fabs(t3));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_ser;
  }
  else if(x > 0.0) {
    return dilog_series(x, result);
  }
  else {
    /* x == 0.0 */
    result->val = 0.0;
    result->err = 0.0;
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
                gsl_sf_result * real_result, gsl_sf_result * imag_result)
{
  double alpha = 1.0 - cos_theta;
  double beta  = sin_theta;
  double ck = cos_theta;
  double sk = sin_theta;
  double rk = r;
  double real_sum = r*ck;
  double imag_sum = r*sk;
  int kmax = 50 + (int)(22.0/(-log(r))); /* tuned for double-precision */
  int k;
  for(k=2; k<kmax; k++) {
    double ck_tmp = ck;
    ck = ck - (alpha*ck + beta*sk);
    sk = sk - (alpha*sk - beta*ck_tmp);
    rk *= r;
    real_sum += rk/((double)k*k) * ck;
    imag_sum += rk/((double)k*k) * sk;
  }

  real_result->val = real_sum;
  real_result->err = 2.0 * kmax * GSL_DBL_EPSILON * fabs(real_sum);
  imag_result->val = imag_sum;
  imag_result->err = 2.0 * kmax * GSL_DBL_EPSILON * fabs(imag_sum);

  return GSL_SUCCESS;
}


/* Evaluate a series for Li_2(z) when |z| is near 1.
 * This is uniformly good away from z=1.
 *
 *   Li_2(z) = Sum[ a^n/n! H_n(theta), {n, 0, Infinity}]
 *
 * where
 *   H_n(theta) = Sum[ e^(i m theta) m^n / m^2, {m, 1, Infinity}]
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
                gsl_sf_result * real_result, gsl_sf_result * imag_result)
{
  double a = log(r);
  double omc = 1.0 - cos_theta;
  double H_re[7];
  double H_im[7];
  double an, nfact;
  double sum_re, sum_im;
  gsl_sf_result Him0;
  int n;

  H_re[0] = M_PI*M_PI/6.0 + 0.25*(theta*theta - 2.0*M_PI*fabs(theta));
  gsl_sf_clausen_e(theta, &Him0);
  H_im[0] = Him0.val;

  H_re[1] = -0.5*log(2.0*omc);
  H_im[1] = -atan2(-sin_theta, omc);

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
    an *= a;
    nfact *= n;
    t = an/nfact;
    sum_re += t * H_re[n];
    sum_im += t * H_im[n];
  }

  real_result->val = sum_re;
  real_result->err = 2.0 * 6.0 * GSL_DBL_EPSILON * fabs(sum_re) + fabs(an/nfact);
  imag_result->val = sum_im;
  imag_result->err = 2.0 * 6.0 * GSL_DBL_EPSILON * fabs(sum_im) + Him0.err + fabs(an/nfact);

  return GSL_SUCCESS;
}


/* complex dilogarithm in the unit disk
 * assumes:  r < 1  and  0 <= theta <= 2Pi
 */
static
int
dilogc_unitdisk(double r, double theta, gsl_sf_result * real_dl, gsl_sf_result * imag_dl)
{
  const double zeta2 = M_PI*M_PI/6.0;
  int stat_dilog;
  gsl_sf_result cos_theta;
  gsl_sf_result sin_theta;
  int stat_cos = gsl_sf_cos_e(theta, &cos_theta);
  int stat_sin = gsl_sf_sin_e(theta, &sin_theta);
  gsl_sf_result x;
  gsl_sf_result y;
  gsl_sf_result x_tmp, y_tmp, r_tmp;
  gsl_sf_result result_re_tmp, result_im_tmp;
  double cos_theta_tmp;
  double sin_theta_tmp;
  x.val = r * cos_theta.val;
  x.err = r * cos_theta.err;
  y.val = r * sin_theta.val;
  y.err = r * sin_theta.err;

  /* Reflect away from z = 1 if
   * we are too close.
   */
  if(x.val > 0.5) {
    x_tmp.val = 1.0 - x.val;
    x_tmp.err = GSL_DBL_EPSILON * (1.0 + fabs(x.val)) + x.err;
    y_tmp.val = -y.val;
    y_tmp.err = y.err;
    r_tmp.val = sqrt(x_tmp.val*x_tmp.val + y_tmp.val*y_tmp.val);
    r_tmp.err = (x_tmp.err*fabs(x_tmp.val) + y_tmp.err*fabs(y_tmp.val))/fabs(r_tmp.val);
  }
  else {
    x_tmp.val = x.val;
    x_tmp.err = x.err;
    y_tmp.val = y.val;
    y_tmp.err = y.err;
    r_tmp.val = r;
    r_tmp.err = r * GSL_DBL_EPSILON;
  }

  cos_theta_tmp = x_tmp.val / r_tmp.val;
  sin_theta_tmp = y_tmp.val / r_tmp.val;

  /* Calculate dilog of the transformed variable.
   */
  if(r_tmp.val < 0.98) {
    stat_dilog = dilogc_series_1(r_tmp.val, cos_theta_tmp, sin_theta_tmp,
                                 &result_re_tmp, &result_im_tmp
                                 );
  }
  else {
    double theta_tmp = atan2(y_tmp.val, x_tmp.val);
    stat_dilog = dilogc_series_2(r_tmp.val, theta_tmp, cos_theta_tmp, sin_theta_tmp,
                                 &result_re_tmp, &result_im_tmp
                                 );
  }

  /* Unwind reflection if necessary.
   *
   * Li2(z) = -Li2(1-z) + zeta(2) - ln(z) ln(1-z)
   */
  if(x.val > 0.5) {
    double lnz    =  log(r);                         /*  log(|z|)   */
    double lnomz  =  log(r_tmp.val);                 /*  log(|1-z|) */
    double argz   =  theta;                          /*  arg(z)     */
    double argomz =  atan2(y_tmp.val, x_tmp.val);    /*  arg(1-z)   */
    real_dl->val  = -result_re_tmp.val + zeta2 - lnz*lnomz + argz*argomz;
    real_dl->err  =  result_re_tmp.err;
    real_dl->err +=  GSL_DBL_EPSILON * (zeta2 + fabs(lnz*lnomz) + fabs(argz*argomz));
    real_dl->err +=  2.0 * GSL_DBL_EPSILON * fabs(real_dl->val);
    imag_dl->val  = -result_im_tmp.val - argz*lnomz - argomz*lnz;
    imag_dl->err  =  result_im_tmp.err;
    imag_dl->err +=  GSL_DBL_EPSILON * (fabs(argz*lnomz) + fabs(argomz*lnz));
    imag_dl->err +=  2.0 * GSL_DBL_EPSILON * fabs(imag_dl->val);
  }
  else {
    real_dl->val = result_re_tmp.val;
    real_dl->err = result_re_tmp.err;
    imag_dl->val = result_im_tmp.val;
    imag_dl->err = result_im_tmp.err;
  }

  return GSL_ERROR_SELECT_3(stat_dilog, stat_sin, stat_cos);
}



/* This is the real part of Li_2(r e^{i theta}).
 */
static
int
dilog_li2_real(double r, double theta, gsl_sf_result * result)
{
  /* Reduce argument to within the unit disk. */
  const double r_tmp = ( r > 1.0 ? 1.0/r : r );
  double result_tmp;

  if(r_tmp < GSL_SQRT_DBL_EPSILON)
  {
    const double c = cos(theta);
    result_tmp = r_tmp * (c + r_tmp * (2.0*c*c - 1.0));
  }
  else if(r_tmp < 0.98)
  {
  }
  else
  {
  }

  return 0;
}



/* This is the function
 *
 *   Cl_2(2 omega) + Cl_2(2 theta) - Cl_2(2 omega + 2 theta)
 *
 * where
 *
 *   theta := atan2(-A sin(omega), 1 - A cos(omega))
 *
 * This function arises as part of the imaginary part of the
 * complex dilogarithm. It is manifestly periodic in omega.
 * It also has the symmetry omega <--> theta, which arises
 * from an underlying relation z <--> 1-z.
 */
static
int
dilog_gfunc(double A, double omega, double theta, double cos_omega, double sin_omega, gsl_sf_result * g)
{
  const double a0 = 1.0 - A * cos_omega;
  const double a1 = - A * sin_omega;
  /* const double theta = atan2(a1, a0); */
  gsl_sf_result cl2_2omega;
  gsl_sf_result cl2_2theta;
  gsl_sf_result cl2_2ot;   
  const int stat_0 = gsl_sf_clausen_e(2.0*omega, &cl2_2omega);
  const int stat_1 = gsl_sf_clausen_e(2.0*theta, &cl2_2theta);
  const int stat_2 = gsl_sf_clausen_e(2.0*(theta + omega), &cl2_2ot);
  g->val = cl2_2omega.val + cl2_2theta.val - cl2_2ot.val;
  g->err = cl2_2omega.err + cl2_2theta.err + cl2_2ot.err;
  return GSL_ERROR_SELECT_3(stat_0, stat_1, stat_2);
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/


int
gsl_sf_dilog_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x >= 0.0) {
    return dilog_xge0(x, result);
  }
  else {
    gsl_sf_result d1, d2;
    int stat_d1 = dilog_xge0( -x, &d1);
    int stat_d2 = dilog_xge0(x*x, &d2);
    result->val  = -d1.val + 0.5 * d2.val;
    result->err  =  d1.err + 0.5 * d2.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_d1, stat_d2);
  }
}


int
gsl_sf_complex_dilog_e(const double r, double theta,
                          gsl_sf_result * real_dl, gsl_sf_result * imag_dl)
{
  /* CHECK_POINTER(real_dl) */
  /* CHECK_POINTER(imag_dl) */

  if(r == 0.0) {
    real_dl->val = 0.0;
    real_dl->err = 0.0;
    imag_dl->val = 0.0;
    imag_dl->err = 0.0;
    return GSL_SUCCESS;
  }

/*
  if(theta < 0.0 || theta > 2.0*M_PI) {
    gsl_sf_angle_restrict_pos_e(&theta);
  }
*/

  /* Trap cases of real-valued argument.
   */
  if(theta == 0.0) {
    int stat_d;
    imag_dl->val = ( r > 1.0 ? -M_PI * log(r) : 0.0 );
    imag_dl->err = 2.0 * GSL_DBL_EPSILON * fabs(imag_dl->val);
    stat_d = gsl_sf_dilog_e(r, real_dl);
    return stat_d;
  }
  if(theta == M_PI) {
    int stat_d;
    imag_dl->val = 0.0;
    imag_dl->err = 0.0;
    stat_d = gsl_sf_dilog_e(-r, real_dl);
    return stat_d;
  }

  /* Trap unit circle case.
   */
  if(r == 1.0) {
    gsl_sf_result theta_restrict;
    int stat_r = gsl_sf_angle_restrict_pos_err_e(theta, &theta_restrict);
    int stat_c;
    const double term1 = theta_restrict.val*theta_restrict.val;
    const double term2 = 2.0*M_PI*fabs(theta_restrict.val);
    const double term1_err = 2.0 * fabs(theta_restrict.val * theta_restrict.err);
    const double term2_err = 2.0*M_PI*fabs(theta_restrict.err);
    real_dl->val  = M_PI*M_PI/6.0 + 0.25*(term1 - term2);
    real_dl->err  = 2.0 * GSL_DBL_EPSILON * (M_PI*M_PI/6.0 + 0.25 * (fabs(term1) + fabs(term2)));
    real_dl->err += 0.25 * (term1_err + term2_err);
    real_dl->err += 2.0 * GSL_DBL_EPSILON * fabs(real_dl->val);
    stat_c = gsl_sf_clausen_e(theta, imag_dl);
    stat_r = 0;  /* discard restrict status */
    return stat_c;
  }

  /* Generic case.
   */
  {
    int stat_dilog;
    double r_tmp, theta_tmp;
    gsl_sf_result result_re_tmp, result_im_tmp;

    /* Reduce argument to unit disk.
     */
    if(r > 1.0) {
      r_tmp     = 1.0 / r;
      theta_tmp = /* 2.0*M_PI */ - theta;
    }
    else {
      r_tmp     = r;
      theta_tmp = theta;
    }

    /* Calculate in the unit disk.
     */
    stat_dilog = dilogc_unitdisk(r_tmp, theta_tmp,
                                 &result_re_tmp, &result_im_tmp
                                 );

    /* Unwind the inversion if necessary. We calculate
     * the imaginary part explicitly if using the inversion
     * because there is no simple relationship between
     * arg(1-z) and arg(1 - 1/z), which is the "omega"
     * term in [Lewin A.2.5 (1)].
     */
    if(r > 1.0) {
      double theta_restricted = theta;
      const int stat_reduct = gsl_sf_angle_restrict_pos_e(&theta_restricted);
      const double zeta2 = M_PI*M_PI/6.0;
      const double cos_theta = cos(theta);
      const double sin_theta = sin(theta);
      const double x = r * cos_theta;
      const double y = r * sin_theta;
      const double omega = atan2(y, 1.0-x);
      const double lnr = log(r);
      const double pmt = M_PI - theta_restricted;

      gsl_sf_result Cl_a, Cl_b, Cl_c;
      const int stat_c1 = gsl_sf_clausen_e(2.0*omega, &Cl_a);
      const int stat_c2 = gsl_sf_clausen_e(2.0*theta, &Cl_b);
      const int stat_c3 = gsl_sf_clausen_e(2.0*(omega+theta), &Cl_c);
      const int stat_c  = GSL_ERROR_SELECT_3(stat_c1, stat_c2, stat_c3);
      /*
      gsl_sf_result Cl_abc;
      const int stat_c = dilog_gfunc(r, theta, cos_theta, sin_theta, &Cl_abc);
      */

      const double r1 = -result_re_tmp.val;
      const double r2 = -0.5*lnr*lnr;
      const double r3 =  0.5*pmt*pmt;
      const double r4 = -zeta2;
      const double r5 =  omega*lnr;
      real_dl->val  = r1 + r2 + r3 + r4;
      real_dl->err  = result_re_tmp.err;
      real_dl->err += GSL_DBL_EPSILON * (fabs(r1) + fabs(r2) + fabs(r3) + fabs(r4));
      real_dl->err += 2.0 * GSL_DBL_EPSILON * fabs(real_dl->val);
      imag_dl->val  = r5 + 0.5*(Cl_a.val + Cl_b.val - Cl_c.val);
      /* imag_dl->val  = r5 + 0.5*Cl_abc.val; */
      imag_dl->err  = GSL_DBL_EPSILON * fabs(r5);
      imag_dl->err += GSL_DBL_EPSILON * 0.5*(fabs(Cl_a.val) + fabs(Cl_b.val) + fabs(Cl_c.val));
      /* imag_dl->err += GSL_DBL_EPSILON * 0.5*fabs(Cl_abc.val); */
      imag_dl->err += 0.5*(Cl_a.err + Cl_b.err + Cl_c.err);
      /* imag_dl->err += 0.5*Cl_abc.err; */
      imag_dl->err += 2.0*GSL_DBL_EPSILON * fabs(imag_dl->val);
      return GSL_ERROR_SELECT_3(stat_dilog, stat_c, stat_reduct);
    }
    else {
      real_dl->val = result_re_tmp.val;
      real_dl->err = result_re_tmp.err;
      imag_dl->val = result_im_tmp.val;
      imag_dl->err = result_im_tmp.err;
      return stat_dilog;
    }
  }
}


int
gsl_sf_complex_spence_e(
  const double A,
  const double omega,
  gsl_sf_result * real_dl,
  gsl_sf_result * imag_dl
  )
{
  const double cos_omega = cos(omega);
  const double sin_omega = sin(omega);
  const double ax = 1.0 - A * cos_omega;
  const double ay = - A * sin_omega;
  const double r = sqrt(ax*ax + ay*ay);
  const double theta = atan2(ay, ax);

  const double imag_term_disc = omega * log(r);
  gsl_sf_result imag_term_periodic;
  const int stat_g = dilog_gfunc(A, omega, theta, cos_omega, sin_omega, &imag_term_periodic);
  imag_dl->val = 0.5 * (imag_term_disc + imag_term_periodic.val);
  imag_dl->err = 0.5 * imag_term_periodic.err;
  imag_dl->err += GSL_DBL_EPSILON * fabs(imag_dl->val);

  return 0; /* GSL_ERROR_SELECT_(stat_g); */
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_dilog(const double x)
{
  EVAL_RESULT(gsl_sf_dilog_e(x, &result));
}

