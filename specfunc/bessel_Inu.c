/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_temme.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#define locMin(a,b)  ((a) < (b) ? (a) : (b))
#define locMax(a,b)  ((a) > (b) ? (a) : (b))


/* Evaluate continued fraction CF1 to obtain I'/I.
 */
static
int
bessel_I_CF1(const double nu, const double x, double * result)
{
  const int max_iter = 5000;
  int i = 0;
  double x_inv = 1.0/x;
  double r = nu * x_inv;
  double b = 2.0*nu*x_inv;
  double d = 0.0;
  double c;
  r = locMax(r, 1.0e-100);
  c = r;
  while(i < max_iter) {
    double del;
    b  += 2.0*x_inv;
    d   = 1.0/(b+d);
    c   = b + 1.0/c;
    del = c*d;
    r  *= del;
    if (fabs(del-1.0) < GSL_MACH_EPS) break;
    ++i;
  }
  
  *result = r;
  if(i == max_iter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}


/* Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
 * to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
 */
static
int
bessel_K_CF2(const double nu, const double x, double * K_nu, double * K_nup1)
{
  const int maxiter = 10000;

  int i = 1;
  double bi = 2.0*(1.0 + x);
  double di = 1.0/bi;
  double delhi = di;
  double hi    = di;

  double qi   = 0.0;
  double qip1 = 1.0;

  double ai = -(0.25 - nu*nu);
  double a1 = ai;
  double ci = -ai;
  double Qi = -ai;

  double s = 1.0 + Qi*delhi;

  for(i=2; i<=maxiter; i++) {
    double dels;
    double tmp;
    ai -= 2.0*(i-1);
    ci  = -ai*ci/i;
    tmp  = (qi - bi*qip1)/ai;
    qi   = qip1;
    qip1 = tmp;
    Qi += ci*qip1;
    bi += 2.0;
    di  = 1.0/(bi + ai*di);
    delhi = (bi*di - 1.0) * delhi;
    hi += delhi;
    dels = Qi*delhi;
    s += dels;
    if(fabs(dels/s) < GSL_MACH_EPS) break;
  }
  
  hi *= -a1;
  
  *K_nu   = sqrt(M_PI/(2.0*x)) / s;
  *K_nup1 = *K_nu * (nu + x + 0.5 - hi)/x;
  if(i == maxiter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}


/*
 * Assumes x > 0  and  nu_min >= 0.
 */
static
int
bessel_InuKnu_steed_scaled(const double nu_min, const double x, const int kmax,
                           double * result_I_array,  double * result_K_array,
                           double * result_Ip_array, double * result_Kp_array)
{
  int n;
  int k;
  int N     = (int)(nu_min + 0.5);
  double mu = nu_min - N;
  double I_small = 100.0*DBL_MIN;
  double I_ratio_numax;
  double I_ratio_numin;
  double I_scale_numin;
  double I_min, Ip_min;
  double I_min_save;
  double K_mu, K_mup1, K_min, K_minp1;
  double Kp_min;
  double K_max, Kp_max;

  /* Evaluate I'/I at nu_max.
   */
  bessel_I_CF1(nu_min + kmax, x, &I_ratio_numax);

  /* Do backward recurrence for I and I' from nu_max to nu_min.
   */
  gsl_sf_bessel_I_recur(nu_min, x, kmax,
                        I_small, I_small*I_ratio_numax,
		        &I_min, &Ip_min,
		        result_I_array, result_Ip_array
		        );
  I_ratio_numin = Ip_min/I_min;

  /* Evaluate K_mu and K_{mu+1}.
   */
  if(x < 2.0) {
    double Kp_mu;
    gsl_sf_bessel_K_scaled_temme(mu, x, &K_mu, &K_mup1, &Kp_mu);
  }
  else {
    bessel_K_CF2(mu, x, &K_mu, &K_mup1);
  }

  /* Do forward recurrence to obtain K and K' at nu_min.
   */
  K_min = K_mu;
  K_minp1 = K_mup1;
  for(n=1; n<=N; n++) {
    double nu = mu + n;
    double K_minp1_save = K_minp1;
    K_minp1 = K_min + 2.0*nu/x*K_minp1;
    K_min   = K_minp1_save;
  }

  /* Use Wronskian relation to obtain correctly
   * normalized values for K,K' I,I' at nu_min,
   * and determine scaling factor that must be
   * applied to the I,I' values.
   */
  Kp_min = nu_min/x * K_min - K_minp1;
  I_min_save = I_min;
  I_min  = 1.0/(x*(I_ratio_numin*K_min - Kp_min));
  Ip_min = I_ratio_numin*I_min;
  I_scale_numin = I_min / I_min_save;

  /* Apply scaling factor to I,I'.
   */
  if(result_I_array != (double *)0) {
    for(k=0; k<=kmax; k++) {
      result_I_array[k] *= I_scale_numin;
    }
  }
  if(result_Ip_array != (double *)0) {
    for(k=0; k<=kmax; k++) {
      result_Ip_array[k] *= I_scale_numin;
    }
  }

  /* Do forward recurrence to obtain K,K' values.
   */
  gsl_sf_bessel_K_recur(nu_min, x, kmax,
                        K_min, Kp_min,
			&K_max, &Kp_max,
			result_K_array, result_Kp_array
			);
 
}

/* Assumes we can apply the uniform asymptotic result
 * to determine starting values for recursion.
 */
static
int
bessel_InuKnu_scaled_asymp(const double nu_min, const double x, const int kmax,
                           double * result_I_array,  double * result_K_array,
                           double * result_Ip_array, double * result_Kp_array)
{
  double nu_max = nu_min + kmax;
  double K_min, K_minp1;
  double I_max, I_maxp1;
  int stat_I   = gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nu_max,     x, &I_max);
  int stat_Ip1 = gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nu_max+1.0, x, &I_maxp1);
  int stat_K   = gsl_sf_bessel_Knu_scaled_asymp_unif_impl(nu_min,     x, &K_min);
  int stat_Kp1 = gsl_sf_bessel_Knu_scaled_asymp_unif_impl(nu_min+1.0, x, &K_minp1);
  double Ip_max = I_maxp1 + nu_max/x * I_max;
  double I_min, Ip_min;
  double K_max, Kp_max;
  double Kp_min = -K_minp1 + nu_min/x * K_min;

  gsl_sf_bessel_I_recur(nu_min, x, kmax,
                        I_max, Ip_max,
	                &I_min, &Ip_min,
	                result_I_array, result_Ip_array
	                );

  gsl_sf_bessel_K_recur(nu_min, x, kmax,
                        K_min, Kp_min,
			&K_max, &Kp_max,
                        result_K_array, result_Kp_array
			);
}



int
gsl_sf_bessel_InuKnu_scaled_impl(double nu, double x, int kmax,
                                 double * result_I,  double * result_K,
                                 double * result_Ip, double * result_Kp)
{
  if(x < 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    int status = 0;
    int k;
    if(result_I != (double *)0) {
      for(k=0; k<=kmax; k++) result_I[k] = 0.0;
      if(nu == 0.0) result_I[k] = 1.0;
    }
    if(result_Ip != (double *)0) {
      for(k=0; k<=kmax; k++) result_Ip[k] = 0.0;
      if(nu < 1.0) status += 1;
      else if(nu == 1.0) result_Ip[0] = 0.5;
    }
    if(result_K != (double *)0) {
      for(k=0; k<=kmax; k++) result_K[k] = 0.0;
      status += 1;
    }
    if(result_Kp != (double *)0) {
      for(k=0; k<=kmax; k++) result_Kp[k] = 0.0;
      status += 1;
    }
    if(status) {
      return GSL_EDOM;
    }
    else {
      return GSL_SUCCESS;
    }
  }
  else if(locMin( 0.29/(nu*nu), 0.5/(nu*nu + x*x) ) < GSL_ROOT3_MACH_EPS) {
    return bessel_InuKnu_scaled_asymp(nu, x, kmax,
                                      result_I,  result_K,
                                      result_Ip, result_Kp
			              );
  }
  else {
    return bessel_InuKnu_steed_scaled(nu, x, kmax,
                                      result_I,  result_K,
                                      result_Ip, result_Kp
			              );
  }
}

int
gsl_sf_bessel_Inu_impl(double nu, double x, double * result)
{
  if(x < 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else if(x*x < 10.0*(nu+1.0)*GSL_ROOT5_MACH_EPS) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 4, &b);
    *result = ex * b;
    return stat;
  }
  else if(x*x < 10.0*(nu+1.0)) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 100, result);
    *result = ex * b;
    return stat;
  }
  else if(locMin( 0.29/(nu*nu), 0.5/(nu*nu + x*x) ) < GSL_ROOT3_MACH_EPS) {
    return gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nu, x, result);
  }
  else {
    /* FIXME */
  }
}


int
gsl_sf_bessel_Inu_e(double nu, double x, double * result)
{
  int status = gsl_sf_bessel_Inu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Inu_e", status);
  }
  return status;
}


double
gsl_sf_bessel_Inu(double nu, double x)
{
  double y;
  int status = gsl_sf_bessel_Inu_impl(nu, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Inu", status);
  }
  return y;
}
