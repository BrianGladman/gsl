/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_temme.h"
#include "gsl_sf_bessel.h"


static
int
bessel_Y_temme(const double nu, const double x,
               double * Y_nu, double * Y_nup1, double * Yp_nu)
{
  const int max_iter = 15000;
  
  const double half_x = 0.5 * x;
  const double ln_half_x = log(half_x);
  const double half_x_nu = exp(nu*ln_half_x);
  const double pi_nu   = M_PI * nu;
  const double alpha   = pi_nu / 2.0;
  const double sigma   = -nu * ln_half_x;
  const double sinrat  = (fabs(pi_nu) < GSL_MACH_EPS ? 1.0 : pi_nu/sin(pi_nu));
  const double sinhrat = (fabs(sigma) < GSL_MACH_EPS ? 1.0 : sinh(sigma)/sigma);
  const double sinhalf = (fabs(alpha) < GSL_MACH_EPS ? 1.0 : sin(alpha)/alpha);
  const double sin_sqr = nu*M_PI*M_PI*0.5 * sinhalf*sinhalf;
  
  double sum0, sum1;
  double fk, pk, qk, hk, ck;
  int k = 0;

  double g_1pnu, g_1mnu, g1, g2;
  gsl_sf_temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2);
  
  fk = 2.0/M_PI * sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
  pk = 1.0/M_PI /half_x_nu * g_1pnu;
  qk = 1.0/M_PI *half_x_nu * g_1mnu;
  hk = pk;
  ck = 1.0;

  sum0 = fk + sin_sqr * qk;
  sum1 = pk;

  while(k < max_iter) {
    double del0;
    double del1;
    double gk;
    k++;
    fk  = (k*fk + pk + qk)/(k*k-nu*nu);
    ck *= -half_x*half_x/k;
    pk /= (k - nu);
    qk /= (k + nu);
    gk  = fk + sin_sqr * qk;
    hk  = -k*gk + pk; 
    del0 = ck * gk;
    del1 = ck * hk;
    sum0 += del0;
    sum1 += del1;
    if (fabs(del0) < (1.0 + fabs(sum0))*GSL_MACH_EPS) break;
  }

  *Y_nu   = -sum0;
  *Y_nup1 = -sum1 * 2.0/x;
  *Yp_nu  = nu/x * *Y_nu - *Y_nup1;
  
/* rjmu=w/(rymup-f*rymu); Equation (6.7.13). */

  if(k == max_iter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}


static
int
bessel_J_CF1(const double nu, const double x, double * result, double * sign)
{
  const int max_iter = 5000;
  const double SMALL = 1.0e-100;

  int i = 1;
  int isign = 1;
  double x_inv = 1.0/x;
  double r = nu * x_inv;
  double b = 2.0*nu*x_inv;
  double d = 0.0;
  double c;
  r = locMax(r, SMALL);
  c = r;
  while(i < max_iter) {
    double del;
    b  += 2.0*x_inv;
    d   = b - d;
    if(fabs(d) < SMALL) d = SMALL;
    c   = b - 1.0/c;
    if(fabs(c) < SMALL) c = SMALL;
    d = 1.0/d;
    del = c*d;
    r *= del;
    if(d < 0.0) isign = -isign;
    if(fabs(del-1.0) < GSL_MACH_EPS) break;
  }
  
  *result = r;
  *sign = isign;
  if(i == max_iter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}


static
int
bessel_J_recur(const double nu_min, const double x, const int kmax,
               const double J_start, const double Jp_start,
	       double * J_end, double * Jp_end,
	       double * J_array, double * Jp_array)
{
  double x_inv  = 1.0/x;
  double nu_max = nu_min + kmax;
  double nu     = nu_max;
  double J_nu  = J_start;
  double Jp_nu = Jp_start;
  int k;

  if(J_array  != (double *)0)  J_array[kmax] = J_start;
  if(Jp_array != (double *)0) Jp_array[kmax] = Jp_start;
  for(k=kmax-1; k>=0; k--) {
    double nu_x_inv = nu*x_inv;
    double J_nu_save = J_nu;
    J_nu  = nu_x_inv*J_nu + Jp_nu;
    Jp_nu = (nu_x_inv - x_inv)*J_nu - J_nu_save;
    if(J_array  != (double *)0)  J_array[k] = J_nu;
    if(Jp_array != (double *)0) Jp_array[k] = Jp_nu;
    nu -= 1.0;
  }
  *J_end  = J_nu;
  *Jp_end = Jp_nu;
  return GSL_SUCCESS;
}


static
int
bessel_Y_CF2(const double nu, const double x,
             double * Y_nu, double * Y_nup1, double * Yp_nu)
{
  const int max_iter = 10000;
  const double SMALL = 1.0e-100;

  int i = 1;
  double gamma;

  double x_inv = 1.0/x;
  double a = 0.25 - nu*nu;
  double p = -0.5*x_inv;
  double q = 1.0;
  double br = 2.0*x;
  double bi = 2.0;
  double fact = a*x_inv/(p*p + q*q);
  double cr = br + q*fact;
  double ci = bi + p*fact;
  double den = br*br + bi*bi;
  double dr = br/den;
  double di = -bi/den;
  double dlr = cr*dr - ci*di;
  double dli = cr*di + ci*dr;
  double temp = p*dlr - q*dli;
  q = p*dli + q*dlr;
  p = temp;
  for (i=2; i<=max_iter; i++) {
    a  += 2*(i-1);
    bi += 2.0;
    dr = a*dr + br;
    di = a*di + bi;
    if(fabs(dr)+fabs(di) < SMALL) dr = SMALL;
    fact = a/(cr*cr+ci*ci);
    cr = br + cr*fact;
    ci = bi - ci*fact;
    if(fabs(cr)+fabs(ci) < SMALL) cr = SMALL;
    den = dr*dr + di*di;
    dr /= den;
    di /= -den;
    dlr = cr*dr - ci*di;
    dli = cr*di + ci*dr;
    temp = p*dlr - q*dli;
    q = p*dli + q*dlr;
    p = temp;
    if(fabs(dlr-1.0)+fabs(dli) < GSL_MACH_EPS) break;
  }

/*
  gamma = (p-f)/q;
  rjmu = sqrt(w/((p-f)*gamma + q));
  rjmu = SIGN(rjmu, rjl);
  rymu = rjmu*gamma;
  rymup = rymu*(p+q/gamma);
  ry1 = xmu*xi*rymu - rymup;
  */

  if(i == max_iter) {
    return GSL_EMAXITER;
  }
  else {
    return GSL_SUCCESS;
  }
}


static
int
bessel_Y_recur(const double nu_min, const double x, const int kmax,
               const double Y_min, const double Y_minp1,
               double * result_Y_array, double * result_Yp_array)
{
  int k;

  if(result_Y_array != (double *)0) {
    result_Y_array[0] = Y_min;
    if(kmax > 0) result_Y_array[1] = Y_minp1;
  }
  if(result_Yp_array != (double *)0) {
    result_Yp_array[0] = -Y_minp1 + nu_min/x * Y_min;
    if(kmax > 0) result_Yp_array[1] = -(nu_min + 1.0)/x * Y_minp1 - Y_min;
  }

  if(kmax > 0) {
    double Y_k;
    double Y_km1 = Y_minp1;
    double Y_km2 = Y_min;
    for(k=2; k<=kmax; k++) {
      double nu   = nu_min + k - 1.0;
      double nuox = nu/x;
      Y_k = -Y_km2 + 2.0*nuox * Y_km1;
      if(result_Y_array  != (double *)0)  result_Y_array[k] = Y_k;
      if(result_Yp_array != (double *)0) result_Yp_array[k] = -nuox*Y_k+Y_km1;
      Y_km2 = Y_km1;
      Y_km1 = Y_k;
    }
  }
}


static
int
bessel_JnuYnu_steed(const double nu_min, const double x, const int kmax,
                    double * result_J_array,  double * result_Y_array,
                    double * result_Jp_array, double * result_Yp_array)
{
  int n;
  int k;
  int N     = (int)(nu_min + 0.5);
  double mu = nu_min - N;
  double J_small = 100.0*DBL_MIN;
  double J_ratio_numax;
  double J_ratio_numin;
  double J_sign_numax;
  double J_scale_numin;
  double J_min, Jp_min;
  double J_min_save;
  double Y_mu, Y_mup1, Yp_mu, Y_min, Y_minp1;
  double Yp_min;

  /* Evaluate J'/J at nu_max.
   */
  bessel_J_CF1(nu_min + kmax, x, &J_ratio_numax, &J_sign_numax);

  /* Do backward recurrence for J and J' from nu_max to nu_min.
   */
  bessel_J_recur(nu_min, x, kmax,
                 J_small, J_small*J_ratio_numax,
		 &J_min, &Jp_min,
		 result_J_array, result_Jp_array
		 );
  J_ratio_numin = Jp_min/J_min;

  /* Evaluate Y_mu and Y_{mu+1}.
   */
  if(x < 2.0) {
    bessel_Y_temme(mu, x, &Y_mu, &Y_mup1, &Yp_mu);
  }
  else {
    bessel_Y_CF2(mu, x, &Y_mu, &Y_mup1, &Yp_mu);
  }

  /* Do forward recurrence to obtain K and K' at nu_min.
   */
  Y_min = Y_mu;
  Y_minp1 = Y_mup1;
  for(n=1; n<=N; n++) {
    double nu = mu + n;
    double Y_minp1_save = Y_minp1;
    Y_minp1 = Y_min + 2.0*nu/x*Y_minp1;
    Y_min   = Y_minp1_save;
  }

  /* Use Wronskian relation to obtain correctly
   * normalized values for Y,Y' J,J' at nu_min,
   * and determine scaling factor that must be
   * applied to the J,J' values.
   */
  Yp_min = nu_min/x * Y_min - Y_minp1;
  J_min_save = J_min;
  J_min  = 1.0/(x*(J_ratio_numin*J_min - Jp_min));
  Jp_min = J_ratio_numin*J_min;
  J_scale_numin = J_min / J_min_save;

  /* Apply scaling factor to J,J'.
   */
  if(result_J_array != (double *)0) {
    for(k=0; k<=kmax; k++) {
      result_J_array[k] *= J_scale_numin;
    }
  }
  if(result_Jp_array != (double *)0) {
    for(k=0; k<=kmax; k++) {
      result_Jp_array[k] *= J_scale_numin;
    }
  }

  /* Do forward recurrence to obtain Y,Y' values.
   */
  bessel_Y_recur(nu_min, x, kmax, Y_min, Y_minp1, result_Y_array, result_Yp_array);
 
}


static
int
bessel_JnuYnu_asymp(const double nu_min, const double x, const int kmax,
                    double * result_J_array,  double * result_Y_array,
                    double * result_Jp_array, double * result_Yp_array)
{
  double nu_max = nu_min + kmax;
  double Y_min, Y_minp1;
  double J_max, J_maxp1;
  int stat_J   = gsl_sf_bessel_Jnu_asymp_Olver_impl(nu_max,     x, &J_max);
  int stat_Jp1 = gsl_sf_bessel_Jnu_asymp_Olver_impl(nu_max+1.0, x, &J_maxp1);
  int stat_Y   = gsl_sf_bessel_Ynu_asymp_Olver_impl(nu_min,     x, &Y_min);
  int stat_Yp1 = gsl_sf_bessel_Ynu_asymp_Olver_impl(nu_min+1.0, x, &Y_minp1);
  double Jp_max = J_maxp1 + nu_max/x * J_max;
  double J_min, Jp_min;

  bessel_J_recur(nu_min, x, kmax,
                 J_max, Jp_max,
	         &J_min, &Jp_min,
	         result_J_array, result_Jp_array
	         );

  bessel_Y_recur(nu_min, x, kmax,
                 Y_min, Y_minp1,
                 result_Y_array, result_Yp_array);
}


int
gsl_sf_bessel_JnuYnu_impl(const double nu, const double x, const int kmax,
                          double * result_J,  double * result_Y,
                          double * result_Jp, double * result_Yp)
{
  if(x < 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    int status = 0;
    int k;
    if(result_J != (double *)0) {
      for(k=0; k<=kmax; k++) result_J[k] = 0.0;
      if(nu == 0.0) result_J[k] = 1.0;
    }
    if(result_Jp != (double *)0) {
      for(k=0; k<=kmax; k++) result_Jp[k] = 0.0;
      if(nu < 1.0) status += 1;
      else if(nu == 1.0) result_Jp[0] = 0.5;
    }
    if(result_Y != (double *)0) {
      for(k=0; k<=kmax; k++) result_Y[k] = 0.0;
      status += 1;
    }
    if(result_Yp != (double *)0) {
      for(k=0; k<=kmax; k++) result_Yp[k] = 0.0;
      status += 1;
    }
    if(status) {
      return GSL_EDOM;
    }
    else {
      return GSL_SUCCESS;
    }
  }
  else if(nu > 30.0) {
    return bessel_JnuYnu_asymp(nu, x, kmax,
                               result_J,  result_Y,
                               result_Jp, result_Yp
			       );
  }
  else {
    return bessel_JnuYnu_steed(nu, x, kmax,
                               result_J,  result_Y,
                               result_Jp, result_Yp
			       );
  }
}


int
gsl_sf_bessel_Jnu_impl(double nu, double x, double * result)
{
  const double nu_cut = 20.0;

  if(x < 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else if(x*x < 10.0*(nu+1.0)*GSL_ROOT5_MACH_EPS) {
    return gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, -1, 4, result);
  }
  else if(x*x < 10.0*(nu+1.0)) {
    return gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, -1, 100, result);
  }
  else if(nu > nu_cut) {
    return gsl_sf_bessel_Jnu_asymp_Olver_impl(nu, x, result);
  }
  else {
    /* Evaluate at large enough nu and apply backward recurrence */

    int n;
    int steps = ceil(nu_cut - nu) + 1;
    double Jnp1;
    double Jn;
    double Jnm1;
    
    gsl_sf_bessel_Jnu_asymp_Olver_impl(nu + steps + 1.0, x, &Jnp1);
    gsl_sf_bessel_Jnu_asymp_Olver_impl(nu + steps      , x, &Jn);
    
    for(n=steps; n>0; n--) {
      Jnm1 = 2.0*(nu+n)/x * Jn - Jnp1;
      Jnp1 = Jn;
      Jn   = Jnm1;
    }
    
    *result = Jnm1;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_bessel_Jnu_e(double nu, double x, double * result)
{
  int status = gsl_sf_bessel_Jnu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Jnu_e", status);
  }
  return status;
}


double
gsl_sf_bessel_Jnu(double nu, double x)
{
  double y;
  int status = gsl_sf_bessel_Jnu_impl(nu, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Jnu", status);
  }
  return y;
}
