/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* l >= 0;  0 <= x <= 1
 *
 * if harvest != 0: 
 *   size of harvest[] = l-m+1
 *   harvest[0] = P_m^m(x)  ...  harvest[l-m] = P_l^m(x)
 */
int gsl_sf_legendre_Pl_impl(const int l,
                            const double x,
                            double * result, double * harvest
			    )
{ 
  if(l < 0 || x < 0. || x > 1.) {
    return GSL_EDOM;
  }
  else if(l == 0) {
    *result = 1.;
    if(harvest != 0) harvest[0] = 1.;
    return GSL_SUCCESS;
  }
  else if(l == 1) {
    *result = x;
    if(harvest != 0) {
      harvest[0] = 1.;
      harvest[1] = x;
    }
    return GSL_SUCCESS;
  }
  else if(l == 2) {
    double P2 = 0.5*(3.*x*x-1.);
    *result = P2;
    if(harvest != 0) {
      harvest[0] = 1.;
      harvest[1] = x;
      harvest[2] = P2;
    }
    return GSL_SUCCESS;
  }
  else {
    int ell;
    double pmm   = 1.;             /* P_m^m(x) */
    double pmmp1 = x * pmm;        /* P_{m+1}^m(x) */
    double p_ell = pmmp1;

    /* Compute P_l, l > 1, upward recurrence on l. */
    if(harvest != 0) {
      harvest[0] = pmm;
      harvest[1] = pmmp1;
      for(ell=2; ell <= l; ell++){
        p_ell = (x*(2*ell-1)*pmmp1 - (ell-1)*pmm) / ell;
        pmm = pmmp1;
        pmmp1 = p_ell;
        harvest[ell] = p_ell;
      }
    }
    else {
      for(ell=2; ell <= l; ell++){
        p_ell = (x*(2*ell-1)*pmmp1 - (ell-1)*pmm) / ell;
        pmm = pmmp1;
        pmmp1 = p_ell;
      }
    }

    *result = p_ell;
    return GSL_SUCCESS;
  }
}

/* l >= m >= 0;  0 <= x <= 1
 *
 * if harvest != 0: 
 *   size of harvest[] = l-m+1
 *   harvest[0] = P_m^m(x)  ...  harvest[l-m] = P_l^m(x)
 */
int gsl_sf_legendre_Plm_impl(const int l, const int m,
                             const double one_m_x, const double one_p_x,
                             double * result, double * harvest
			     )
{
  double x = 0.5 * (one_p_x - one_m_x);
  double pmm;                 /* P_m^m(x) */
  double pmmp1;               /* P_{m+1}^m(x) */

  /* If l is large and m is large, then we have to worry
   * about overflow. Calculate an approximate exponent which
   * measures the normalization of this thing.
   */
  double dif = l-m;
  double sum = l+m;
  double exp_check = 0.5 * log(2.*l+1.) 
                     + 0.5 * dif * (log(dif)-1.)
                     - 0.5 * sum * (log(sum)-1.);

  /* check args */
  if(m < 0 || l < m || x < 0. || x > 1.) {
    return GSL_EDOM;
  }
  
  /* Bail out if it looks like overflow. */  
  if(exp_check < GSL_LOG_DBL_MIN + 10.){
    return GSL_EOVRFLW;
  }

  /* Calculate P_m^m from the analytic result:
   *          P_0^0(x) = 1
   *          P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2) , m > 0
   */
  pmm = 1.;
  if(m > 0){
    int i;
    double circ = sqrt(one_m_x)*sqrt(one_p_x);
    double fact = 1.;
    for(i=1; i<=m; i++){
      pmm  *= -fact * circ;
      fact += 2.;
    }
  }

  /* and calculate P_{m+1}^m */
  pmmp1 = x * (2*m + 1) * pmm;
  
  if(harvest != 0) harvest[0] = pmm;
  
  if(l == m){
    *result = pmm;
  }
  else{
    int ell;
    double p_ell = pmmp1;
    
    /* Compute P_l^m, l > m+1, upward recurrence on l. */
    if(harvest != 0) {
      harvest[1] = pmmp1;
      for(ell=m+2; ell <= l; ell++){
        p_ell = (x*(2*ell-1)*pmmp1 - (ell+m-1)*pmm) / (ell-m);
        pmm = pmmp1;
        pmmp1 = p_ell;
        harvest[ell-m] = p_ell;
      }
    }
    else {
      for(ell=m+2; ell <= l; ell++){
        p_ell = (x*(2*ell-1)*pmmp1 - (ell+m-1)*pmm) / (ell-m);
        pmm = pmmp1;
        pmmp1 = p_ell;
      }
    }

    *result = p_ell;
    return GSL_SUCCESS;
  }
}

/* l >= m >= 0; 0 <= x <= 1
 *
 * if harvest != 0: 
 *   size of harvest[] = l-m+1
 *   harvest[0] = sphP_m^m(x)  ...  harvest[l-m] = sphP_l^m(x)
 */
int gsl_sf_legendre_sphPlm_impl(const int l, int m,
                                const double one_m_x, const double one_p_x,
                                double * result, double * harvest
				)
{
  int i;
  double x = 0.5 * (one_p_x - one_m_x);

  /* Starting value for recursion, Y_m^m(x).
   * We include part of the normalization factor here.
   */
  double ymm = sqrt((2.*(double)l+1.) / (4.*M_PI));

  /* check args */
  if(m < 0 || l < m || x < 0. || x > 1.) {
    return GSL_EDOM;
  }

  /* If m > 0, then calculate Y_m^m from the analytic result.
   * If m==0, then we don't have to do anything; ymm is ready to go.
   * 
   *          Y_m^m(x) = sqrt(1/(2m)!) (-1)^m (2m-1)!! (1-x^2)^(m/2)
   */
  if(m > 0){
    double circ = sqrt(one_m_x*one_p_x);
    double fact1 = 1.;
    double fact2 = 1. / M_SQRT2;
    for(i=1; i<=m; i++){
      ymm   *= -fact1 * fact2 * circ;
      fact1 += 2.;
      fact2  = 1. / sqrt(fact1 * (fact1 + 1.));
    }
  }

  if(harvest != 0) harvest[0] = ymm;

  if(l == m){
    *result = ymm;
    return GSL_SUCCESS;
  }
  else{
    int ell;
    double ymmp1 = x * sqrt(2.*(double)m + 1) * ymm;  /* Y_{m+1}^m. */
    double y_ell = ymmp1;

    /* Compute Y_l^m, l > m+1, downward recursion on l. */
    if(harvest != 0) {
      harvest[1] = ymmp1;
      for(ell=m+2; ell <= l; ell++){
        double factor1 = sqrt((double)(ell-m) / (double)(ell+m));
        double factor2 = factor1 * sqrt((double)(ell-m-1.) / (double)(ell+m-1.));
        y_ell = (x*(2.*ell-1)*ymmp1*factor1 - (ell+m-1.)*ymm*factor2) / (ell-m);
        ymm   = ymmp1;
        ymmp1 = y_ell;
        harvest[ell-m] = y_ell;
      }
    }
    else {
      for(ell=m+2; ell <= l; ell++){
        double factor1 = sqrt((double)(ell-m) / (double)(ell+m));
        double factor2 = factor1 * sqrt((double)(ell-m-1.) / (double)(ell+m-1.));
        y_ell = (x*(2.*ell-1)*ymmp1*factor1 - (ell+m-1.)*ymm*factor2) / (ell-m);
        ymm   = ymmp1;
        ymmp1 = y_ell;
      }
    }
    
    *result = y_ell;
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_legendre_Pl_e(const int l, const double x, double * result)
{
  int status = gsl_sf_legendre_Pl_impl(l, x, result, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Pl_e", status);
  }
  return status;
}

int gsl_sf_legendre_Plm_e(const int l, const int m, const double x, double * result)
{
  int status = gsl_sf_legendre_Plm_impl(l, m, 1.-x, 1.+x, result, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Plm_e", status);
  }
  return status;
}

int gsl_sf_legendre_sphPlm_e(const int l, const int m, const double x, double * result)
{
  int status = gsl_sf_legendre_sphPlm_impl(l, m, 1.-x, 1.+x, result, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_sphPlm_e", status);
  }
  return status;
}

int gsl_sf_legendre_Pl_array_e(const int l, const double x, double * result_array)
{
  double y;
  int status = gsl_sf_legendre_Pl_impl(l, x, &y, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Pl_array_e", status);
  }
  return status;
}

int gsl_sf_legendre_Plm_array_e(const int lmax, const int m, const double x, double * result_array)
{
  double y;
  int status = gsl_sf_legendre_Plm_impl(lmax, m, 1.-x, 1.+x, &y, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Plm_array_e", status);
  }
  return status;
}

int gsl_sf_legendre_sphPlm_array_e(const int lmax, const int m, const double x, double * result_array)
{
  double y;
  int status = gsl_sf_legendre_sphPlm_impl(lmax, m, 1.-x, 1.+x, &y, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_sphPlm_array_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_legendre_array_size(const int lmax, const int m)
{
  return lmax-m+1;
}

double gsl_sf_legendre_Pl(const int l, const double x)
{
  double y;
  int status = gsl_sf_legendre_Pl_impl(l, x, &y, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Pl", status);
  }
  return y;
}

double gsl_sf_legendre_Plm(const int l, const int m, const double x)
{
  double y;
  int status = gsl_sf_legendre_Plm_impl(l, m, 1.-x, 1.+x, &y, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Plm", status);
  }
  return y;
}

double gsl_sf_legendre_sphPlm(const int l, const int m, const double x)
{
  double y;
  int status = gsl_sf_legendre_sphPlm_impl(l, m, 1.-x, 1.+x, &y, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_sphPlm", status);
  }
  return y;
}
