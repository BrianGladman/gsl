/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"

/* FIXME: These are not good for |x| very near 1, due to
 * the 1-x*x problem.
 */

/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* l >= m >= 0; |x| <= 1; l >= 0 */
int gsl_sf_legendre_Plm_impl(const int l, const int m, const double x, double * result, double * harvest)
{
  int i;
  double pmm;  /* Starting value for recursion, P_m^m(x) */

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
  if(m < 0 || m > l || x < 0. || l < 0) {
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
    double circ = sqrt((1.-x)*(1.+x));
    double fact = 1.;
    for(i=1; i<=m; i++){
      pmm  *= -fact * circ;
      fact += 2.;
    }
  }

  if(harvest != 0) harvest[0] = pmm;
  
  if(l == m){
    *result = pmm;
    return GSL_SUCCESS;
  }
  else{
    double pmmp1 = x * (2*m + 1) * pmm;    /* P_{m+1}^m. */

    if(harvest != 0) harvest[1] = pmmp1;
    
    if(l == (m+1)){
      *result = pmmp1;  
      return GSL_SUCCESS;
    }
    else {
      /* Otherwise compute P_l^m, l > m+1, recursively. */
      double p_ell;
      int ell;

      if(harvest != 0) {
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
}

/* l >= |m| >= 0; |x| <= 1; l >= 0 */
int gsl_sf_legendre_sphPlm_impl(const int l, int m, const double x, double * result, double * harvest)
{
  int i;
  double sgn_factor;

  /* Starting value for recursion, Y_m^m(x).
   * We include part of the normalization factor here.
   */
  double ymm = sqrt((2.*(double)l+1.) / (4.*M_PI));

  /* check args */
  if(abs(m) > l || l < 0 || x < 0.) {
    return GSL_EDOM;
  }
  
  /* strip sign of m */
  if(m < 0) {
    m = -m;
    sgn_factor = (GSL_IS_ODD(m) ? -1. : 1.);  /* (-1)^m */
  }
  else {
    sgn_factor = 1.;
  }

  /* If m > 0, then calculate Y_m^m from the analytic result.
   * If m==0, then we don't have to do anything; ymm is ready to go.
   * 
   *          Y_m^m(x) = sqrt(1/(2m)!) (-1)^m (2m-1)!! (1-x^2)^(m/2)
   */
  if(m > 0){
    double circ = sqrt((1.-x)*(1.+x));
    double fact1 = 1.;
    double fact2 = 1. / sqrt(2.);
    for(i=1; i<=m; i++){
      ymm   *= -fact1 * fact2 * circ;
      fact1 += 2.;
      fact2  = 1. / sqrt(fact1 * (fact1 + 1.));
    }
  }

  if(harvest != 0) harvest[0] = sgn_factor * ymm;

  if(l == m){
    *result = sgn_factor * ymm;
    return GSL_SUCCESS;
  }
  else{
    double ymmp1 = x * sqrt(2.*(double)m + 1) * ymm;  /* Y_{m+1}^m. */

    if(harvest != 0) harvest[1] = sgn_factor * ymmp1;

    if(l == (m+1)){
      *result = sgn_factor * ymmp1;
      return GSL_SUCCESS;
    }
    else {
      /* Otherwise compute Y_l^m, l > m+1, recursively. */
      double y_ell;
      int ell;
      
      if(harvest != 0) {
        for(ell=m+2; ell <= l; ell++){
	  double factor1 = sqrt((double)(ell-m) / (double)(ell+m));
	  double factor2 = factor1 * sqrt((double)(ell-m-1.) / (double)(ell+m-1.));
	  y_ell = (x*(2.*ell-1)*ymmp1*factor1 - (ell+m-1.)*ymm*factor2) / (ell-m);
	  ymm   = ymmp1;
	  ymmp1 = y_ell;
	  harvest[ell-m] = sgn_factor * y_ell;
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
      
      *result = sgn_factor * y_ell;
      return GSL_SUCCESS;
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_legendre_Plm_e(const int l, const int m, const double x, double * result)
{
  int status = gsl_sf_legendre_Plm_impl(l, m, x, result, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Plm_e", status);
  }
  return status;
}

int gsl_sf_legendre_sphPlm_e(const int l, const int m, const double x, double * result)
{
  int status = gsl_sf_legendre_sphPlm_impl(l, m, x, result, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_sphPlm_e", status);
  }
  return status;
}

int gsl_sf_legendre_Plm_array_e(const int lmax, const int m, const double x, double * result_array)
{
  double y;
  int status = gsl_sf_legendre_Plm_impl(lmax, m, x, &y, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Plm_array_e", status);
  }
  return status;
}

int gsl_sf_legendre_sphPlm_array_e(const int lmax, const int m, const double x, double * result_array)
{
  double y;
  int status = gsl_sf_legendre_sphPlm_impl(lmax, m, x, &y, result_array);
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

double gsl_sf_legendre_Plm(const int l, const int m, const double x)
{
  double y;
  int status = gsl_sf_legendre_Plm_impl(l, m, x, &y, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Plm");
  }
  return y;
}

double gsl_sf_legendre_sphPlm(const int l, const int m, const double x)
{
  double y;
  int status = gsl_sf_legendre_sphPlm_impl(l, m, x, &y, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_sphPlm");
  }
  return y;
}
