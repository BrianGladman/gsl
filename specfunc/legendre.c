#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"

double gsl_sf_legendre(int l, int m, double x)
{
  int i;
  double dif, sum, exp_check;
  
  /* Starting value for recursion, P_m^m(x) */
  double pmm = 1.0;

  /* Strip sign of m. */
  int am = m > 0 ? m : -m;

  /* Check arguments. */
  if(m < 0){
    char buff[100];
    sprintf(buff, "gsl_sf_legendre: negative m= %d  will be sign-stripped", m);
    GSL_ERROR_MESSAGE(buff, GSL_EDOM);
  }
  if(am > l){
    char buff[100];
    sprintf(buff, "gsl_sf_legendre: bad m argument for legendre: m= %d  l= %d",m,l);
    GSL_ERROR_MESSAGE(buff, GSL_EDOM);
    return 0.;
  }
  if(fabs(x) > 1.){
    char buff[100];
    sprintf(buff, "gsl_sf_legendre: bad x argument for legendre: x= %g", x);
    GSL_ERROR_MESSAGE(buff, GSL_EDOM);
    return 0.;
  }
  
  /* If l is large and m is large, then we have to worry
   * about overflow. Calculate an approximate exponent which
   * measures the normalization of this thing.
   */
  dif = l-am;
  sum = l+am;
  exp_check = 0.5 * log(2.*l+1.) 
    + 0.5 * dif * (log(dif)-1.)
      - 0.5 * sum * (log(sum)-1.);

  /* Bail out if if looks bad.
   */  
  if(exp_check < GSL_LOG_DBL_MIN + 10.){
    char buff[100];
    sprintf(buff, "gsl_sf_legendre: apparent overflow condition: l= %d  m= %d",l,m);
    GSL_ERROR_MESSAGE(buff, GSL_EDOM);
    return 0.;
  }

  /* If m > 0, then calculate P_m^m from the analytic result:
   * 
   *          P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2)
   */
  if(m > 0){
    double somx2 = sqrt((1.-x)*(1.+x));
    double fact = 1.;

    for(i=1; i<=m; i++){
      pmm *= -fact * somx2;
      fact += 2.;
    }
  }

  if(l == m){
    /* If we are already there, then we are done. */
    return pmm;
  }
  else{
    
    /* Otherwise calculate P_{m+1}^m. */
    double pmmp1 = x * (2*m + 1) * pmm;

    if(l == (m+1)){
      /* If we are just at the penultimate step, then we are done. */
      return pmmp1;
    }
    else {
      /* Otherwise compute P_l^m, l > m+1, recursively. */
      double pll;
      int ll;

      for(ll=m+2; ll <= l; ll++){
	pll = (x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm) / (ll-m);
	pmm = pmmp1;
	pmmp1 = pll;
      }
      
      return pll;
    }
  }
}


double gsl_sf_Ylm_legendre(int l, int m, double x)
{
  int i;
  double sgn_m, sgn_factor;

  /* Starting value for recursion, Y_m^m(x).
   * We include part of the normalization factor here.
   */
  double ymm = sqrt((2.*(double)l+1.) / (4.*M_PI));

  /* Check arguments. */
  if(abs(m) > l){
    char buff[100];
    sprintf(buff, "Ylm_legendre: bad m argument: m= %d  l= %d", m, l);
    GSL_ERROR_MESSAGE(buff, GSL_EDOM);
    return 0.;
  }
  if(fabs(x) > 1.){
    char buff[100];
    sprintf(buff, "Ylm_legendre: bad x argument: x= %g", x);
    GSL_ERROR_MESSAGE(buff, GSL_EDOM);
    return 0.;
  }

  /* strip sign of m */
  sgn_m = (m < 0 ? -1. : 1.);
  m *= sgn_m;
  sgn_factor = (sgn_m < 0. ? pow_int(-1.,m) : 1.);

  /* If m > 0, then calculate Y_m^m from the analytic result.
   * If m==0, then we don't have to do anything; ymm is ready to go.
   * 
   *          Y_m^m(x) = sqrt(1/(2m)!) (-1)^m (2m-1)!! (1-x^2)^(m/2)
   */
  if(m > 0){
    double somx2 = sqrt((1.-x)*(1.+x));
    double fact1 = 1.;
    double fact2 = 1. / sqrt(2.);

    for(i=1; i<=m; i++){
      ymm *= -fact1 * fact2 * somx2;
      fact1 += 2.;
      fact2 = 1. / sqrt(fact1 * (fact1 + 1.));
    }
  }

  if(l == m){
    /* If we are already there, then we are done. */
    return sgn_factor * ymm;
  }
  else{
    /* Otherwise calculate Y_{m+1}^m. */
    double ymmp1 = x * sqrt(2.*(double)m + 1) * ymm;

    if(l == (m+1)){
      /* If we are just at the penultimate step, then we are done. */
      return sgn_factor * ymmp1;
    }
    else {
      /* Otherwise compute Y_l^m, l > m+1, recursively. */
      double yll;
      int ll;
      
      for(ll=m+2; ll <= l; ll++){
	double factor1 = sqrt((double)(ll-m) / (double)(ll+m));
	double factor2 = factor1 * sqrt((double)(ll-m-1.) / (double)(ll+m-1.));
	yll = (x*(2.*ll-1)*ymmp1*factor1 - (ll+m-1.)*ymm*factor2) / (ll-m);
	ymm = ymmp1;
	ymmp1 = yll;
      }
      
      return sgn_factor * yll;
    }
  }
}
