/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_log.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/


int
gsl_sf_legendre_Pl_impl(const int l, const double x, double * result)
{ 
  if(l < 0 || x < -1.0 || x > 1.0) {
    return GSL_EDOM;
  }
  else if(l == 0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(l == 1) {
    *result = x;
    return GSL_SUCCESS;
  }
  else if(l == 2) {
    *result = 0.5 * (3.0*x*x - 1.0);
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(x == -1.0) {
    *result = ( GSL_IS_ODD(l) ? -1.0 : 1.0 );
    return GSL_SUCCESS;
  }
  else if(l < 1.0e+06) {
    /* Compute by upward recurrence on l.
     */
    double pmm   = 1.0;     /* P_0(x) */
    double pmmp1 = x;	    /* P_1(x) */
    double p_ell = pmmp1;
    int ell;

    for(ell=2; ell <= l; ell++){
      p_ell = (x*(2*ell-1)*pmmp1 - (ell-1)*pmm) / ell;
      pmm = pmmp1;
      pmmp1 = p_ell;
    }

    *result = p_ell;
    return GSL_SUCCESS;
  }
  else {
    /* Asymptotic expansion.
     * FIXME: need another term or two here
     */
    double th  = acos(x);
    double pre = sqrt(th/sin(th));
    double J0;
    int stat_J0 = gsl_sf_bessel_J0_impl((l+0.5)*th, &J0);
    *result = pre * J0;
    return stat_J0;
  }
}


int
gsl_sf_legendre_Pl_array_impl(const int lmax, const double x, double * result_array)
{
  if(lmax < 0 || x < -1.0 || x > 1.0) {
    return GSL_EDOM;
  }
  else if(lmax == 0) {
    result_array[0] = 1.0;
    return GSL_SUCCESS;
  }
  else if(lmax == 1) {
    result_array[0] = 1.0;
    result_array[1] = x;
    return GSL_SUCCESS;
  }
  else {
    double pmm   = 1.0;     /* P_0(x) */
    double pmmp1 = x;	    /* P_1(x) */
    double p_ell = pmmp1;
    int ell;

    result_array[0] = 1.0;
    result_array[1] = x;

    for(ell=2; ell <= lmax; ell++){
      p_ell = (x*(2*ell-1)*pmmp1 - (ell-1)*pmm) / ell;
      pmm = pmmp1;
      pmmp1 = p_ell;
      result_array[ell] = p_ell;
    }

    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Plm_impl(const int l, const int m, const double x, double * result)
{
  /* If l is large and m is large, then we have to worry
   * about overflow. Calculate an approximate exponent which
   * measures the normalization of this thing.
   */
  double dif = l-m;
  double sum = l+m;
  double exp_check = 0.5 * log(2.0*l+1.0) 
                     + 0.5 * dif * (log(dif)-1.0)
                     - 0.5 * sum * (log(sum)-1.0);

  if(m < 0 || l < m || x < -1.0 || x > 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(exp_check < GSL_LOG_DBL_MIN + 10.0){
    /* Bail out.
     */
    *result = 0.0;
    return GSL_EOVRFLW;
  }
  else {
    double pmm;                 /* P_m^m(x) */
    double pmmp1;               /* P_{m+1}^m(x) */

    /* Calculate P_m^m from the analytic result:
     *          P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2) , m > 0
     */
    pmm = 1.0;
    if(m > 0){
      double circ = sqrt(1.0-x)*sqrt(1.0+x);
      double fact = 1.0;
      int i;
      for(i=1; i<=m; i++) {
        pmm  *= -fact * circ;
        fact += 2.0;
      }
    }

    /* Calculate P_{m+1}^m. */
    pmmp1 = x * (2*m + 1) * pmm;

    if(l == m){
      *result = pmm;
      return GSL_SUCCESS;
    }
    else if(l == m + 1) {
      *result = pmmp1;
      return GSL_SUCCESS;
    }
    else{
      double p_ell;
      int ell;
    
      /* Compute P_l^m, l > m+1 by upward recurrence on l. */
      for(ell=m+2; ell <= l; ell++){
        p_ell = (x*(2*ell-1)*pmmp1 - (ell+m-1)*pmm) / (ell-m);
        pmm = pmmp1;
        pmmp1 = p_ell;
      }

      *result = p_ell;
      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_legendre_Plm_array_impl(const int lmax, const int m, const double x, double * result_array)
{
  /* If l is large and m is large, then we have to worry
   * about overflow. Calculate an approximate exponent which
   * measures the normalization of this thing.
   */
  double dif = lmax-m;
  double sum = lmax+m;
  double exp_check = 0.5 * log(2.0*lmax+1.0) 
                     + 0.5 * dif * (log(dif)-1.0)
                     - 0.5 * sum * (log(sum)-1.0);

  if(m < 0 || lmax < m || x < -1.0 || x > 1.0) {
    return GSL_EDOM;
  }
  else if(m > 0 && (x == 1.0 || x == -1.0)) {
    int ell;
    for(ell=m; ell<=lmax; ell++) result_array[ell-m] = 0.0;
    return GSL_SUCCESS;
  }
  else if(exp_check < GSL_LOG_DBL_MIN + 10.0){
    /* Bail out.
     */
    return GSL_EOVRFLW;
  }
  else {
    double pmm;                 /* P_m^m(x)     */
    double pmmp1;               /* P_{m+1}^m(x) */

    /* Calculate P_m^m from the analytic result:
     *          P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2) , m > 0
     */
    pmm = 1.0;
    if(m > 0){
      int i;
      double circ = sqrt(1.0-x)*sqrt(1.0+x);
      double fact = 1.0;
      for(i=1; i<=m; i++){
        pmm  *= -fact * circ;
        fact += 2.0;
      }
    }

    /* Calculate P_{m+1}^m. */
    pmmp1 = x * (2*m + 1) * pmm;

    if(lmax == m){
      result_array[0] = pmm;
      return GSL_SUCCESS;
    }
    else if(lmax == m + 1) {
      result_array[0] = pmm;
      result_array[1] = pmmp1;
      return GSL_SUCCESS;
    }
    else{
      double p_ell;
      int ell;

      result_array[0] = pmm;
      result_array[1] = pmmp1;

      /* Compute P_l^m, l >= m+2, by upward recursion on l. */
      for(ell=m+2; ell <= lmax; ell++){
        p_ell = (x*(2*ell-1)*pmmp1 - (ell+m-1)*pmm) / (ell-m);
        pmm = pmmp1;
        pmmp1 = p_ell;
        result_array[ell-m] = p_ell;
      }

      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_legendre_sphPlm_impl(const int l, int m, const double x, double * result)
{
  if(m < 0 || l < m || x < -1.0 || x > 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(m == 0) {
    double P;
    int stat_P = gsl_sf_legendre_Pl_impl(l, x, &P);
    *result = sqrt((2.0*l + 1.0)/(4.0*M_PI)) * P;
    return stat_P;
  }
  else if(x == 1.0 || x == -1.0) {
    /* m > 0 here */
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* m > 0 and |x| < 1 here */

    /* Starting value for recursion.
     * Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) ) (-1)^m (1-x^2)^(m/2) / pi^(1/4)
     */
    double lncirc;
    double lnpoch;
    double lnpre;
    double sgn = ( GSL_IS_ODD(m) ? -1.0 : 1.0);
    double ymm;
    double ymmp1;
    gsl_sf_log_1plusx_impl(-x*x, &lncirc);
    gsl_sf_lnpoch_impl(m, 0.5, &lnpoch);  /* Gamma(m+1/2)/Gamma(m) */
    lnpre = -0.25*M_LNPI + 0.5 * (lnpoch + m*lncirc);
    ymm   = sqrt((2.0+1.0/m)/(4.0*M_PI)) * sgn * exp(lnpre);
    ymmp1 = x * sqrt(2.0*m + 3.0) * ymm;

    if(l == m){
      *result = ymm;
      return GSL_SUCCESS;
    }
    else if(l == m + 1) {
      *result = ymmp1;
      return GSL_SUCCESS;
    }
    else{
      double y_ell;
      int ell;

      /* Compute Y_l^m, l > m+1, upward recursion on l. */
      for(ell=m+2; ell <= l; ell++){
        double rat1 = (double)(ell-m)/(double)(ell+m);
	double rat2 = (ell-m-1.0)/(ell+m-1.0);
        double factor1 = sqrt(rat1*(2*ell+1)*(2*ell-1));
        double factor2 = sqrt(rat1*rat2*(2*ell+1)/(2*ell-3));
        y_ell = (x*ymmp1*factor1 - (ell+m-1)*ymm*factor2) / (ell-m);
        ymm   = ymmp1;
        ymmp1 = y_ell;
      }

      *result = y_ell;
      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_legendre_sphPlm_array_impl(const int lmax, int m, const double x, double * result_array)
{
  if(m < 0 || lmax < m || x < -1.0 || x > 1.0) {
    return GSL_EDOM;
  }
  else if(m > 0 && (x == 1.0 || x == -1.0)) {
    int ell;
    for(ell=m; ell<=lmax; ell++) result_array[ell-m] = 0.0;
    return GSL_SUCCESS;
  }
  else {
    double ymm;
    double ymmp1;

    if(m == 0) {
      ymm   = 0.5/M_SQRTPI;          /* Y00 = 1/sqrt(4pi) */
      ymmp1 = x * M_SQRT3 * ymm;
    }
    else {
      /* |x| < 1 here */

      double lncirc;
      double lnpoch;
      double lnpre;
      double sgn = ( GSL_IS_ODD(m) ? -1.0 : 1.0);
      gsl_sf_log_1plusx_impl(-x*x, &lncirc);
      gsl_sf_lnpoch_impl(m, 0.5, &lnpoch);  /* Gamma(m+1/2)/Gamma(m) */
      lnpre = -0.25*M_LNPI + 0.5 * (lnpoch + m*lncirc);
      ymm   = sqrt((2.0+1.0/m)/(4.0*M_PI)) * sgn * exp(lnpre);
      ymmp1 = x * sqrt(2.0*m + 3.0) * ymm;
    }

    if(lmax == m){
      result_array[0] = ymm;
      return GSL_SUCCESS;
    }
    else if(lmax == m + 1) {
      result_array[0] = ymm;
      result_array[1] = ymmp1;
      return GSL_SUCCESS;
    }
    else{
      double y_ell;
      int ell;

      result_array[0] = ymm;
      result_array[1] = ymmp1;

      /* Compute Y_l^m, l > m+1, upward recursion on l. */
      for(ell=m+2; ell <= lmax; ell++){
        double rat1 = (double)(ell-m)/(double)(ell+m);
	double rat2 = (ell-m-1.0)/(ell+m-1.0);
        double factor1 = sqrt(rat1*(2*ell+1)*(2*ell-1));
        double factor2 = sqrt(rat1*rat2*(2*ell+1)/(2*ell-3));
        y_ell = (x*ymmp1*factor1 - (ell+m-1)*ymm*factor2) / (ell-m);
        ymm   = ymmp1;
        ymmp1 = y_ell;
	result_array[ell-m] = y_ell;
      }
    }

    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_legendre_Pl_e(const int l, const double x, double * result)
{
  int status = gsl_sf_legendre_Pl_impl(l, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Pl_e", status);
  }
  return status;
}

int gsl_sf_legendre_Plm_e(const int l, const int m, const double x, double * result)
{
  int status = gsl_sf_legendre_Plm_impl(l, m, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Plm_e", status);
  }
  return status;
}

int gsl_sf_legendre_sphPlm_e(const int l, const int m, const double x, double * result)
{
  int status = gsl_sf_legendre_sphPlm_impl(l, m, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_sphPlm_e", status);
  }
  return status;
}

int gsl_sf_legendre_Pl_array_e(const int lmax, const double x, double * result_array)
{
  int status = gsl_sf_legendre_Pl_array_impl(lmax, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Pl_array_e", status);
  }
  return status;
}

int gsl_sf_legendre_Plm_array_e(const int lmax, const int m, const double x, double * result_array)
{
  int status = gsl_sf_legendre_Plm_array_impl(lmax, m, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Plm_array_e", status);
  }
  return status;
}

int gsl_sf_legendre_sphPlm_array_e(const int lmax, const int m, const double x, double * result_array)
{
  int status = gsl_sf_legendre_sphPlm_array_impl(lmax, m, x, result_array);
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
  int status = gsl_sf_legendre_Pl_impl(l, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Pl", status);
  }
  return y;
}

double gsl_sf_legendre_Plm(const int l, const int m, const double x)
{
  double y;
  int status = gsl_sf_legendre_Plm_impl(l, m, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Plm", status);
  }
  return y;
}

double gsl_sf_legendre_sphPlm(const int l, const int m, const double x)
{
  double y;
  int status = gsl_sf_legendre_sphPlm_impl(l, m, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_sphPlm", status);
  }
  return y;
}


double gsl_sf_legendre_P1(double x) { return x; }
double gsl_sf_legendre_P2(double x) { return 0.5*(3.0*x*x - 1.0); }
double gsl_sf_legendre_P3(double x) { return 0.5*x*(5.0*x*x - 3.0); }
double gsl_sf_legendre_P4(double x) { double x2 = x*x; return (35.0*x2*x2 -30.0*x2 + 3.0)/8.0; }
double gsl_sf_legendre_P5(double x) { double x2 = x*x; return x*(63.0*x2*x2 -70.0*x2 + 15.0)/8.0; }
