/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_log.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_legendre_P1_impl(double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
}

int
gsl_sf_legendre_P2_impl(double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else {
    result->val = 0.5*(3.0*x*x - 1.0);
    result->err = GSL_DBL_EPSILON * (fabs(3.0*x*x) + 1.0);
    return GSL_SUCCESS;
  }
}

int
gsl_sf_legendre_P3_impl(double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else {
    result->val = 0.5*x*(5.0*x*x - 3.0);
    result->err = GSL_DBL_EPSILON * (fabs(result->val) + 0.5 * fabs(x) * (fabs(5.0*x*x) + 3.0));
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Pl_impl(const int l, const double x, gsl_sf_result * result)
{ 
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(l < 0 || x < -1.0 || x > 1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(l == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l == 1) {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l == 2) {
    result->val = 0.5 * (3.0*x*x - 1.0);
    result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x == -1.0) {
    result->val = ( GSL_IS_ODD(l) ? -1.0 : 1.0 );
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l < 100000) {
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

    result->val = p_ell;
    result->err = (0.5 * ell + 1.0) * GSL_DBL_EPSILON * fabs(p_ell);
    return GSL_SUCCESS;
  }
  else {
    /* Asymptotic expansion.
     * [Olver, p. 473]
     */
    double u  = l + 0.5;
    double th = acos(x);
    gsl_sf_result J0;
    gsl_sf_result Jm1;
    int stat_J0  = gsl_sf_bessel_J0_impl(u*th, &J0);
    int stat_Jm1 = gsl_sf_bessel_Jn_impl(-1, u*th, &Jm1);
    double pre;
    double B00;
    double c1;

    /* B00 = 1/8 (1 - th cot(th) / th^2
     * pre = sqrt(th/sin(th))
     */
    if(th < GSL_ROOT4_MACH_EPS) {
      B00 = (1.0 + th*th/15.0)/24.0;
      pre = 1.0 + th*th/12.0;
    }
    else {
      double sin_th = sqrt(1.0 - x*x);
      double cot_th = x / sin_th;
      B00 = 1.0/8.0 * (1.0 - th * cot_th) / (th*th);
      pre = sqrt(th/sin_th);
    }

    c1 = th/u * B00;

    result->val  = pre * (J0.val + c1 * Jm1.val);
    result->err  = pre * (J0.err + fabs(c1) * Jm1.err);
    result->err += GSL_SQRT_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_J0, stat_Jm1);
  }
}


int
gsl_sf_legendre_Pl_array_impl(const int lmax, const double x, double * result_array)
{
  if(result_array == 0) {
    return GSL_EFAULT;
  }
  else if(lmax < 0 || x < -1.0 || x > 1.0) {
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
gsl_sf_legendre_Plm_impl(const int l, const int m, const double x, gsl_sf_result * result)
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

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(m < 0 || l < m || x < -1.0 || x > 1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(exp_check < GSL_LOG_DBL_MIN + 10.0){
    /* Bail out.
     */
    result->val = 0.0;
    result->err = 0.0;
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
      result->val = pmm;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(pmm);
      return GSL_SUCCESS;
    }
    else if(l == m + 1) {
      result->val = pmmp1;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(pmmp1);
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

      result->val = p_ell;
      result->err = (0.5*(l-m) + 1.0) * GSL_DBL_EPSILON * fabs(p_ell);

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

  if(result_array == 0) {
    return GSL_EFAULT;
  }
  else if(m < 0 || lmax < m || x < -1.0 || x > 1.0) {
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
gsl_sf_legendre_sphPlm_impl(const int l, int m, const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(m < 0 || l < m || x < -1.0 || x > 1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(m == 0) {
    gsl_sf_result P;
    int stat_P = gsl_sf_legendre_Pl_impl(l, x, &P);
    double pre = sqrt((2.0*l + 1.0)/(4.0*M_PI));
    result->val  = pre * P.val;
    result->err  = pre * P.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_P;
  }
  else if(x == 1.0 || x == -1.0) {
    /* m > 0 here */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* m > 0 and |x| < 1 here */

    /* Starting value for recursion.
     * Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) ) (-1)^m (1-x^2)^(m/2) / pi^(1/4)
     */
    gsl_sf_result lncirc;
    gsl_sf_result lnpoch;
    double lnpre;
    double sgn = ( GSL_IS_ODD(m) ? -1.0 : 1.0);
    double ymm, ymmerr;
    double ymmp1_factor = x * sqrt(2.0*m + 3.0);
    double ymmp1;
    gsl_sf_log_1plusx_impl(-x*x, &lncirc);
    gsl_sf_lnpoch_impl(m, 0.5, &lnpoch);  /* Gamma(m+1/2)/Gamma(m) */
    lnpre = -0.25*M_LNPI + 0.5 * (lnpoch.val + m*lncirc.val);
    ymm   = sqrt((2.0+1.0/m)/(4.0*M_PI)) * sgn * exp(lnpre);
    ymmp1 = ymmp1_factor * ymm;
    ymmerr = (lnpoch.err + fabs(m*lncirc.err)) * GSL_DBL_EPSILON * fabs(ymm);

    if(l == m){
      result->val  = ymm;
      result->err  = ymmerr;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(ymm);
      return GSL_SUCCESS;
    }
    else if(l == m + 1) {
      result->val  = ymmp1;
      result->err  = fabs(ymmp1_factor) * ymmerr;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(ymmp1);
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

      result->val  = y_ell;
      result->err  = (0.5*(l-m) + 1.0) * GSL_DBL_EPSILON * fabs(y_ell);
      result->err += fabs(ymmerr/ymm) * fabs(y_ell);

      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_legendre_sphPlm_array_impl(const int lmax, int m, const double x, double * result_array)
{
  if(result_array == 0) {
    return GSL_EFAULT;
  }
  else if(m < 0 || lmax < m || x < -1.0 || x > 1.0) {
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

      gsl_sf_result lncirc;
      gsl_sf_result lnpoch;
      double lnpre;
      double sgn = ( GSL_IS_ODD(m) ? -1.0 : 1.0);
      gsl_sf_log_1plusx_impl(-x*x, &lncirc);
      gsl_sf_lnpoch_impl(m, 0.5, &lnpoch);  /* Gamma(m+1/2)/Gamma(m) */
      lnpre = -0.25*M_LNPI + 0.5 * (lnpoch.val + m*lncirc.val);
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

int gsl_sf_legendre_P1_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_legendre_P1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_P1_e", status);
  }
  return status;
}

int gsl_sf_legendre_P2_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_legendre_P2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_P2_e", status);
  }
  return status;
}

int gsl_sf_legendre_P3_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_legendre_P3_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_P3_e", status);
  }
  return status;
}

int gsl_sf_legendre_Pl_e(const int l, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_legendre_Pl_impl(l, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Pl_e", status);
  }
  return status;
}

int gsl_sf_legendre_Plm_e(const int l, const int m, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_legendre_Plm_impl(l, m, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Plm_e", status);
  }
  return status;
}

int gsl_sf_legendre_sphPlm_e(const int l, const int m, const double x, gsl_sf_result * result)
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
