/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_test.h>
#include <gsl_errno.h>
#include "gsl_sum.h"

#define locMAX(a,b)  ((a) > (b) ? (a) : (b))


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/


int
gsl_sum_levin_u_step_impl(const double term,
                          const unsigned int n,
                          double * q_num,
                          double * q_den,
                          double * sum_accel,
                          double * sum_plain)
{
  if(term == 0.0) {
    /* This is actually harmless when treated in
     * this way. A term which is exactly zero is
     * simply ignored; the state is not changed.
     * We return GSL_EZERODIV as an indicator that
     * this occured.
     */
    return GSL_EZERODIV;
  }
  else if(n == 0) {
    *sum_accel = term;
    *sum_plain = term;
    q_den[0] = 1.0/term;
    q_num[0] = 1.0;
    return GSL_SUCCESS;
  }
  else {
    double factor = 1.0;
    double ratio  = (double)n/(n+1.0);
    int j;

    *sum_plain += term;
    q_den[n] = 1.0/(term*(n+1.0)*(n+1.0));
    q_num[n] = *sum_plain * q_den[n];

    for(j=n-1; j>=0; j--) {
      double c = factor*(j+1)/(n+1);
      factor *= ratio;
      q_den[j] = q_den[j+1] - c*q_den[j];
      q_num[j] = q_num[j+1] - c*q_num[j];
    }

    *sum_accel = q_num[0] / q_den[0];
    return GSL_SUCCESS;
  }
}


int
gsl_sum_levin_u_con_impl(const double * array, const unsigned int array_size,
                         const unsigned int min_terms,
                         const unsigned int max_terms,
                         double * q_num,
                         double * q_den,
                         double * sum_accel,
                         double * sum_plain,
		         double * precision)
{
  if(array_size == 0) {
    *sum_accel = 0.0;
    *sum_plain = 0.0;
    return GSL_SUCCESS;
  }
  else if(array_size == 1) {
    *sum_accel = array[0];
    *sum_plain = array[0];
    return GSL_SUCCESS;
  }
  else {
    const double SMALL = 0.01;
    const int nmax = locMAX(max_terms, array_size) - 1;
    double noise_n  = 0.0,  noise_nm1 = 0.0;
    double trunc_n  = 0.0,  trunc_nm1 = 0.0;
    double result_n = 0.0, result_nm1 = 0.0;
    int n;
    int better = 0;
    int before = 0;
    int converging = 0;
    double least_trunc = DBL_MAX;
    double result_least_trunc;

    /* Calculate specified minimum number of terms.
     * No convergence tests are made, and no
     * truncation information is stored.
     */
    for(n=0; n<min_terms; n++) {
      const double t = array[n];
      int status;

      noise_nm1  = noise_n;
      noise_n    = 3.0*GSL_MACH_EPS*fabs(t);

      result_nm1 = result_n;
      status = gsl_sum_levin_u_step_impl(t, n, q_num, q_den, &result_n, sum_plain);

      if(status != GSL_SUCCESS && status != GSL_EZERODIV) {
        *sum_accel = result_nm1;
        return status;
      }
    }

    /* Assume the result after the minimum
     * calculation is the best.
     */
    result_least_trunc = result_n;

    /* Calculate up to maximum number of terms.
     * Check truncation condition.
     */
    for(; n<=nmax; n++) {
      const double t = array[n];
      int status;

      noise_nm1  = noise_n;
      noise_n    = 3.0*GSL_MACH_EPS*fabs(t);

      result_nm1 = result_n;
      status = gsl_sum_levin_u_step_impl(t, n, q_num, q_den, &result_n, sum_plain);

      if(status != GSL_SUCCESS && status != GSL_EZERODIV) {
        *sum_accel = result_nm1;
        return status;
      }

      trunc_nm1  = trunc_n;
      trunc_n    = fabs(result_n - result_nm1);

      /* Determine if we are in the
       * convergence region.
       */
      better     = ( trunc_n < trunc_nm1 || trunc_n < SMALL*fabs(result_n) );
      converging = converging || (better && before);
      before     = better;

      if(converging) {
        if(trunc_n < least_trunc) {
	  /* Found a low truncation point in
	   * the convergence region. Save it.
	   */
          least_trunc = trunc_n;
          result_least_trunc = result_n;
        }

	if(fabs(trunc_n/result_n) < 10.0*GSL_MACH_EPS) break;
      }
    }

    if(converging) {
      /* Stopped in the convergence region.
       * Return result and error estimate.
       */
      *sum_accel = result_least_trunc;
      *precision = fabs(least_trunc / *sum_accel);
      return GSL_SUCCESS;
    }
    else {
      /* Never reached the convergence region.
       * Use the last calculated values.
       */
      *sum_accel = result_n;
      *precision = fabs(result_n/trunc_n);
      return GSL_SUCCESS;
    }
  }
}


int
gsl_sum_levin_u_impl(const double * array, const unsigned int array_size,
                     double * q_num,
                     double * q_den,
                     double * sum_accel,
                     double * sum_plain,
		     double * precision)
{
  return gsl_sum_levin_u_con_impl(array, array_size,
                                  0, array_size-1,
                                  q_num, q_den,
                                  sum_accel, sum_plain, precision);
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sum_levin_u_step_e(const double t,
                       const unsigned int n,
                       double * q_num,
                       double * q_den,
                       double * sum_accel,
                       double * sum_plain)
{
  int status = gsl_sum_levin_u_step_impl(t, n, q_num, q_den, sum_accel, sum_plain);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sum_levin_u_step_e", status);
  }
  return status;
}


int
gsl_sum_levin_u_con_e(const double * array,
                      const unsigned int array_size,
		      const int min_terms, const int max_terms,
                      double * q_num,
                      double * q_den,
                      double * sum_accel,
                      double * sum_plain,
                      double * precision)
{
  int status = gsl_sum_levin_u_con_impl(array, array_size,
                                        min_terms, max_terms,
                                        q_num, q_den,
                                        sum_accel, sum_plain, precision);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sum_levin_u_con_e", status);
  }
  return status;
}


int
gsl_sum_levin_u_e(const double * array,
                  const unsigned int array_size,
                  double * q_num,
                  double * q_den,
                  double * sum_accel,
                  double * sum_plain,
                  double * precision)
{
  int status = gsl_sum_levin_u_impl(array, array_size,
                                    q_num, q_den,
                                    sum_accel, sum_plain, precision);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sum_levin_u_e", status);
  }
  return status;
}
