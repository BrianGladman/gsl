/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SUM_H_
#define GSL_SUM_H_


/* Basic Levin-u acceleration method.
 * No derivative information.
 * Based on WHIZ(), TOMS-602.
 *
 *   array       = array of series elements
 *   array_size  = size of array
 *   q_num       = backward diagonal of numerator; length >= array_size
 *   q_den       = backward diagonal of numerator; length >= array_size
 *   sum_accel   = result of summation acceleration
 *   sum_plain   = simple sum of series
 *   precision   = precision estimate
 *
 * See [Fessler et al., ACM TOMS 9, 346 (1983) and TOMS-602]
 */
int gsl_sum_levin_u_impl(const double * array, unsigned int array_size,
                         double * q_num,
                         double * q_den,
                         double * sum_accel,
                         double * sum_plain,
                         double * precision);
int gsl_sum_levin_u_e(const double * array, unsigned int array_size,
                      double * q_num,
                      double * q_den,
                      double * sum_accel,
                      double * sum_plain,
                      double * precision);


/* Basic Levin-u acceleration method
 * with constraints on the terms used.
 * No derivative information.
 * Based on WHIZ(), TOMS-602.
 *
 *   array       = array of series elements
 *   array_size  = size of array
 *   q_num       = backward diagonal of numerator; length >= array_size
 *   min_terms   = minimum number of terms to sum
 *   max_terms   = maximum number of terms to sum
 *   q_den       = backward diagonal of numerator; length >= array_size
 *   sum_accel   = result of summation acceleration
 *   sum_plain   = simple sum of series
 *   precision   = precision estimate
 *
 * See [Fessler et al., ACM TOMS 9, 346 (1983) and TOMS-602]
 */
int gsl_sum_levin_u_con_impl(const double * array, unsigned int array_size,
                             unsigned int min_terms, unsigned int max_terms,
                             double * q_num,
                             double * q_den,
                             double * sum_accel,
                             double * sum_plain,
                             double * precision);
int gsl_sum_levin_u_con_e(const double * array, unsigned int array_size,
                          int min_terms, int max_terms,
                          double * q_num,
                          double * q_den,
                          double * sum_accel,
                          double * sum_plain,
                          double * precision);


/* Basic Levin-u step w/o reference to the array of terms.
 * We only need to specify the value of the current term
 * to execute the step. Based on version of "WHIZ()" modified
 * by Goano, TOMS-745 [really closer to original WHIZ1()].
 *
 * sum = t0 + ... + t_{n-1} + term;  term = t_{n}
 *
 *   term   = value of the series term to be added
 *   n      = position of term in series (starting from 0)
 *   q_num  = backward diagonal of numerator, assumed to have size >= n + 1
 *   q_den  = backward diagonal of numerator, assumed to have size >= n + 1
 *   sum_accel = result of summation acceleration
 *   sum_plain = simple sum of series
 */
int
gsl_sum_levin_u_step_impl(double term,
                          unsigned int n,
                          double * q_num,
                          double * q_den,
                          double * sum_accel,
                          double * sum_plain);
int
gsl_sum_levin_u_step_e(double term,
                       unsigned int n,
                       double * q_num,
                       double * q_den,
                       double * sum_accel,
                       double * sum_plain);


#endif  /* !GSL_SUM_H_ */
