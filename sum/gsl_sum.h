/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SUM_H
#define GSL_SUM_H


/* Basic Levin-u acceleration method.
 * No derivative information.
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

int gsl_sum_levin_u_accel (const double * array, size_t array_size,
			   double * q_num, double * q_den,
			   double * sum_accel, double * sum_plain,
			   double * precision);

/* Basic Levin-u acceleration method
 * with constraints on the terms used.
 * No derivative information.
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

int gsl_sum_levin_u_accel_minmax (const double * array, size_t array_size,
				  size_t min_terms, size_t max_terms,
				  double * q_num,
				  double * q_den,
				  double * sum_accel,
				  double * sum_plain,
				  double * precision)
;
/* Basic Levin-u step w/o reference to the array of terms.
 * We only need to specify the value of the current term
 * to execute the step. See TOMS-745.
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

int gsl_sum_levin_u_step(double term,
			 size_t n,
			 double * q_num,
			 double * q_den,
			 double * sum_accel,
			 double * sum_plain);

int
gsl_sum_levin_u_with_derivs_step (const double term,
				  const size_t n,
				  const size_t nmax,
				  double *q_num,
				  double *q_den,
				  double *dq_num,
				  double *dq_den,
				  double *dsum,
				  double *sum_accel,
				  double *sum_plain);

int gsl_sum_levin_u_with_derivs_accel (const double * array, 
				       size_t array_size,
				       double * q_num,
				       double * q_den,
				       double * dq_num,
				       double * dq_den,
				       double * dsum,
				       double * sum_accel,
				       double * sum_plain,
				       double * precision);

int gsl_sum_levin_u_with_derivs_accel_minmax (const double * array, 
					      size_t array_size,
					      size_t min_terms, 
					      size_t max_terms,
					      double * q_num,
					      double * q_den,
					      double * dq_num,
					      double * dq_den,
					      double * dsum,
					      double * sum_accel,
					      double * sum_plain,
					      double * precision);


#endif  /* GSL_SUM_H */



