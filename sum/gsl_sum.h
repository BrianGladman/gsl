/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SUM_H__
#define __GSL_SUM_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/* Basic Levin-u acceleration method.
 *
 *   array       = array of series elements
 *   array_size  = size of array
 *   q_num       = backward diagonal of numerator; length >= array_size
 *   q_den       = backward diagonal of denominator; length >= array_size
 *   dq_num      = table of numerator derivatives; length >= array_size**2
 *   dq_den      = table of denominator derivatives; length >= array_size**2
 *   dsum        = derivative of sum wrt term i; length >= array_size
 *   sum_accel   = result of summation acceleration
 *   sum_plain   = simple sum of series
 *   precision   = precision estimate
 *
 * See [Fessler et al., ACM TOMS 9, 346 (1983) and TOMS-602]
 */

int gsl_sum_levin_u_accel (const double * array, 
			   size_t array_size,
			   double * q_num,
			   double * q_den,
			   double * dq_num,
			   double * dq_den,
			   double * dsum,
			   double * sum_accel,
			   size_t * n_used,
			   double * sum_plain,
			   double * precision);

int gsl_sum_levin_u_trunc_accel (const double * array, size_t array_size,
				 double * q_num, double * q_den,
				 double * sum_accel, 
				 size_t * n_used,
				 double * sum_plain,
				 double * precision);

/* Basic Levin-u acceleration method with constraints on the terms
 * used.
 *
 *   array       = array of series elements
 *   array_size  = size of array
 *   q_num       = backward diagonal of numerator; length >= array_size
 *   min_terms   = minimum number of terms to sum
 *   max_terms   = maximum number of terms to sum
 *   q_den       = backward diagonal of numerator; length >= array_size
 *   dq_num      = table of numerator derivatives; length >= array_size**2
 *   dq_den      = table of denominator derivatives; length >= array_size**2
 *   dsum        = derivative of sum wrt term i; length >= array_size
 *   sum_accel   = result of summation acceleration
 *   sum_plain   = simple sum of series
 *   precision   = precision estimate
 *
 * See [Fessler et al., ACM TOMS 9, 346 (1983) and TOMS-602] */


int gsl_sum_levin_u_accel_minmax (const double * array, 
				  size_t array_size,
				  size_t min_terms, 
				  size_t max_terms,
				  double * q_num,
				  double * q_den,
				  double * dq_num,
				  double * dq_den,
				  double * dsum,
				  double * sum_accel,
				  size_t * n_used,
				  double * sum_plain,
				  double * precision);

int gsl_sum_levin_u_trunc_accel_minmax (const double * array, 
					size_t array_size,
					size_t min_terms, 
					size_t max_terms,
					double * q_num,
					double * q_den,
					double * sum_accel,
					size_t * n_used,
					double * sum_plain,
					double * precision);

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
 *   dq_num = table of numerator derivatives; length >= array_size**2
 *   dq_den = table of denominator derivatives; length >= array_size**2
 *   dsum   = derivative of sum wrt term i; length >= array_size
 *   sum_accel = result of summation acceleration
 *   sum_plain = simple sum of series
 */

int
gsl_sum_levin_u_step (double term,
		      size_t n,
		      size_t nmax,
		      double *q_num,
		      double *q_den,
		      double *dq_num,
		      double *dq_den,
		      double *dsum,
		      double *sum_accel,
		      double *sum_plain);

int gsl_sum_levin_u_trunc_step(double term,
			       size_t n,
			       double * q_num,
			       double * q_den,
			       double * sum_accel,
			       double * sum_plain);

__END_DECLS

#endif /* __GSL_SUM_H__ */


