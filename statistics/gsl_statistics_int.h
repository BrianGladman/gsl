/* This is an automatically generated file created by
   convert_double_to_int.pl -- do not edit                */

#ifndef _GSL_STATISTICS_INT_H
#define _GSL_STATISTICS_INT_H

#include <stddef.h>

double gsl_stats_int_mean (const int data[], const size_t n);
double gsl_stats_int_variance (const int data[], const size_t n);
double gsl_stats_int_sd (const int data[], const size_t n);
double gsl_stats_int_est_variance (const int data[], const size_t n);
double gsl_stats_int_est_sd (const int data[], const size_t n);
double gsl_stats_int_absdev (const int data[], const size_t n);

double gsl_stats_int_skew (const int data[], const size_t n);
double gsl_stats_int_kurtosis (const int data[], const size_t n);

double 
gsl_stats_int_variance_with_mean (const int data[], 
			      const size_t n, 
			      const double mean);
double 
gsl_stats_int_sd_with_mean (const int data[], 
			const size_t n, 
			const double mean);
double 
gsl_stats_int_est_variance_with_mean (const int data[], 
				  const size_t n, 
				  const double mean);
double 
gsl_stats_int_est_sd_with_mean (const int data[], 
			    const size_t n, 
			    const double mean);
double 
gsl_stats_int_absdev_with_mean (const int data[], 
			    const size_t n, 
			    const double mean);
double 
gsl_stats_int_skew_with_mean_and_sd (const int data[], 
				 const size_t n, 
				 const double mean, 
				 const double sd);
double 
gsl_stats_int_kurtosis_with_mean_and_sd (const int data[],
				     const size_t n,
				     const double mean,
				     const double sd);

double gsl_stats_int_pvariance (const int data1[], const int data2[],
			    const size_t n1, const size_t n2);
double gsl_stats_int_ttest (const int data1[], const int data2[],
			const size_t n1, const size_t n2);

int /* BASE type */
gsl_stats_int_max (const int data[], const size_t n);
int /* BASE type */
gsl_stats_int_min (const int data[], const size_t n);

size_t gsl_stats_int_max_index (const int data[], const size_t n);
size_t gsl_stats_int_min_index (const int data[], const size_t n);

int gsl_stats_int_sort_data (int data[], const size_t n) ;

double gsl_stats_int_median_from_sorted_data (const int sorted_data[],
					  const size_t n) ;
double gsl_stats_int_percentile_from_sorted_data (const int sorted_data[],
					      const size_t n, const double f) ;

#endif /* _GSL_STATISTICS_INT_H */
