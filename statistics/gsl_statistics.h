#ifndef _GSL_STATISTICS_H
#define _GSL_STATISTICS_H

#include <stddef.h>

double gsl_stats_mean (const double data[], size_t n);
double gsl_stats_variance (const double data[], size_t n);
double gsl_stats_sd (const double data[], size_t n);
double gsl_stats_est_variance (const double data[], size_t n);
double gsl_stats_est_sd (const double data[], size_t n);
double gsl_stats_absdev (const double data[], size_t n);

double gsl_stats_skew (const double data[], size_t n);
double gsl_stats_kurtosis (const double data[], size_t n);

double gsl_stats_variance_with_mean (const double data[], size_t n, 
				     double mean);
double gsl_stats_sd_with_mean (const double data[], size_t n, double mean);
double gsl_stats_est_variance_with_mean (const double data[], size_t n,
				  double mean);
double gsl_stats_est_sd_with_mean (const double data[], size_t n, double mean);
double gsl_stats_absdev_with_mean (const double data[], size_t n, double mean);
double gsl_stats_skew_with_mean_and_sd (const double data[], size_t n, 
					double mean, double sd);
double gsl_stats_kurtosis_with_mean_and_sd (const double data[], size_t n,
					    double mean,  double sd);

double gsl_stats_pvariance (const double data1[], const double data2[],
			    size_t n1, size_t n2);
double gsl_stats_ttest (const double data1[], const double data2[],
			size_t n1, size_t n2);

double gsl_stats_max (const double data[], size_t n);
double gsl_stats_min (const double data[], size_t n);

size_t gsl_stats_max_index (const double data[], size_t n);
size_t gsl_stats_min_index (const double data[], size_t n);

int gsl_stats_sort_data (double data[], size_t n) ;

double gsl_stats_median_from_sorted_data (const double sorted_data[],
					  size_t n) ;
double gsl_stats_percentile_from_sorted_data (const double sorted_data[],
					      size_t n, const double f) ;

#endif /* _GSL_STATISTICS_H */
