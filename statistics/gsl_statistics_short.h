#ifndef GSL_STATISTICS_SHORT_H
#define GSL_STATISTICS_SHORT_H

#include <stddef.h>

double gsl_stats_short_mean (const short data[], size_t n);
double gsl_stats_short_variance (const short data[], size_t n);
double gsl_stats_short_sd (const short data[], size_t n);
double gsl_stats_short_est_variance (const short data[], size_t n);
double gsl_stats_short_est_sd (const short data[], size_t n);
double gsl_stats_short_absdev (const short data[], size_t n);

double gsl_stats_short_skew (const short data[], size_t n);
double gsl_stats_short_kurtosis (const short data[], size_t n);
double gsl_stats_short_lag1_autocorrelation (const short data[], size_t n);

double gsl_stats_short_variance_with_mean (const short data[], size_t n, 
					 double mean);
double gsl_stats_short_sd_with_mean (const short data[], size_t n, double mean);
double gsl_stats_short_est_variance_with_mean (const short data[], size_t n,
					     double mean);
double gsl_stats_short_est_sd_with_mean (const short data[], size_t n, double mean);
double gsl_stats_short_absdev_with_mean (const short data[], size_t n, double mean);
double gsl_stats_short_skew_with_mean_and_sd (const short data[], size_t n, 
					    double mean, double sd);
double gsl_stats_short_kurtosis_with_mean_and_sd (const short data[], size_t n,
						double mean,  double sd);
double gsl_stats_short_lag1_autocorrelation_with_mean (const short data[], 
						     size_t n, double mean);

double gsl_stats_short_pvariance (const short data1[], const short data2[],
				size_t n1, size_t n2);
double gsl_stats_short_ttest (const short data1[], const short data2[],
			    size_t n1, size_t n2);

short gsl_stats_short_max (const short data[], size_t n);
short gsl_stats_short_min (const short data[], size_t n);

size_t gsl_stats_short_max_index (const short data[], size_t n);
size_t gsl_stats_short_min_index (const short data[], size_t n);

void gsl_stats_short_sort_data (short data[], size_t n) ;

double gsl_stats_short_median_from_sorted_data (const short sorted_data[],
					      size_t n) ;
double gsl_stats_short_percentile_from_sorted_data (const short sorted_data[],
						  size_t n, const double f) ;

#endif /* GSL_STATISTICS_SHORT_H */
