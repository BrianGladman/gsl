#ifndef GSL_STATISTICS_UCHAR_H
#define GSL_STATISTICS_UCHAR_H

#include <stddef.h>

double gsl_stats_uchar_mean (const unsigned char data[], size_t n);
double gsl_stats_uchar_variance (const unsigned char data[], size_t n);
double gsl_stats_uchar_sd (const unsigned char data[], size_t n);
double gsl_stats_uchar_est_variance (const unsigned char data[], size_t n);
double gsl_stats_uchar_est_sd (const unsigned char data[], size_t n);
double gsl_stats_uchar_absdev (const unsigned char data[], size_t n);

double gsl_stats_uchar_skew (const unsigned char data[], size_t n);
double gsl_stats_uchar_kurtosis (const unsigned char data[], size_t n);
double gsl_stats_uchar_lag1_autocorrelation (const unsigned char data[], size_t n);

double gsl_stats_uchar_variance_with_mean (const unsigned char data[], size_t n, 
					 double mean);
double gsl_stats_uchar_sd_with_mean (const unsigned char data[], size_t n, double mean);
double gsl_stats_uchar_est_variance_with_mean (const unsigned char data[], size_t n,
					     double mean);
double gsl_stats_uchar_est_sd_with_mean (const unsigned char data[], size_t n, double mean);
double gsl_stats_uchar_absdev_with_mean (const unsigned char data[], size_t n, double mean);
double gsl_stats_uchar_skew_with_mean_and_sd (const unsigned char data[], size_t n, 
					    double mean, double sd);
double gsl_stats_uchar_kurtosis_with_mean_and_sd (const unsigned char data[], size_t n,
						double mean,  double sd);
double gsl_stats_uchar_lag1_autocorrelation_with_mean (const unsigned char data[], 
						     size_t n, double mean);

double gsl_stats_uchar_pvariance (const unsigned char data1[], const unsigned char data2[],
				size_t n1, size_t n2);
double gsl_stats_uchar_ttest (const unsigned char data1[], const unsigned char data2[],
			    size_t n1, size_t n2);

unsigned char gsl_stats_uchar_max (const unsigned char data[], size_t n);
unsigned char gsl_stats_uchar_min (const unsigned char data[], size_t n);

size_t gsl_stats_uchar_max_index (const unsigned char data[], size_t n);
size_t gsl_stats_uchar_min_index (const unsigned char data[], size_t n);

void gsl_stats_uchar_sort_data (unsigned char data[], size_t n) ;

double gsl_stats_uchar_median_from_sorted_data (const unsigned char sorted_data[],
					      size_t n) ;
double gsl_stats_uchar_percentile_from_sorted_data (const unsigned char sorted_data[],
						  size_t n, const double f) ;

#endif /* GSL_STATISTICS_UCHAR_H */
