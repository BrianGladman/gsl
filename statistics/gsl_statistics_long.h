#ifndef GSL_STATISTICS_LONG_H
#define GSL_STATISTICS_LONG_H

#include <stddef.h>

double gsl_stats_long_mean (const long data[], size_t stride, size_t n);
double gsl_stats_long_variance (const long data[], size_t stride, size_t n);
double gsl_stats_long_sd (const long data[], size_t stride, size_t n);
double gsl_stats_long_est_variance (const long data[], size_t stride, size_t n);
double gsl_stats_long_est_sd (const long data[], size_t stride, size_t n);
double gsl_stats_long_absdev (const long data[], size_t stride, size_t n);

double gsl_stats_long_skew (const long data[], size_t stride, size_t n);
double gsl_stats_long_kurtosis (const long data[], size_t stride, size_t n);
double gsl_stats_long_lag1_autocorrelation (const long data[], size_t stride, size_t n);

double gsl_stats_long_variance_with_mean (const long data[], size_t stride, size_t n, double mean);
double gsl_stats_long_sd_with_mean (const long data[], size_t stride, size_t n, double mean);
double gsl_stats_long_est_variance_with_mean (const long data[], size_t stride, size_t n, double mean);
double gsl_stats_long_est_sd_with_mean (const long data[], size_t stride, size_t n, double mean);
double gsl_stats_long_absdev_with_mean (const long data[], size_t stride, size_t n, double mean);
double gsl_stats_long_skew_with_mean_and_sd (const long data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_long_kurtosis_with_mean_and_sd (const long data[], size_t stride, size_t n, double mean,  double sd);
double gsl_stats_long_lag1_autocorrelation_with_mean (const long data[], size_t stride, size_t n, double mean);

double gsl_stats_long_pvariance (const long data1[], const long data2[], size_t stride1, size_t n1, const size_t stride2, size_t n2);
double gsl_stats_long_ttest (const long data1[], const long data2[], size_t stride1, size_t n1, size_t stride2, size_t n2);

long gsl_stats_long_max (const long data[], size_t stride, size_t n);
long gsl_stats_long_min (const long data[], size_t stride, size_t n);
void gsl_stats_long_minmax (long * min, long * max, const long data[], size_t stride, size_t n);

size_t gsl_stats_long_max_index (const long data[], size_t stride, size_t n);
size_t gsl_stats_long_min_index (const long data[], size_t stride, size_t n);
void gsl_stats_long_minmax_index (size_t * min_index, size_t * max_index, const long data[], size_t stride, size_t n);

void gsl_stats_long_sort_data (long data[], size_t stride, size_t n) ;

double gsl_stats_long_median_from_sorted_data (const long sorted_data[], size_t stride, size_t n) ;
double gsl_stats_long_quantile_from_sorted_data (const long sorted_data[], size_t stride, size_t n, const double f) ;

#endif /* GSL_STATISTICS_LONG_H */
