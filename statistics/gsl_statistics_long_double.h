#ifndef GSL_STATISTICS_LONG_DOUBLE_H
#define GSL_STATISTICS_LONG_DOUBLE_H

#include <stddef.h>

double gsl_stats_long_double_mean (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_variance (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_sd (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_est_variance (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_est_sd (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_absdev (const long double data[], size_t stride, size_t n);

double gsl_stats_long_double_skew (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_kurtosis (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_lag1_autocorrelation (const long double data[], size_t stride, size_t n);

double gsl_stats_long_double_variance_with_mean (const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_sd_with_mean (const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_est_variance_with_mean (const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_est_sd_with_mean (const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_absdev_with_mean (const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_skew_with_mean_and_sd (const long double data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_long_double_kurtosis_with_mean_and_sd (const long double data[], size_t stride, size_t n, double mean,  double sd);
double gsl_stats_long_double_lag1_autocorrelation_with_mean (const long double data[], size_t stride, size_t n, double mean);

double gsl_stats_long_double_pvariance (const long double data1[], const long double data2[], size_t stride1, size_t n1, const size_t stride2, size_t n2);
double gsl_stats_long_double_ttest (const long double data1[], const long double data2[], size_t stride1, size_t n1, size_t stride2, size_t n2);

long double gsl_stats_long_double_max (const long double data[], size_t stride, size_t n);
long double gsl_stats_long_double_min (const long double data[], size_t stride, size_t n);
void gsl_stats_long_double_minmax (long double * min, long double * max, const long double data[], size_t stride, size_t n);

size_t gsl_stats_long_double_max_index (const long double data[], size_t stride, size_t n);
size_t gsl_stats_long_double_min_index (const long double data[], size_t stride, size_t n);
void gsl_stats_long_double_minmax_index (size_t * min_index, size_t * max_index, const long double data[], size_t stride, size_t n);

void gsl_stats_long_double_sort_data (long double data[], size_t stride, size_t n) ;

double gsl_stats_long_double_median_from_sorted_data (const long double sorted_data[], size_t stride, size_t n) ;
double gsl_stats_long_double_quantile_from_sorted_data (const long double sorted_data[], size_t stride, size_t n, const double f) ;

#endif /* GSL_STATISTICS_LONG_DOUBLE_H */
