#ifndef GSL_STATISTICS_FLOAT_H
#define GSL_STATISTICS_FLOAT_H

#include <stddef.h>

double gsl_stats_float_mean (const float data[], size_t stride, size_t n);
double gsl_stats_float_variance (const float data[], size_t stride, size_t n);
double gsl_stats_float_sd (const float data[], size_t stride, size_t n);
double gsl_stats_float_est_variance (const float data[], size_t stride, size_t n);
double gsl_stats_float_est_sd (const float data[], size_t stride, size_t n);
double gsl_stats_float_absdev (const float data[], size_t stride, size_t n);

double gsl_stats_float_skew (const float data[], size_t stride, size_t n);
double gsl_stats_float_kurtosis (const float data[], size_t stride, size_t n);
double gsl_stats_float_lag1_autocorrelation (const float data[], size_t stride, size_t n);

double gsl_stats_float_variance_with_mean (const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_sd_with_mean (const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_est_variance_with_mean (const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_est_sd_with_mean (const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_absdev_with_mean (const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_skew_with_mean_and_sd (const float data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_float_kurtosis_with_mean_and_sd (const float data[], size_t stride, size_t n, double mean,  double sd);
double gsl_stats_float_lag1_autocorrelation_with_mean (const float data[], size_t stride, size_t n, double mean);

double gsl_stats_float_pvariance (const float data1[], const float data2[], size_t stride1, size_t n1, const size_t stride2, size_t n2);
double gsl_stats_float_ttest (const float data1[], const float data2[], size_t stride1, size_t n1, size_t stride2, size_t n2);

float gsl_stats_float_max (const float data[], size_t stride, size_t n);
float gsl_stats_float_min (const float data[], size_t stride, size_t n);
void gsl_stats_float_minmax (float * min, float * max, const float data[], size_t stride, size_t n);

size_t gsl_stats_float_max_index (const float data[], size_t stride, size_t n);
size_t gsl_stats_float_min_index (const float data[], size_t stride, size_t n);
void gsl_stats_float_minmax_index (size_t * min_index, size_t * max_index, const float data[], size_t stride, size_t n);

void gsl_stats_float_sort_data (float data[], size_t stride, size_t n) ;

double gsl_stats_float_median_from_sorted_data (const float sorted_data[], size_t stride, size_t n) ;
double gsl_stats_float_quantile_from_sorted_data (const float sorted_data[], size_t stride, size_t n, const double f) ;

#endif /* GSL_STATISTICS_FLOAT_H */
