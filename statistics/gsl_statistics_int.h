#ifndef GSL_STATISTICS_INT_H
#define GSL_STATISTICS_INT_H

#include <stddef.h>

double gsl_stats_int_mean (const int data[], size_t stride, size_t n);
double gsl_stats_int_variance (const int data[], size_t stride, size_t n);
double gsl_stats_int_sd (const int data[], size_t stride, size_t n);
double gsl_stats_int_est_variance (const int data[], size_t stride, size_t n);
double gsl_stats_int_est_sd (const int data[], size_t stride, size_t n);
double gsl_stats_int_absdev (const int data[], size_t stride, size_t n);

double gsl_stats_int_skew (const int data[], size_t stride, size_t n);
double gsl_stats_int_kurtosis (const int data[], size_t stride, size_t n);
double gsl_stats_int_lag1_autocorrelation (const int data[], size_t stride, size_t n);

double gsl_stats_int_variance_with_mean (const int data[], size_t stride, size_t n, double mean);
double gsl_stats_int_sd_with_mean (const int data[], size_t stride, size_t n, double mean);
double gsl_stats_int_est_variance_with_mean (const int data[], size_t stride, size_t n, double mean);
double gsl_stats_int_est_sd_with_mean (const int data[], size_t stride, size_t n, double mean);
double gsl_stats_int_absdev_with_mean (const int data[], size_t stride, size_t n, double mean);
double gsl_stats_int_skew_with_mean_and_sd (const int data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_int_kurtosis_with_mean_and_sd (const int data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_int_lag1_autocorrelation_with_mean (const int data[], size_t stride, size_t n, double mean);


double gsl_stats_int_pvariance (const int data1[], size_t stride1, size_t n1, const int data2[], const size_t stride2, size_t n2);
double gsl_stats_int_ttest (const int data1[], size_t stride1, size_t n1, const int data2[], size_t stride2, size_t n2);

int gsl_stats_int_max (const int data[], size_t stride, size_t n);
int gsl_stats_int_min (const int data[], size_t stride, size_t n);
void gsl_stats_int_minmax (int * min, int * max, const int data[], size_t stride, size_t n);

size_t gsl_stats_int_max_index (const int data[], size_t stride, size_t n);
size_t gsl_stats_int_min_index (const int data[], size_t stride, size_t n);
void gsl_stats_int_minmax_index (size_t * min_index, size_t * max_index, const int data[], size_t stride, size_t n);

double gsl_stats_int_median_from_sorted_data (const int sorted_data[], size_t stride, size_t n) ;
double gsl_stats_int_quantile_from_sorted_data (const int sorted_data[], size_t stride, size_t n, const double f) ;

#endif /* GSL_STATISTICS_INT_H */
