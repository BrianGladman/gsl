#ifndef GSL_STATISTICS_ULONG_H
#define GSL_STATISTICS_ULONG_H

#include <stddef.h>

double gsl_stats_ulong_mean (const unsigned long data[], size_t stride, size_t n);
double gsl_stats_ulong_variance (const unsigned long data[], size_t stride, size_t n);
double gsl_stats_ulong_sd (const unsigned long data[], size_t stride, size_t n);
double gsl_stats_ulong_est_variance (const unsigned long data[], size_t stride, size_t n);
double gsl_stats_ulong_est_sd (const unsigned long data[], size_t stride, size_t n);
double gsl_stats_ulong_absdev (const unsigned long data[], size_t stride, size_t n);

double gsl_stats_ulong_skew (const unsigned long data[], size_t stride, size_t n);
double gsl_stats_ulong_kurtosis (const unsigned long data[], size_t stride, size_t n);
double gsl_stats_ulong_lag1_autocorrelation (const unsigned long data[], size_t stride, size_t n);

double gsl_stats_ulong_variance_with_mean (const unsigned long data[], size_t stride, size_t n, double mean);
double gsl_stats_ulong_sd_with_mean (const unsigned long data[], size_t stride, size_t n, double mean);
double gsl_stats_ulong_est_variance_with_mean (const unsigned long data[], size_t stride, size_t n, double mean);
double gsl_stats_ulong_est_sd_with_mean (const unsigned long data[], size_t stride, size_t n, double mean);
double gsl_stats_ulong_absdev_with_mean (const unsigned long data[], size_t stride, size_t n, double mean);
double gsl_stats_ulong_skew_with_mean_and_sd (const unsigned long data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_ulong_kurtosis_with_mean_and_sd (const unsigned long data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_ulong_lag1_autocorrelation_with_mean (const unsigned long data[], size_t stride, size_t n, double mean);


double gsl_stats_ulong_pvariance (const unsigned long data1[], size_t stride1, size_t n1, const unsigned long data2[], const size_t stride2, size_t n2);
double gsl_stats_ulong_ttest (const unsigned long data1[], size_t stride1, size_t n1, const unsigned long data2[], size_t stride2, size_t n2);

unsigned long gsl_stats_ulong_max (const unsigned long data[], size_t stride, size_t n);
unsigned long gsl_stats_ulong_min (const unsigned long data[], size_t stride, size_t n);
void gsl_stats_ulong_minmax (unsigned long * min, unsigned long * max, const unsigned long data[], size_t stride, size_t n);

size_t gsl_stats_ulong_max_index (const unsigned long data[], size_t stride, size_t n);
size_t gsl_stats_ulong_min_index (const unsigned long data[], size_t stride, size_t n);
void gsl_stats_ulong_minmax_index (size_t * min_index, size_t * max_index, const unsigned long data[], size_t stride, size_t n);

double gsl_stats_ulong_median_from_sorted_data (const unsigned long sorted_data[], size_t stride, size_t n) ;
double gsl_stats_ulong_quantile_from_sorted_data (const unsigned long sorted_data[], size_t stride, size_t n, const double f) ;

#endif /* GSL_STATISTICS_ULONG_H */
