#ifndef GSL_STATISTICS_USHORT_H
#define GSL_STATISTICS_USHORT_H

#include <stddef.h>

double gsl_stats_ushort_mean (const unsigned short data[], size_t stride, size_t n);
double gsl_stats_ushort_variance (const unsigned short data[], size_t stride, size_t n);
double gsl_stats_ushort_sd (const unsigned short data[], size_t stride, size_t n);
double gsl_stats_ushort_est_variance (const unsigned short data[], size_t stride, size_t n);
double gsl_stats_ushort_est_sd (const unsigned short data[], size_t stride, size_t n);
double gsl_stats_ushort_absdev (const unsigned short data[], size_t stride, size_t n);

double gsl_stats_ushort_skew (const unsigned short data[], size_t stride, size_t n);
double gsl_stats_ushort_kurtosis (const unsigned short data[], size_t stride, size_t n);
double gsl_stats_ushort_lag1_autocorrelation (const unsigned short data[], size_t stride, size_t n);

double gsl_stats_ushort_variance_with_mean (const unsigned short data[], size_t stride, size_t n, double mean);
double gsl_stats_ushort_sd_with_mean (const unsigned short data[], size_t stride, size_t n, double mean);
double gsl_stats_ushort_est_variance_with_mean (const unsigned short data[], size_t stride, size_t n, double mean);
double gsl_stats_ushort_est_sd_with_mean (const unsigned short data[], size_t stride, size_t n, double mean);
double gsl_stats_ushort_absdev_with_mean (const unsigned short data[], size_t stride, size_t n, double mean);
double gsl_stats_ushort_skew_with_mean_and_sd (const unsigned short data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_ushort_kurtosis_with_mean_and_sd (const unsigned short data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_ushort_lag1_autocorrelation_with_mean (const unsigned short data[], size_t stride, size_t n, double mean);


double gsl_stats_ushort_pvariance (const unsigned short data1[], size_t stride1, size_t n1, const unsigned short data2[], const size_t stride2, size_t n2);
double gsl_stats_ushort_ttest (const unsigned short data1[], size_t stride1, size_t n1, const unsigned short data2[], size_t stride2, size_t n2);

unsigned short gsl_stats_ushort_max (const unsigned short data[], size_t stride, size_t n);
unsigned short gsl_stats_ushort_min (const unsigned short data[], size_t stride, size_t n);
void gsl_stats_ushort_minmax (unsigned short * min, unsigned short * max, const unsigned short data[], size_t stride, size_t n);

size_t gsl_stats_ushort_max_index (const unsigned short data[], size_t stride, size_t n);
size_t gsl_stats_ushort_min_index (const unsigned short data[], size_t stride, size_t n);
void gsl_stats_ushort_minmax_index (size_t * min_index, size_t * max_index, const unsigned short data[], size_t stride, size_t n);

double gsl_stats_ushort_median_from_sorted_data (const unsigned short sorted_data[], size_t stride, size_t n) ;
double gsl_stats_ushort_quantile_from_sorted_data (const unsigned short sorted_data[], size_t stride, size_t n, const double f) ;

#endif /* GSL_STATISTICS_USHORT_H */
