#ifndef GSL_STATISTICS_CHAR_H
#define GSL_STATISTICS_CHAR_H

#include <stddef.h>

double gsl_stats_char_mean (const char data[], size_t stride, size_t n);
double gsl_stats_char_variance (const char data[], size_t stride, size_t n);
double gsl_stats_char_sd (const char data[], size_t stride, size_t n);
double gsl_stats_char_est_variance (const char data[], size_t stride, size_t n);
double gsl_stats_char_est_sd (const char data[], size_t stride, size_t n);
double gsl_stats_char_absdev (const char data[], size_t stride, size_t n);

double gsl_stats_char_skew (const char data[], size_t stride, size_t n);
double gsl_stats_char_kurtosis (const char data[], size_t stride, size_t n);
double gsl_stats_char_lag1_autocorrelation (const char data[], size_t stride, size_t n);

double gsl_stats_char_variance_with_mean (const char data[], size_t stride, size_t n, double mean);
double gsl_stats_char_sd_with_mean (const char data[], size_t stride, size_t n, double mean);
double gsl_stats_char_est_variance_with_mean (const char data[], size_t stride, size_t n, double mean);
double gsl_stats_char_est_sd_with_mean (const char data[], size_t stride, size_t n, double mean);
double gsl_stats_char_absdev_with_mean (const char data[], size_t stride, size_t n, double mean);
double gsl_stats_char_skew_with_mean_and_sd (const char data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_char_kurtosis_with_mean_and_sd (const char data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_char_lag1_autocorrelation_with_mean (const char data[], size_t stride, size_t n, double mean);


double gsl_stats_char_pvariance (const char data1[], size_t stride1, size_t n1, const char data2[], const size_t stride2, size_t n2);
double gsl_stats_char_ttest (const char data1[], size_t stride1, size_t n1, const char data2[], size_t stride2, size_t n2);

char gsl_stats_char_max (const char data[], size_t stride, size_t n);
char gsl_stats_char_min (const char data[], size_t stride, size_t n);
void gsl_stats_char_minmax (char * min, char * max, const char data[], size_t stride, size_t n);

size_t gsl_stats_char_max_index (const char data[], size_t stride, size_t n);
size_t gsl_stats_char_min_index (const char data[], size_t stride, size_t n);
void gsl_stats_char_minmax_index (size_t * min_index, size_t * max_index, const char data[], size_t stride, size_t n);

double gsl_stats_char_median_from_sorted_data (const char sorted_data[], size_t stride, size_t n) ;
double gsl_stats_char_quantile_from_sorted_data (const char sorted_data[], size_t stride, size_t n, const double f) ;

#endif /* GSL_STATISTICS_CHAR_H */
