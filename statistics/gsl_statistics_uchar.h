#ifndef __GSL_STATISTICS_UCHAR_H__
#define __GSL_STATISTICS_UCHAR_H__

#include <stddef.h>

double gsl_stats_uchar_mean (const unsigned char data[], size_t stride, size_t n);
double gsl_stats_uchar_variance (const unsigned char data[], size_t stride, size_t n);
double gsl_stats_uchar_sd (const unsigned char data[], size_t stride, size_t n);
double gsl_stats_uchar_variance_with_fixed_mean (const unsigned char data[], size_t stride, size_t n, double mean);
double gsl_stats_uchar_sd_with_fixed_mean (const unsigned char data[], size_t stride, size_t n, double mean);
double gsl_stats_uchar_absdev (const unsigned char data[], size_t stride, size_t n);
double gsl_stats_uchar_skew (const unsigned char data[], size_t stride, size_t n);
double gsl_stats_uchar_kurtosis (const unsigned char data[], size_t stride, size_t n);
double gsl_stats_uchar_lag1_autocorrelation (const unsigned char data[], size_t stride, size_t n);

double gsl_stats_uchar_variance_m (const unsigned char data[], size_t stride, size_t n, double mean);
double gsl_stats_uchar_sd_m (const unsigned char data[], size_t stride, size_t n, double mean);
double gsl_stats_uchar_absdev_m (const unsigned char data[], size_t stride, size_t n, double mean);
double gsl_stats_uchar_skew_m_sd (const unsigned char data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_uchar_kurtosis_m_sd (const unsigned char data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_uchar_lag1_autocorrelation_m (const unsigned char data[], size_t stride, size_t n, double mean);


double gsl_stats_uchar_pvariance (const unsigned char data1[], size_t stride1, size_t n1, const unsigned char data2[], const size_t stride2, size_t n2);
double gsl_stats_uchar_ttest (const unsigned char data1[], size_t stride1, size_t n1, const unsigned char data2[], size_t stride2, size_t n2);

unsigned char gsl_stats_uchar_max (const unsigned char data[], size_t stride, size_t n);
unsigned char gsl_stats_uchar_min (const unsigned char data[], size_t stride, size_t n);
void gsl_stats_uchar_minmax (unsigned char * min, unsigned char * max, const unsigned char data[], size_t stride, size_t n);

size_t gsl_stats_uchar_max_index (const unsigned char data[], size_t stride, size_t n);
size_t gsl_stats_uchar_min_index (const unsigned char data[], size_t stride, size_t n);
void gsl_stats_uchar_minmax_index (size_t * min_index, size_t * max_index, const unsigned char data[], size_t stride, size_t n);

double gsl_stats_uchar_median_from_sorted_data (const unsigned char sorted_data[], size_t stride, size_t n) ;
double gsl_stats_uchar_quantile_from_sorted_data (const unsigned char sorted_data[], size_t stride, size_t n, const double f) ;

#endif /* __GSL_STATISTICS_UCHAR_H__ */
