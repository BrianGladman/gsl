#ifndef GSL_STATISTICS_LONG_DOUBLE_H
#define GSL_STATISTICS_LONG_DOUBLE_H

#include <stddef.h>

double gsl_stats_long_double_mean (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_variance (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_sd (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_variance_with_fixed_mean (const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_sd_with_fixed_mean (const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_absdev (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_skew (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_kurtosis (const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_lag1_autocorrelation (const long double data[], size_t stride, size_t n);

double gsl_stats_long_double_variance_m (const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_sd_m (const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_absdev_m (const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_skew_m_sd (const long double data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_long_double_kurtosis_m_sd (const long double data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_long_double_lag1_autocorrelation_m (const long double data[], size_t stride, size_t n, double mean);

/* DEFINED FOR FLOATING POINT TYPES ONLY */

double gsl_stats_long_double_wmean (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_wvariance (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_wsd (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_wvariance_with_fixed_mean (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_wsd_with_fixed_mean (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n, double mean);
double gsl_stats_long_double_wabsdev (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_wskew (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n);
double gsl_stats_long_double_wkurtosis (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n);

double gsl_stats_long_double_wvariance_m (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n, double wmean);
double gsl_stats_long_double_wsd_m (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n, double wmean);
double gsl_stats_long_double_wabsdev_m (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n, double wmean);
double gsl_stats_long_double_wskew_m_sd (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n, double wmean, double wsd);
double gsl_stats_long_double_wkurtosis_m_sd (const long double w[], size_t wstride, const long double data[], size_t stride, size_t n, double wmean, double wsd);

/* END OF FLOATING POINT TYPES */

double gsl_stats_long_double_pvariance (const long double data1[], size_t stride1, size_t n1, const long double data2[], const size_t stride2, size_t n2);
double gsl_stats_long_double_ttest (const long double data1[], size_t stride1, size_t n1, const long double data2[], size_t stride2, size_t n2);

long double gsl_stats_long_double_max (const long double data[], size_t stride, size_t n);
long double gsl_stats_long_double_min (const long double data[], size_t stride, size_t n);
void gsl_stats_long_double_minmax (long double * min, long double * max, const long double data[], size_t stride, size_t n);

size_t gsl_stats_long_double_max_index (const long double data[], size_t stride, size_t n);
size_t gsl_stats_long_double_min_index (const long double data[], size_t stride, size_t n);
void gsl_stats_long_double_minmax_index (size_t * min_index, size_t * max_index, const long double data[], size_t stride, size_t n);

double gsl_stats_long_double_median_from_sorted_data (const long double sorted_data[], size_t stride, size_t n) ;
double gsl_stats_long_double_quantile_from_sorted_data (const long double sorted_data[], size_t stride, size_t n, const double f) ;

#endif /* GSL_STATISTICS_LONG_DOUBLE_H */
