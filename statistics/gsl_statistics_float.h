#ifndef __GSL_STATISTICS_FLOAT_H__
#define __GSL_STATISTICS_FLOAT_H__

#include <stddef.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

double gsl_stats_float_mean (const float data[], size_t stride, size_t n);
double gsl_stats_float_variance (const float data[], size_t stride, size_t n);
double gsl_stats_float_sd (const float data[], size_t stride, size_t n);
double gsl_stats_float_variance_with_fixed_mean (const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_sd_with_fixed_mean (const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_absdev (const float data[], size_t stride, size_t n);
double gsl_stats_float_skew (const float data[], size_t stride, size_t n);
double gsl_stats_float_kurtosis (const float data[], size_t stride, size_t n);
double gsl_stats_float_lag1_autocorrelation (const float data[], size_t stride, size_t n);

double gsl_stats_float_variance_m (const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_sd_m (const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_absdev_m (const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_skew_m_sd (const float data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_float_kurtosis_m_sd (const float data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_float_lag1_autocorrelation_m (const float data[], size_t stride, size_t n, double mean);

/* DEFINED FOR FLOATING POINT TYPES ONLY */

double gsl_stats_float_wmean (const float w[], size_t wstride, const float data[], size_t stride, size_t n);
double gsl_stats_float_wvariance (const float w[], size_t wstride, const float data[], size_t stride, size_t n);
double gsl_stats_float_wsd (const float w[], size_t wstride, const float data[], size_t stride, size_t n);
double gsl_stats_float_wvariance_with_fixed_mean (const float w[], size_t wstride, const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_wsd_with_fixed_mean (const float w[], size_t wstride, const float data[], size_t stride, size_t n, double mean);
double gsl_stats_float_wabsdev (const float w[], size_t wstride, const float data[], size_t stride, size_t n);
double gsl_stats_float_wskew (const float w[], size_t wstride, const float data[], size_t stride, size_t n);
double gsl_stats_float_wkurtosis (const float w[], size_t wstride, const float data[], size_t stride, size_t n);

double gsl_stats_float_wvariance_m (const float w[], size_t wstride, const float data[], size_t stride, size_t n, double wmean);
double gsl_stats_float_wsd_m (const float w[], size_t wstride, const float data[], size_t stride, size_t n, double wmean);
double gsl_stats_float_wabsdev_m (const float w[], size_t wstride, const float data[], size_t stride, size_t n, double wmean);
double gsl_stats_float_wskew_m_sd (const float w[], size_t wstride, const float data[], size_t stride, size_t n, double wmean, double wsd);
double gsl_stats_float_wkurtosis_m_sd (const float w[], size_t wstride, const float data[], size_t stride, size_t n, double wmean, double wsd);

/* END OF FLOATING POINT TYPES */

double gsl_stats_float_pvariance (const float data1[], size_t stride1, size_t n1, const float data2[], const size_t stride2, size_t n2);
double gsl_stats_float_ttest (const float data1[], size_t stride1, size_t n1, const float data2[], size_t stride2, size_t n2);

float gsl_stats_float_max (const float data[], size_t stride, size_t n);
float gsl_stats_float_min (const float data[], size_t stride, size_t n);
void gsl_stats_float_minmax (float * min, float * max, const float data[], size_t stride, size_t n);

size_t gsl_stats_float_max_index (const float data[], size_t stride, size_t n);
size_t gsl_stats_float_min_index (const float data[], size_t stride, size_t n);
void gsl_stats_float_minmax_index (size_t * min_index, size_t * max_index, const float data[], size_t stride, size_t n);

double gsl_stats_float_median_from_sorted_data (const float sorted_data[], size_t stride, size_t n) ;
double gsl_stats_float_quantile_from_sorted_data (const float sorted_data[], size_t stride, size_t n, const double f) ;

__END_DECLS

#endif /* __GSL_STATISTICS_FLOAT_H__ */
