#ifndef __GSL_STATISTICS_DOUBLE_H__
#define __GSL_STATISTICS_DOUBLE_H__

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

double gsl_stats_mean (const double data[], size_t stride, size_t n);
double gsl_stats_variance (const double data[], size_t stride, size_t n);
double gsl_stats_sd (const double data[], size_t stride, size_t n);
double gsl_stats_variance_with_fixed_mean (const double data[], size_t stride, size_t n, double mean);
double gsl_stats_sd_with_fixed_mean (const double data[], size_t stride, size_t n, double mean);
double gsl_stats_absdev (const double data[], size_t stride, size_t n);
double gsl_stats_skew (const double data[], size_t stride, size_t n);
double gsl_stats_kurtosis (const double data[], size_t stride, size_t n);
double gsl_stats_lag1_autocorrelation (const double data[], size_t stride, size_t n);

double gsl_stats_variance_m (const double data[], size_t stride, size_t n, double mean);
double gsl_stats_sd_m (const double data[], size_t stride, size_t n, double mean);
double gsl_stats_absdev_m (const double data[], size_t stride, size_t n, double mean);
double gsl_stats_skew_m_sd (const double data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_kurtosis_m_sd (const double data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_lag1_autocorrelation_m (const double data[], size_t stride, size_t n, double mean);

/* DEFINED FOR FLOATING POINT TYPES ONLY */

double gsl_stats_wmean (const double w[], size_t wstride, const double data[], size_t stride, size_t n);
double gsl_stats_wvariance (const double w[], size_t wstride, const double data[], size_t stride, size_t n);
double gsl_stats_wsd (const double w[], size_t wstride, const double data[], size_t stride, size_t n);
double gsl_stats_wvariance_with_fixed_mean (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double mean);
double gsl_stats_wsd_with_fixed_mean (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double mean);
double gsl_stats_wabsdev (const double w[], size_t wstride, const double data[], size_t stride, size_t n);
double gsl_stats_wskew (const double w[], size_t wstride, const double data[], size_t stride, size_t n);
double gsl_stats_wkurtosis (const double w[], size_t wstride, const double data[], size_t stride, size_t n);

double gsl_stats_wvariance_m (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double wmean);
double gsl_stats_wsd_m (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double wmean);
double gsl_stats_wabsdev_m (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double wmean);
double gsl_stats_wskew_m_sd (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double wmean, double wsd);
double gsl_stats_wkurtosis_m_sd (const double w[], size_t wstride, const double data[], size_t stride, size_t n, double wmean, double wsd);

/* END OF FLOATING POINT TYPES */

double gsl_stats_pvariance (const double data1[], size_t stride1, size_t n1, const double data2[], const size_t stride2, size_t n2);
double gsl_stats_ttest (const double data1[], size_t stride1, size_t n1, const double data2[], size_t stride2, size_t n2);

double gsl_stats_max (const double data[], size_t stride, size_t n);
double gsl_stats_min (const double data[], size_t stride, size_t n);
void gsl_stats_minmax (double * min, double * max, const double data[], size_t stride, size_t n);

size_t gsl_stats_max_index (const double data[], size_t stride, size_t n);
size_t gsl_stats_min_index (const double data[], size_t stride, size_t n);
void gsl_stats_minmax_index (size_t * min_index, size_t * max_index, const double data[], size_t stride, size_t n);

double gsl_stats_median_from_sorted_data (const double sorted_data[], size_t stride, size_t n) ;
double gsl_stats_quantile_from_sorted_data (const double sorted_data[], size_t stride, size_t n, const double f) ;

__END_DECLS

#endif /* __GSL_STATISTICS_DOUBLE_H__ */
