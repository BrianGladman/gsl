#ifndef __GSL_STATISTICS_UINT_H__
#define __GSL_STATISTICS_UINT_H__

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

double gsl_stats_uint_mean (const unsigned int data[], size_t stride, size_t n);
double gsl_stats_uint_variance (const unsigned int data[], size_t stride, size_t n);
double gsl_stats_uint_sd (const unsigned int data[], size_t stride, size_t n);
double gsl_stats_uint_variance_with_fixed_mean (const unsigned int data[], size_t stride, size_t n, double mean);
double gsl_stats_uint_sd_with_fixed_mean (const unsigned int data[], size_t stride, size_t n, double mean);
double gsl_stats_uint_absdev (const unsigned int data[], size_t stride, size_t n);
double gsl_stats_uint_skew (const unsigned int data[], size_t stride, size_t n);
double gsl_stats_uint_kurtosis (const unsigned int data[], size_t stride, size_t n);
double gsl_stats_uint_lag1_autocorrelation (const unsigned int data[], size_t stride, size_t n);

double gsl_stats_uint_variance_m (const unsigned int data[], size_t stride, size_t n, double mean);
double gsl_stats_uint_sd_m (const unsigned int data[], size_t stride, size_t n, double mean);
double gsl_stats_uint_absdev_m (const unsigned int data[], size_t stride, size_t n, double mean);
double gsl_stats_uint_skew_m_sd (const unsigned int data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_uint_kurtosis_m_sd (const unsigned int data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_uint_lag1_autocorrelation_m (const unsigned int data[], size_t stride, size_t n, double mean);


double gsl_stats_uint_pvariance (const unsigned int data1[], size_t stride1, size_t n1, const unsigned int data2[], const size_t stride2, size_t n2);
double gsl_stats_uint_ttest (const unsigned int data1[], size_t stride1, size_t n1, const unsigned int data2[], size_t stride2, size_t n2);

unsigned int gsl_stats_uint_max (const unsigned int data[], size_t stride, size_t n);
unsigned int gsl_stats_uint_min (const unsigned int data[], size_t stride, size_t n);
void gsl_stats_uint_minmax (unsigned int * min, unsigned int * max, const unsigned int data[], size_t stride, size_t n);

size_t gsl_stats_uint_max_index (const unsigned int data[], size_t stride, size_t n);
size_t gsl_stats_uint_min_index (const unsigned int data[], size_t stride, size_t n);
void gsl_stats_uint_minmax_index (size_t * min_index, size_t * max_index, const unsigned int data[], size_t stride, size_t n);

double gsl_stats_uint_median_from_sorted_data (const unsigned int sorted_data[], size_t stride, size_t n) ;
double gsl_stats_uint_quantile_from_sorted_data (const unsigned int sorted_data[], size_t stride, size_t n, const double f) ;

__END_DECLS

#endif /* __GSL_STATISTICS_UINT_H__ */
