#ifndef __GSL_STATISTICS_ULONG_H__
#define __GSL_STATISTICS_ULONG_H__

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

double gsl_stats_ulong_mean (const unsigned long data[], const size_t stride, const size_t n);
double gsl_stats_ulong_variance (const unsigned long data[], const size_t stride, const size_t n);
double gsl_stats_ulong_sd (const unsigned long data[], const size_t stride, const size_t n);
double gsl_stats_ulong_variance_with_fixed_mean (const unsigned long data[], const size_t stride, const size_t n, const double mean);
double gsl_stats_ulong_sd_with_fixed_mean (const unsigned long data[], const size_t stride, const size_t n, const double mean);
double gsl_stats_ulong_absdev (const unsigned long data[], const size_t stride, const size_t n);
double gsl_stats_ulong_skew (const unsigned long data[], const size_t stride, const size_t n);
double gsl_stats_ulong_kurtosis (const unsigned long data[], const size_t stride, const size_t n);
double gsl_stats_ulong_lag1_autocorrelation (const unsigned long data[], const size_t stride, const size_t n);

double gsl_stats_ulong_variance_m (const unsigned long data[], const size_t stride, const size_t n, const double mean);
double gsl_stats_ulong_sd_m (const unsigned long data[], const size_t stride, const size_t n, const double mean);
double gsl_stats_ulong_absdev_m (const unsigned long data[], const size_t stride, const size_t n, const double mean);
double gsl_stats_ulong_skew_m_sd (const unsigned long data[], const size_t stride, const size_t n, const double mean, const double sd);
double gsl_stats_ulong_kurtosis_m_sd (const unsigned long data[], const size_t stride, const size_t n, const double mean, const double sd);
double gsl_stats_ulong_lag1_autocorrelation_m (const unsigned long data[], const size_t stride, const size_t n, const double mean);


double gsl_stats_ulong_pvariance (const unsigned long data1[], const size_t stride1, const size_t n1, const unsigned long data2[], const size_t stride2, const size_t n2);
double gsl_stats_ulong_ttest (const unsigned long data1[], const size_t stride1, const size_t n1, const unsigned long data2[], const size_t stride2, const size_t n2);

unsigned long gsl_stats_ulong_max (const unsigned long data[], const size_t stride, const size_t n);
unsigned long gsl_stats_ulong_min (const unsigned long data[], const size_t stride, const size_t n);
void gsl_stats_ulong_minmax (unsigned long * min, unsigned long * max, const unsigned long data[], const size_t stride, const size_t n);

size_t gsl_stats_ulong_max_index (const unsigned long data[], const size_t stride, const size_t n);
size_t gsl_stats_ulong_min_index (const unsigned long data[], const size_t stride, const size_t n);
void gsl_stats_ulong_minmax_index (size_t * min_index, size_t * max_index, const unsigned long data[], const size_t stride, const size_t n);

double gsl_stats_ulong_median_from_sorted_data (const unsigned long sorted_data[], const size_t stride, const size_t n) ;
double gsl_stats_ulong_quantile_from_sorted_data (const unsigned long sorted_data[], const size_t stride, const size_t n, const double f) ;

__END_DECLS

#endif /* __GSL_STATISTICS_ULONG_H__ */
