#ifndef __GSL_STATISTICS_INT_H__
#define __GSL_STATISTICS_INT_H__

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

double gsl_stats_int_mean (const int data[], size_t stride, size_t n);
double gsl_stats_int_variance (const int data[], size_t stride, size_t n);
double gsl_stats_int_sd (const int data[], size_t stride, size_t n);
double gsl_stats_int_variance_with_fixed_mean (const int data[], size_t stride, size_t n, double mean);
double gsl_stats_int_sd_with_fixed_mean (const int data[], size_t stride, size_t n, double mean);
double gsl_stats_int_absdev (const int data[], size_t stride, size_t n);
double gsl_stats_int_skew (const int data[], size_t stride, size_t n);
double gsl_stats_int_kurtosis (const int data[], size_t stride, size_t n);
double gsl_stats_int_lag1_autocorrelation (const int data[], size_t stride, size_t n);

double gsl_stats_int_variance_m (const int data[], size_t stride, size_t n, double mean);
double gsl_stats_int_sd_m (const int data[], size_t stride, size_t n, double mean);
double gsl_stats_int_absdev_m (const int data[], size_t stride, size_t n, double mean);
double gsl_stats_int_skew_m_sd (const int data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_int_kurtosis_m_sd (const int data[], size_t stride, size_t n, double mean, double sd);
double gsl_stats_int_lag1_autocorrelation_m (const int data[], size_t stride, size_t n, double mean);


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

__END_DECLS

#endif /* __GSL_STATISTICS_INT_H__ */
