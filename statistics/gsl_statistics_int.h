/* This is an automatically generated file created by
   convert_double_to_int.pl -- do not edit                */

#ifndef _GSL_STATISTICS_INT_H
#define _GSL_STATISTICS_INT_H

double gsl_stats_int_mean (const int data[], unsigned int n);
double gsl_stats_int_variance (const int data[], unsigned int n);
double gsl_stats_int_sd (const int data[], unsigned int n);
double gsl_stats_int_est_variance (const int data[], unsigned int n);
double gsl_stats_int_est_sd (const int data[], unsigned int n);
double gsl_stats_int_absdev (const int data[], unsigned int n);

double gsl_stats_int_skew (const int data[], unsigned int n);
double gsl_stats_int_kurtosis (const int data[], unsigned int n);

double gsl_stats_int_variance_with_mean (const int data[], unsigned int n, double mean);
double gsl_stats_int_sd_with_mean (const int data[], unsigned int n, double mean);
double gsl_stats_int_est_variance_with_mean (const int data[], unsigned int n, double mean);
double gsl_stats_int_est_sd_with_mean (const int data[], unsigned int n, double mean);
double gsl_stats_int_absdev_with_mean (const int data[], unsigned int n, double mean);
double gsl_stats_int_skew_with_mean_and_sd (const int data[], unsigned int n, double mean, double sd);
double gsl_stats_int_kurtosis_with_mean_and_sd (const int data[], unsigned int n, double mean, double sd);

double gsl_stats_int_pvariance (const int data1[], const int data2[],
			    unsigned int n1, unsigned int n2);
double gsl_stats_int_ttest (const int data1[], const int data2[],
			unsigned int n1, unsigned int n2);

int
gsl_stats_int_max (const int data[], unsigned int n);
int
gsl_stats_int_min (const int data[], unsigned int n);

unsigned int gsl_stats_int_max_index (const int data[], unsigned int n);
unsigned int gsl_stats_int_min_index (const int data[], unsigned int n);

int gsl_stats_int_sort_data (int data[], unsigned int n) ;

double gsl_stats_int_median_from_sorted_data (const int sorted_data[],
					  unsigned int n) ;
double gsl_stats_int_percentile_from_sorted_data (const int sorted_data[],
					      unsigned int n, double f) ;

#endif /* _GSL_STATISTICS_INT_H */
