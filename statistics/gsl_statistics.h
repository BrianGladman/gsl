#ifndef _GSL_STATISTICS_H
#define _GSL_STATISTICS_H

double gsl_stats_mean (const double data[], unsigned int n);
double gsl_stats_variance (const double data[], unsigned int n);
double gsl_stats_stddev (const double data[], unsigned int n);
double gsl_stats_est_variance (const double data[], unsigned int n);
double gsl_stats_est_stddev (const double data[], unsigned int n);
double gsl_stats_pvariance (const double data1[], const double data2[], 
			    unsigned int n1, unsigned int n2);
double gsl_stats_ttest (const double data1[], const double data2[],
			 unsigned int n1, unsigned int n2);

double gsl_stats_max (const double data[], unsigned int n);
double gsl_stats_min (const double data[], unsigned int n);

double gsl_stats_imean (const int data[], unsigned int n);
double gsl_stats_ivariance (const int data[], unsigned int n);
double gsl_stats_istddev (const int data[], unsigned int n);
double gsl_stats_iest_variance (const int data[], unsigned int n);
double gsl_stats_iest_stddev (const int data[], unsigned int n);
double gsl_stats_ipvariance (const int data1[], const int data2[], unsigned int n1, unsigned int n2);
double gsl_stats_ittest (const int data1[], const int data2[], unsigned int n1, unsigned int n2);
int gsl_stats_imax (const int data[], unsigned int n);
int gsl_stats_imin (const int data[], unsigned int n);

#endif /* _GSL_STATISTICS_H */
