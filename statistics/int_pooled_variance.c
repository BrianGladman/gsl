/* This is an automatically generated file created by
   convert_double_to_int.pl -- do not edit                */

#include <gsl_statistics_int.h>

double 
gsl_stats_int_pvariance (const int data1[], const int data2[],
		     const size_t n1, const size_t n2)
{
  /* Find the pooled variance of two integer datasets */

  const double var1 = gsl_stats_int_est_variance (data1, n1);
  const double var2 = gsl_stats_int_est_variance (data2, n2);

  /* calculate the pooled variance */

  const double pooled_variance = 
    (((n1 - 1) * var1) + ((n2 - 1) * var2)) / (n1 + n2 - 2);

  return pooled_variance;
}

