#include <gsl_statistics.h>

double 
gsl_stats_pvariance (const double data1[], const double data2[],
		     const size_t n1, const size_t n2)
{
  /* Find the pooled variance of two datasets */

  const double var1 = gsl_stats_est_variance (data1, n1);
  const double var2 = gsl_stats_est_variance (data2, n2);

  /* calculate the pooled variance */

  const double pooled_variance = 
    (((n1 - 1) * var1) + ((n2 - 1) * var2)) / (n1 + n2 - 2);

  return pooled_variance;
}

