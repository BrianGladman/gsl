#include <gsl_statistics.h>

double 
gsl_stats_ipvariance (const int data1[], const int data2[],
		     unsigned int n1, unsigned int n2)
{
  /* Find the pooled variance of two integer datasets */

  double var1, var2, pooled_variance;

  pooled_variance = 0;

  /* find the variances */
  var1 = gsl_stats_iest_variance (data1, n1);
  var2 = gsl_stats_iest_variance (data2, n2);

  /* calculate the pooled variance */
  pooled_variance = (((n1 - 1) * var1) + ((n2 - 1) * var2)) / (n1 + n2 - 2);

  return pooled_variance;
}

