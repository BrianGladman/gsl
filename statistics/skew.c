#include <math.h>
#include <gsl_statistics.h>

double 
gsl_stats_skew (const double data[], unsigned int n)
{
  double mean = gsl_stats_mean(data, n);
  double est_sd = gsl_stats_est_sd_with_mean(data, n, mean);
  return gsl_stats_skew_with_mean_and_sd(data, n, mean, est_sd);
}
    
double 
gsl_stats_skew_with_mean_and_sd (const double data[], unsigned int n,
				 double mean, double sd)
{
  /* takes a dataset and finds the skewness */

  double sum = 0, skew;
  unsigned int i;

  /* find the sum of the cubed deviations, normalized by the sd */
  for (i = 0; i < n; i++)
    {
      double x = (data[i] - mean) / sd;
      sum += x * x * x;
    }

  skew = sum / n;

  return skew;
}

