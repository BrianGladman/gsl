#include <math.h>
#include <gsl_statistics_int.h>

double 
gsl_stats_int_kurtosis (const int data[], unsigned int n)
{
  double mean = gsl_stats_int_mean(data, n);
  double est_sd = gsl_stats_int_est_sd_with_mean(data, n, mean);
  return gsl_stats_int_kurtosis_with_mean_and_sd(data, n, mean, est_sd);
}
    
double 
gsl_stats_int_kurtosis_with_mean_and_sd (const int data[], unsigned int n,
				     double mean, double sd)
{
  /* takes an integer dataset and finds the kurtosis */

  double sum = 0, kurtosis;
  unsigned int i;

  /* find the fourth moment the deviations, normalized by the sd */
  for (i = 0; i < n; i++)
    {
      double x = (data[i] - mean) / sd;
      sum += x * x * x * x;
    }

  kurtosis = (sum / n) - 3.0;  /* makes kurtosis zero for a Gaussian */

  return kurtosis;
}

