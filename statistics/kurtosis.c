#include <math.h>
#include <gsl_statistics.h>

double 
gsl_stats_kurtosis (const double data[], const size_t n)
{
  const double mean = gsl_stats_mean(data, n);
  const double est_sd = gsl_stats_est_sd_with_mean(data, n, mean);
  return gsl_stats_kurtosis_with_mean_and_sd(data, n, mean, est_sd);
}
    
double 
gsl_stats_kurtosis_with_mean_and_sd (const double data[], const size_t n,
				     const double mean, const double sd)
{
  /* takes a dataset and finds the kurtosis */

  double sum = 0, kurtosis;
  size_t i;

  /* find the fourth moment the deviations, normalized by the sd */
  for (i = 0; i < n; i++)
    {
      const double x = (data[i] - mean) / sd;
      sum += x * x * x * x;
    }

  kurtosis = (sum / n) - 3.0;  /* makes kurtosis zero for a Gaussian */

  return kurtosis;
}

