/* This is an automatically generated file created by
   convert_double_to_int.pl -- do not edit                */

#include <math.h>
#include <gsl_statistics_int.h>

double 
gsl_stats_int_kurtosis (const int data[], const size_t n)
{
  const double mean = gsl_stats_int_mean(data, n);
  const double est_sd = gsl_stats_int_est_sd_with_mean(data, n, mean);
  return gsl_stats_int_kurtosis_with_mean_and_sd(data, n, mean, est_sd);
}
    
double 
gsl_stats_int_kurtosis_with_mean_and_sd (const int data[], const size_t n,
				     const double mean, const double sd)
{
  /* takes an integer dataset and finds the kurtosis */

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

