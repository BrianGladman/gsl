/* This is an automatically generated file created by
   convert_double_to_int.pl -- do not edit                */

#include <math.h>
#include <gsl_statistics_int.h>

double 
gsl_stats_int_skew (const int data[], unsigned int n)
{
  double mean = gsl_stats_int_mean(data, n);
  double est_sd = gsl_stats_int_est_sd_with_mean(data, n, mean);
  return gsl_stats_int_skew_with_mean_and_sd(data, n, mean, est_sd);
}
    
double 
gsl_stats_int_skew_with_mean_and_sd (const int data[], unsigned int n,
				 double mean, double sd)
{
  /* takes an integer dataset and finds the skewness */

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

