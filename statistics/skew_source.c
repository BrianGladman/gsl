#include <math.h>

#include "source.h"

double 
FUNCTION(gsl_stats,skew) (const BASE data[], const size_t n)
{
  const double mean = FUNCTION(gsl_stats,mean)(data, n);
  const double est_sd = FUNCTION(gsl_stats,est_sd_with_mean)(data, n, mean);
  return FUNCTION(gsl_stats,skew_with_mean_and_sd)(data, n, mean, est_sd);
}
    
double 
FUNCTION(gsl_stats,skew_with_mean_and_sd) (const BASE data[], const size_t n,
					   const double mean, const double sd)
{
  /* takes a dataset and finds the skewness */

  double sum = 0, skew;
  size_t i;

  /* find the sum of the cubed deviations, normalized by the sd */
  for (i = 0; i < n; i++)
    {
      const double x = (data[i] - mean) / sd;
      sum += x * x * x;
    }

  skew = sum / n;

  return skew;
}

