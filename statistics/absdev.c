#include <math.h>
#include <gsl_statistics.h>

double 
gsl_stats_absdev (const double data[], const size_t n)
{
  const double mean = gsl_stats_mean(data, n);
  return gsl_stats_absdev_with_mean(data, n, mean);
}
    
double 
gsl_stats_absdev_with_mean (const double data[], 
			    const size_t n, 
			    const double mean)
{
  /* takes a dataset and finds the absolute deviation */

  double sum = 0, absdev;
  size_t i;

  /* find the sum of the absolute deviations */
  for (i = 0; i < n; i++)
    {
      const double delta = fabs(data[i] - mean);
      sum += delta;
    }

  absdev = sum / n;

  return absdev;
}

