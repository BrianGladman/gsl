#include <gsl_statistics.h>

double 
gsl_stats_mean (const double data[], const size_t size)
{
  /* Compute the arithmetic mean of a dataset */

  double sum = 0, mean;
  size_t i;

  for (i = 0; i < size; i++)
    {
      sum += data[i];
    }

  mean = sum / size;

  return mean;
}
