#include <gsl_statistics_int.h>

double 
gsl_stats_int_mean (const int data[], unsigned int size)
{
  /* Compute the arithmetic mean of an integer dataset */

  double sum = 0, mean;
  unsigned int i;

  for (i = 0; i < size; i++)
    {
      sum += data[i];
    }

  mean = sum / size;

  return mean;
}
