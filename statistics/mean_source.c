#include "source.h"

double 
FUNCTION(gsl_stats,mean) (const BASE data[], const size_t size)
{
  /* Compute the arithmetic mean of a dataset using the recurrence relation 
     mean_(n) = mean(n-1) + (data[n] - mean(n-1))/(n+1)   */

  double mean = 0;
  size_t i;

  for (i = 0; i < size; i++)
    {
      mean += (data[i] - mean)/((double)(i + 1));
    }

  return mean;
}
