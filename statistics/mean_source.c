#include "source.h"
#include <stdio.h>
double 
FUNCTION(gsl_stats,mean) (const BASE data[], const size_t size)
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
