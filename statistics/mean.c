#include <gsl_statistics.h>

double gsl_stats_mean (const double data[], unsigned int size)
{
  /* Compute the arithmetic mean of a data set */
  
  double sum = 0, mean;
  int i;
  
  for (i = 0; i < size; i++)
    {
      sum += data[i];
    }

  mean = sum / size;

  return mean;
}

double gsl_stats_imean (const int data[], unsigned int size)
{
  /* Compute the arithmetic mean of an integer data set */
  
  double sum = 0, mean;
  int i;
  
  for (i = 0; i < size; i++)
    {
      sum += data[i];
    }

  mean = sum / size;

  return mean;
}


