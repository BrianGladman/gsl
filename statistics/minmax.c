#include <gsl_statistics.h>

double gsl_stats_max (const double data[], unsigned int n)
{
  /* finds the largest member of a data set*/

  double max = data[0] ;
  int i;
  
  for (i = 0; i < n; i++)
    {
      if (data[i] > max) max = data[i];
    }

  return max;
}

double gsl_stats_min (const double data[], unsigned int n)
{
  /* finds the smallest member of a data set */
  
  double min = data[0];
  int i;
    
  for (i = 0; i < n; i++)
    {
      if (data[i] <min) min = data[i];
    }
  
  return min;
}


int gsl_stats_imax (const int data[], unsigned int n)
{
  /* finds the largest member of an integer data set */
  int max = data[0];
  int i;
    
  for (i = 0; i < n; i++)
    {
      if (data[i] > max) max = data[i];
    }

  return max;
}

int gsl_stats_imin (const int data[], unsigned int n)
{
  /* finds the smallest member of an integer data */
  
  int min = data[0] ;
  int i;

  for (i = 0; i < n; i++)
    {
      if (data[i] < min) min = data[i];
    }

  return min;
}

