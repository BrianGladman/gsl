#include <gsl_statistics.h>

int
gsl_stats_imax (const int data[], unsigned int n)
{
  /* finds the largest member of an integer dataset */

  int max = data[0];
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      if (data[i] > max)
	max = data[i];
    }

  return max;
}

int
gsl_stats_imin (const int data[], unsigned int n)
{
  /* finds the smallest member of an integer dataset */

  int min = data[0];
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      if (data[i] < min)
	min = data[i];
    }

  return min;
}

