#include <gsl_statistics.h>

double 
gsl_stats_max (const double data[], unsigned int n)
{
  /* finds the largest member of a dataset */

  double max = data[0];
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      if (data[i] > max)
	max = data[i];
    }

  return max;
}

double 
gsl_stats_min (const double data[], unsigned int n)
{
  /* finds the smallest member of a dataset */

  double min = data[0];
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      if (data[i] < min)
	min = data[i];
    }

  return min;

}

unsigned int
gsl_stats_max_index (const double data[], unsigned int n)
{
  /* finds the index of the largest member of a dataset */

  double max = data[0];
  unsigned int i, max_index = 0;

  for (i = 0; i < n; i++)
    {
      if (data[i] > max)
	{
	  max = data[i];
	  max_index = i ;
	}
    }

  return max_index;
}

unsigned int
gsl_stats_min_index (const double data[], unsigned int n)
{
  /* finds the index of the smallest member of a dataset */

  double min = data[0];
  unsigned int i, min_index = 0;

  for (i = 0; i < n; i++)
    {
      if (data[i] < min)
	{
	  min = data[i];
	  min_index = i ;
	}
    }

  return min_index;
}

