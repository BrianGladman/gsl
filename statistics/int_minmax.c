/* This is an automatically generated file created by
   convert_double_to_int.pl -- do not edit                */

#include <gsl_statistics_int.h>

int
gsl_stats_int_max (const int data[], unsigned int n)
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
gsl_stats_int_min (const int data[], unsigned int n)
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

unsigned int
gsl_stats_int_max_index (const int data[], unsigned int n)
{
  /* finds the index of the largest member of an integer dataset */

  int max = data[0];
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
gsl_stats_int_min_index (const int data[], unsigned int n)
{
  /* finds the index of the smallest member of an integer dataset */

  int min = data[0];
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

