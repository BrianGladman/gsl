
BASE 
FUNCTION(gsl_stats,max) (const BASE data[], const size_t n)
{
  /* finds the largest member of a dataset */

  BASE max = data[0];
  size_t i;

  for (i = 0; i < n; i++)
    {
      if (data[i] > max)
	max = data[i];
    }

  return max;
}

BASE
FUNCTION(gsl_stats,min) (const BASE data[], const size_t n)
{
  /* finds the smallest member of a dataset */

  BASE min = data[0];
  size_t i;

  for (i = 0; i < n; i++)
    {
      if (data[i] < min)
	min = data[i];
    }

  return min;

}

size_t
FUNCTION(gsl_stats,max_index) (const BASE data[], const size_t n)
{
  /* finds the index of the largest member of a dataset */
  /* if there is more than one largest value then we choose the first */

  BASE max = data[0];
  size_t i, max_index = 0;

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

size_t
FUNCTION(gsl_stats,min_index) (const BASE data[], const size_t n)
{
  /* finds the index of the smallest member of a dataset */
  /* if there is more than one largest value then we choose the first  */

  BASE min = data[0];
  size_t i, min_index = 0;

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

