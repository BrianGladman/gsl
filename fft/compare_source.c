#include <config.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "compare.h"

extern int verbose;
extern int tests;
extern int passed;
extern int failed;

int
FUNCTION(compare_complex,results) (const char *name_a, const BASE a[],
				   const char *name_b, const BASE b[],
				   size_t n, size_t stride,
				   const double allowed_ticks)
{
  size_t i;
  double ticks, max_ticks = 0;
  double dr, di;
  const char *flag;

  for (i = 0; i < n; i++)
    {
      dr = b[2*stride*i] - a[2*stride*i];
      di = b[2*stride*i+1] - a[2*stride*i+1];
      ticks = (fabs (dr) + fabs (di)) / BASE_EPSILON;
      if (ticks > max_ticks)
	{
	  max_ticks = ticks;
	}
    }

  if (max_ticks < allowed_ticks)
    {
      return 0;
    }

  printf ("\n%s vs %s : max_ticks = %f\n", name_a, name_b, max_ticks);

  for (i = 0; i < n; i++)
    {
      dr = b[2*stride*i] - a[2*stride*i];
      di = b[2*stride*i+1] - a[2*stride*i+1];
      ticks = (fabs (dr) + fabs (di)) / BASE_EPSILON;

      if (ticks > 1000)
	{
	  flag = "***";
	}
      else
	{
	  flag = "";
	}

      printf ("%15s: %d  %.16f %.16f %s\n", name_a, i,
	      a[2*stride*i], a[2*stride*i+1], flag);
      printf ("%15s: %d  %.16f %.16f %e %s\n", name_b, i,
	      b[2*stride*i], b[2*stride*i+1], ticks, flag);
    }

  return -1;
}
