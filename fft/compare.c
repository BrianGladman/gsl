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
compare_complex_results (const char *name_a, const double a[],
			 const char *name_b, const double b[],
			 size_t n,
			 const double allowed_ticks)
{
  size_t i;
  double ticks, max_ticks = 0;
  double dr, di;
  const char *flag;

  for (i = 0; i < n; i++)
    {
      dr = b[2*i] - a[2*i];
      di = b[2*i+1] - a[2*i+1];
      ticks = (fabs (dr) + fabs (di)) / DBL_EPSILON;
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
      dr = b[2*i] - a[2*i];
      di = b[2*i+1] - a[2*i+1];
      ticks = (fabs (dr) + fabs (di)) / DBL_EPSILON;

      if (ticks > 1000)
	{
	  flag = "***";
	}
      else
	{
	  flag = "";
	}

      printf ("%15s: %d  %.16f %.16f %s\n", name_a, i,
	      a[2*i], a[2*i+1], flag);
      printf ("%15s: %d  %.16f %.16f %e %s\n", name_b, i,
	      b[2*i], b[2*i+1], ticks, flag);
    }

  return -1;
}
