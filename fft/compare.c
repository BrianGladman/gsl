#include <math.h>
#include <float.h>
#include <stdio.h>

#include <compare.h>

extern int verbose;
extern int tests;
extern int passed;
extern int failed;

int
compare_complex_results (const char *name_a, const complex a[],
			 const char *name_b, const complex b[],
			 unsigned int n,
			 const double allowed_ticks)
{
  int i;
  double ticks, max_ticks = 0;
  double dr, di;
  const char *flag;

  for (i = 0; i < n; i++)
    {
      dr = b[i].real - a[i].real;
      di = b[i].imag - a[i].imag;
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
      dr = b[i].real - a[i].real;
      di = b[i].imag - a[i].imag;
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
	      a[i].real, a[i].imag, flag);
      printf ("%15s: %d  %.16f %.16f %e %s\n", name_b, i,
	      b[i].real, b[i].imag, ticks, flag);
    }

  return -1;
}
