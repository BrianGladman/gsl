#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

static inline double gaussian (const gsl_rng * r);

static inline double
gaussian (const gsl_rng * r)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * gsl_rng_uniform (r);
      y = -1 + 2 * gsl_rng_uniform (r);

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  return y * sqrt (-2.0 * log (r2) / r2);	/* Box-Muller transform */
}
