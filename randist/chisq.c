#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

double
gsl_ran_chisq (const gsl_rng * r, unsigned int n)
{
  double sum = 0 ;
  unsigned int i = 0 ;

  for (i = 0; i < n/2; i++)
    {
      double u = gsl_rng_uniform_pos (r);
      sum += -2 * log (u) ;
    }

  if (n % 2 == 1)
    {
      double x, y, r2 ;

      do
	{
	  x = -1 + 2 * gsl_rng_uniform (r);
	  y = -1 + 2 * gsl_rng_uniform (r);
	  r2 = x * x + y * y;
	}
      while (r2 > 1.0 || r2 == 0);
      
      sum += -2.0 * x * x * log(r2) / r2 ;
    }

  /* FIXME: use gamma distribution instead */

  return sum ;
}

