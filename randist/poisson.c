#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

unsigned int
gsl_ran_poisson (const gsl_rng * r, double mu)
{
  double emu;
  double prod = 1.0;
  unsigned int n = 0;

  emu = exp (-mu);		/* This method works well when mu is small */

  while (prod > emu)
    {
      prod *= gsl_rng_uniform (r);
      n++;
    }

  return n - 1;			/* FIXME, could be a problem if mu = 0 (or mu < 0) ? */
}

void
gsl_ran_poisson_array (const gsl_rng * r, size_t n, unsigned int array[],
		       double mu)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      array[i] = gsl_ran_poisson (r, mu);
    }

  return;
}
