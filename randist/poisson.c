#include <config.h>
#include <math.h>
#include <gsl_sf.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The poisson distribution has the form

   p(n) = (mu^n / n!) exp(-mu) 

   for n = 0, 1, 2, ... . The method used here is the one from Knuth. */

unsigned int
gsl_ran_poisson (const gsl_rng * r, double mu)
{
  double emu;
  double prod = 1.0;
  unsigned int n = 0;

  while (mu > 10)
    {
      unsigned int m = mu * (7.0 / 8.0);

      double X = gsl_ran_gamma_int (r, m);

      if (X >= mu)
	{
	  return gsl_ran_binomial (r, mu / X, m - 1);
	}
      else
	{
	  n += m;
	  mu -= X;
	}
    }

  /* This following method works well when mu is small */

  emu = exp (-mu);

  do
    {
      prod *= gsl_rng_uniform (r);
      n++;
    }
  while (prod > emu);

  return n - 1;

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

double
gsl_ran_poisson_pdf (const unsigned int n, const double mu)
{
  double lf = gsl_sf_lnfact (n); 
  double p = exp (log (mu) * n - lf - mu);
  return p;
}
