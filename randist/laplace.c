#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The two-sided exponential probability distribution is  

   p(x) dx = (1/(2 mu)) * exp(-|x/mu|) dx

   for -infty < x < infty. It is also known as the Laplace distribution.  */

double
gsl_ran_laplace (const gsl_rng * r, const double mu)
{
  double u;
  do
    {
      u = 2 * gsl_rng_uniform (r) - 1.0;
    }
  while (u == 0.0);

  if (u < 0)
    {
      return mu * log (-u);
    }
  else
    {
      return -mu * log (u);
    }
}

double
gsl_ran_laplace_pdf (const double x, const double mu)
{
  double p = (1/(2*mu)) * exp (-fabs (x)/mu);
  return p;
}

