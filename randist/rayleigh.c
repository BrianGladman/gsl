#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Rayleigh distribution has the form

   p(x) dx = (x / sigma^2) exp(-x^2/(2 sigma^2)) dx

   for x = 0 ... +infty */

double
gsl_ran_rayleigh (const gsl_rng * r, const double sigma)
{
  double u = gsl_rng_uniform_pos (r);

  return sigma * sqrt(-2.0 * log (u));
}

double
gsl_ran_rayleigh_pdf (const double x, const double sigma)
{
  if (x < 0)
    {
      return 0 ;
    }
  else
    {
      double u = x / sigma ;
      double p = (u / sigma) * exp(-u * u / 2.0) ;
      
      return p;
    }
}

/* The Rayleigh tail distribution has the form

   p(x) dx = (x / sigma^2) exp((a^2 - x^2)/(2 sigma^2)) dx

   for x = a ... +infty */

double
gsl_ran_rayleigh_tail (const gsl_rng * r, const double a, const double sigma)
{
  double u = gsl_rng_uniform_pos (r);

  return sqrt(a * a - 2.0 * sigma * sigma * log (u));
}

double
gsl_ran_rayleigh_tail_pdf (const double x, const double a, const double sigma)
{
  if (x < a)
    {
      return 0 ;
    }
  else
    {
      double u = x / sigma ;
      double v = a / sigma ;

      double p = (u / sigma) * exp((v + u) * (v - u) / 2.0) ;
      
      return p;
    }
}
