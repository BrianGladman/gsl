#include <config.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Cauchy probability distribution is 

   p(x) dx = (1/(pi mu)) (1 + (x/mu)^2)^(-1) dx

   It is also known as the Lorentzian probability distribution */

double
gsl_ran_cauchy (const gsl_rng * r, const double mu)
{
  double u;
  do
    {
      u = gsl_rng_uniform (r);
    }
  while (u == 0.5);

  return mu * tan (M_PI * u);
}

double
gsl_ran_cauchy_pdf (const double x, const double mu)
{
  double u = x / mu;
  double p = (1 / (M_PI * mu)) / (1 + u * u);
  return p;
}
