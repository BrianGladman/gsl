#include <math.h>
#include <gsl_math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

double
gsl_ran_cauchy (const gsl_rng * r, double mu)
{
  double u;
  do
    {
      u = gsl_rng_uniform (r);
    }
  while (u == 0.5);

  return mu * tan (M_PI * u);
}

/* The Cauchy probability distribution is 

   p(x) = (1/pi) (1 + (x/mu)^2)^(-1)

   It is also known as the Lorentzian probability distribution */
