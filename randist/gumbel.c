#include <config.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Type I Gumbel distribution has the form,

   p(x) dx = a b exp(-(b exp(-ax) + ax)) dx

   and the Type II Gumbel distribution has the form,

   p(x) dx = b a x^-(a+1) exp(-b x^-a)) dx

 */

double
gsl_ran_gumbel1 (const gsl_rng * r, const double a, const double b)
{
  double x = gsl_rng_uniform_pos (r);

  double z = (log(b) - log(-log(x))) / a;

  return z;
}

double
gsl_ran_gumbel1_pdf (const double x, const double a, const double b)
{
  double p = a * b *  exp (-(b * exp(-a * x) + a * x));
  return p;
}

double
gsl_ran_gumbel2 (const gsl_rng * r, const double a, const double b)
{
  double x = gsl_rng_uniform_pos (r);

  double z = pow(-b / log(x), 1/a);

  return z;
}

double
gsl_ran_gumbel2_pdf (const double x, const double a, const double b)
{
  if (x <= 0)
    {
      return 0 ;
    }
  else
    {
      double p = b * a *  pow(x,-(a+1)) * exp (-b * pow(x, -a));
      return p;
    }
}




