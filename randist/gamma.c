#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

#define GSL_LOGINFINITY 300.0
static double gamma_int (const gsl_rng * r, int a);
static double gamma_large (const gsl_rng * r, double a);
static double gamma_frac (const gsl_rng * r, double a);

/* The Gamma distribution of order a>0 is defined by:

   F(x) = \frac{1}{\Gamma(a)}\int_0^x t^{a-1} e^{-t} dt

   for x>0.  If X and Y are independent gamma-distributed random
   variables of order a and b, then X+Y has gamma distribution of
   order a+b.

   The algorithms below are from Knuth, vol 2, 2nd ed, p. 129. */

double
gsl_ran_gamma (const gsl_rng * r, double a)
{
  /* assume a > 0 */
  int na;
  na = floor (a);
  if (a == na)
    {
      return gamma_int (r, na);
    }
  else if (na == 0)
    {
      return gamma_frac (r, a);
    }
  else
    {
      return gamma_int (r, na) + gamma_frac (r, a - na);
    }
}

static double
gamma_int (const gsl_rng * r, int a)
{
  if (a < 12)
    {
      int i;
      double prod = 1.0;
      for (i = 0; i < a; ++i)
	prod *= gsl_rng_uniform (r);
      if (prod == 0)
	{
	  return GSL_LOGINFINITY;
	}
      else
	{
	  return -log (prod);
	}
    }
  else
    {
      return gamma_large (r, (double) a);
    }
}

static double
gamma_large (const gsl_rng * r, double a)
{
  /* Works only if a > 1, and is most efficient if a is large

     This algorithm, reported in Knuth, is attributed to Ahrens.  A
     faster one, we are told, can be found in: J. H. Ahrens and
     U. Dieter, Computing 12 (1974) 223-246.  */

  double sqa, x, y, v;
  sqa = sqrt (2 * a - 1);
  do
    {
      do
	{
	  y = tan (M_PI * gsl_rng_uniform (r));
	  x = sqa * y + a - 1;
	}
      while (x <= 0);
      v = gsl_rng_uniform (r);
    }
  while (v > (1 + y * y) * exp ((a - 1) * log (x / (a - 1)) - sqa * y));

  return x;
}

static double
gamma_frac (const gsl_rng * r, double a)
{
  /* This is exercise 16 from Knuth; see page 135, and the solution is
     on page 551.  */

  double p, q, x, u, v;
  p = M_E / (a + M_E);
  do
    {
      u = gsl_rng_uniform (r);
      do
	{
	  v = gsl_rng_uniform (r);
	}
      while (v == 0);
      if (u < p)
	{
	  x = exp ((1 / a) * log (v));
	  q = exp (-x);
	}
      else
	{
	  x = 1 - log (v);
	  q = exp ((a - 1) * log (x));
	}
    }
  while (gsl_rng_uniform (r) >= q);

  return x;
}
