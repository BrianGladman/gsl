#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Of the two methods provided below, I think the Polar method is more
 * efficient, but only when you are actually producing two random
 * deviates.  We don't produce two, because then we'd have to save one
 * in a static variable for the next call, and that would screws up
 * re-entrant or threaded code, so we only produce one.  This makes
 * the Ratio method suddenly more appealing.  There are further tests
 * one can make if the log() is slow.  See Knuth for details */

/* Both methods pass the statistical tests; but the polar method
 * seems to be a touch faster on my home Pentium, EVEN though we
 * are only using half of the available random deviates!
 */

/* Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122 */

double
gsl_ran_gaussian (const gsl_rng * r, const double sigma)
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

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

/* Ratio method (Kinderman-Monahan); see Knuth v2, 3rd ed, p130 */
/* K+M, ACM Trans Math Software 3 (1977) 257-260. */

double
gsl_ran_gaussian_ratio_method (const gsl_rng * r, const double sigma)
{
  double u, v, x;

  do
    {
      v = gsl_rng_uniform (r);
      do
	{
	  u = gsl_rng_uniform (r);
	}
      while (u == 0);
      /* Const 1.715... = sqrt(8/e) */
      x = 1.71552776992141359295 * (v - 0.5) / u;
    }
  while (x * x > -4.0 * log (u));

  return sigma * x;
}

double
gsl_ran_gaussian_pdf (const double x, const double sigma)
{
  double u = x / fabs (sigma);
  double p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-u * u / 2);
  return p;
}

double
gsl_ran_ugaussian (const gsl_rng * r)
{
  return gsl_ran_gaussian (r, 1.0);
}

double
gsl_ran_ugaussian_ratio_method (const gsl_rng * r)
{
  return gsl_ran_gaussian_ratio_method (r, 1.0);
}

double
gsl_ran_ugaussian_pdf (const double x)
{
  return gsl_ran_gaussian_pdf (x, 1.0);
}
