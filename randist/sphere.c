#include <config.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

void
gsl_ran_dir_2d (const gsl_rng * r, double *x, double *y)
{
  double u, v, d, d2;

  do
    {
      u = 2 * gsl_rng_uniform (r) - 1;
      v = 2 * gsl_rng_uniform (r) - 1;

      d2 = u * u + v * v;
    }
  while (d2 > 1.0 || d2 == 0);
  
  d = sqrt(d2) ;

  *x = u / d;
  *y = v / d;
}


void
gsl_ran_dir_3d (const gsl_rng * r, double *x, double *y, double *z)
{
  double u, v, w, d, d2;

  do
    {
      u = 2 * gsl_rng_uniform (r) - 1;
      v = 2 * gsl_rng_uniform (r) - 1;
      w = 2 * gsl_rng_uniform (r) - 1;

      d2 = u * u + v * v + w * w;
    }
  while (d2 > 1.0 || d2 == 0);

  d = sqrt(d2) ;

  *x = u / d;
  *y = v / d;
  *z = w / d;
}
