#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

void
gsl_ran_dir_2d (const gsl_rng * r, double * x, double * y)
{
  double u, v, r2;

  do 
    {
      u = 2 * gsl_rng_uniform (r) - 1;
      v = 2 * gsl_rng_uniform (r) - 1;
      
      r2 = u * u + v * v ;
    }
  while (r2 > 1.0 || r2 == 0) ;

  *x = u / sqrt(r2) ;
  *y = v / sqrt(r2) ;
}


void
gsl_ran_dir_3d (const gsl_rng * r, double * x, double * y, double * z)
{
  double u, v, w, r2;

  do 
    {
      u = 2 * gsl_rng_uniform (r) - 1;
      v = 2 * gsl_rng_uniform (r) - 1;
      w = 2 * gsl_rng_uniform (r) - 1;
      
      r2 = u * u + v * v + w * w;
    }
  while (r2 > 1.0 || r2 == 0) ;

  *x = u / sqrt(r2) ;
  *y = v / sqrt(r2) ;
  *z = w / sqrt(r2) ;
}
