#include <config.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

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

double
gsl_ran_gaussian_pdf (const double x, const double sigma)
{
  double u = x / fabs(sigma) ;
  double p = (1 / (sqrt (2 * M_PI) * fabs(sigma)) ) * exp (-u * u / 2);
  return p;
}

double
gsl_ran_ugaussian (const gsl_rng * r)
{
  return gsl_ran_gaussian (r, 1.0) ;
}

double
gsl_ran_ugaussian_pdf (const double x)
{
  return gsl_ran_gaussian_pdf (x, 1.0) ;
}

void
gsl_ran_bivariate_gaussian (const gsl_rng * r, 
			    double sigma_x, double sigma_y, double rho,
			    double *x, double *y)
{
  double u, v, r2, scale;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      u = -1 + 2 * gsl_rng_uniform (r);
      v = -1 + 2 * gsl_rng_uniform (r);

      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
  while (r2 > 1.0 || r2 == 0);

  scale = sqrt (-2.0 * log (r2) / r2);

  *x = sigma_x * u * scale;
  *y = sigma_y * (rho * u + sqrt(1 - rho*rho) * v) * scale;
}

double
gsl_ran_bivariate_gaussian_pdf (const double x, const double y, 
				const double sigma_x, const double sigma_y,
				const double rho)
{
  double u = x / sigma_x ;
  double v = y / sigma_y ;
  double c = 1 - rho*rho ;
  double p = (1 / (2 * M_PI * sigma_x * sigma_y * sqrt(c))) 
    * exp (-(u * u - 2 * rho * u * v + v * v) / (2 * c));
  return p;
}




