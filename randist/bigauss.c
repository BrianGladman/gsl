#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* The Bivariate Gaussian probability distribution is 

   p(x,y) dxdy = (1/(2 pi sigma_x sigma_y sqrt(r))) 
                    exp(-(x^2 + y^2 - 2 r x y)/(2c)) dxdy     

*/

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
