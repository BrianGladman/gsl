#include <config.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_sf.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The exponential power probability distribution is  

   p(x) dx = (1/(2 mu Gamma(1+1/a))) * exp(-|x/mu|^a) dx

   for -infty < x < infty. For a = 1 it reduces to the Laplace
   distribution. 

   The exponential power distribution is related to the gamma
   distribution by E = mu * pow(G(1/a),1/a), where E is an exponential
   power variate and G is a gamma variate.

   We use this relation for a < 1. For a >=1 we use rejection methods
   based on the laplace and gaussian distributions which should be
   faster.

   See P. R. Tadikamalla, "Random Sampling from the Exponential Power
   Distribution", Journal of the American Statistical Association,
   September 1980, Volume 75, Number 371, pages 683-686.
   
*/

double
gsl_ran_exppow (const gsl_rng * r, const double mu, const double a)
{
  if (a < 1) 
    {
      double u = gsl_rng_uniform (r) ;
      double v = gsl_ran_gamma (r, 1/a, 1.0) ;
      double z = mu * pow(v, 1/a) ;

      if (u > 0.5) 
	{
	  return z ;
	} 
      else 
	{
	  return -z ;
	}
    }
  else if (a == 1) 
    {
      /* Laplace distribution */
      return gsl_ran_laplace (r, mu) ;
    }
  else if (a < 2) 
    {
      /* Use laplace distribution for rejection method */

      double x, y, h, ratio, u ;

      /* Scale factor chosen by upper bound on ratio at a = 2 */

      double s = 1.4489 ; 
      do 
	{
	  x = gsl_ran_laplace (r, mu) ;
	  y = gsl_ran_laplace_pdf (x,mu) ;
	  h = gsl_ran_exppow_pdf (x,mu,a) ;
	  ratio = h/(s * y) ;
	  u = gsl_rng_uniform (r) ;
	} 
      while (u > ratio) ;
      
      return x ;
    }
  else if (a == 2)
    {
      /* Gaussian distribution */
      return gsl_ran_gaussian (r, mu/sqrt(2.0)) ;
    }
  else
    {
      /* Use gaussian for rejection method */

      double x, y, h, ratio, u ;
      const double sigma = mu / sqrt(2.0) ;

      /* Scale factor chosen by upper bound on ratio at a = infinity.
	 This could be improved by using a rational function
	 approximation to the bounding curve. */

      double s = 2.4091 ;  /* this is sqrt(pi) e / 2 */

      do 
	{
	  x = gsl_ran_gaussian (r, sigma) ;
	  y = gsl_ran_gaussian_pdf (x, sigma) ;
	  h = gsl_ran_exppow_pdf (x, mu, a) ;
	  ratio = h/(s * y) ;
	  u = gsl_rng_uniform (r) ;
	} 
      while (u > ratio) ;

      return x;
    }
}

double
gsl_ran_exppow_pdf (const double x, const double mu, const double a)
{
  double lg = 0, p ;
  gsl_sf_lngamma_impl (1+1/a, &lg) ;
  p = (1/(2*mu)) * exp(-pow(fabs(x/mu),a) - lg);
  return p;
}

