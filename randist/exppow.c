#include <math.h>
#include <gsl_math.h>
#include <gsl_sf.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The exponential power probability distribution is  

   p(x) dx = (1/(2 mu Gamma(1+1/a))) * exp(-|x/mu|^a) dx

   for -infty < x < infty. For a = 1 it reduces to the Laplace
   distribution. */

double
gsl_ran_exppow (const gsl_rng * r, const double mu, const double a)
{
  if (a < 1) 
    {
      abort () ; /* FIXME */
    }
  else if (a == 1) 
    {
      /* laplace distribution */
      return gsl_ran_laplace (r, mu) ;
    }
  else if (a < 2) 
    {
      /* use laplace distribution for rejection method */

      double x, y, h, ratio, u ;

      /* scale factor chosen by upper bound on ratio at a = 2 */

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
      /* gaussian distribution */
      return gsl_ran_gaussian (r, mu/sqrt(2.0)) ;
    }
  else
    {
      /* use gaussian for rejection method */

      double x, y, h, ratio, u ;
      const double sigma = mu / sqrt(2.0) ;

      /* scale factor chosen by upper bound on ratio at a = infinity */

      double s = 4.81804 ;
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
  double lg = 0 ;
  int status = gsl_sf_lngamma_impl (1+1/a, &lg) ;
  double p = (1/(2*mu)) * exp(-pow(fabs(x/mu),a) - lg);
  return p;
}

