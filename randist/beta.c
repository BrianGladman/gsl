#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

double
gsl_ran_beta (const gsl_rng * r, double a, double b)
{
  /* The beta distribution has the form

     p(x) dx = (Gamma(a + b)/(Gamma(a) Gamma(b))) x^(a-1) (1-x)^(b-1) dx

     The method used here is the one described in Knuth */

  double x1 = gsl_ran_gamma (r, a) ;
  double x2 = gsl_ran_gamma (r, b) ;

  return x1/(x1 + x2) ;
}
