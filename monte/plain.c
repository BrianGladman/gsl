/* Plain Monte-Carlo. */

/* Author: MJB */
/* RCS: $Id$ */

#define PFAC 0.1
#define TINY 1.0e-30
#define BIG 1.0e30

#define GSL_MONTE_MAX_DIM 10

#define myMAX(a,b) ((a) >= (b) ? (a) : (b))

#include <math.h>
#include <gsl_monte_plain.h>
#include <gsl_rng.h>

int gsl_monte_plain(const gsl_rng *r, const gsl_monte_f_T fun, 
		    const double* xl, const double* xu, const size_t num_dim, 
		    const size_t calls, double* avg, double* var)
{
  int status = 0;
  double sum, sum2;
  double fval;
  double x[GSL_MONTE_MAX_DIM];
  double vol;
  size_t n, i;
  
  vol = 1;
  for (i = 0; i < num_dim; i++) 
    vol *= xu[i]-xl[i];

  sum = sum2 = 0.0;
  
  for (n = 1; n <= calls; n++) {
    for (i = 0; i < num_dim; i++) 
      x[i] = xl[i] + gsl_rng_uniform(r)*(xu[i] - xl[i]);
    fval = (*fun)(x);
    sum += fval;
    sum2 += fval * fval;
  }
  *avg = vol * sum/calls;
  *var = vol * sqrt(myMAX(TINY, (sum2-sum*sum/calls)/(calls*calls)));

  return status;
}
