/* Plain Monte-Carlo. */

/* Author: MJB */
/* RCS: $Id$ */

#define TINY DBL_MIN

#define GSL_MONTE_MAX_DIM 10

#define myMAX(a,b) ((a) >= (b) ? (a) : (b))

#include <math.h>
#include <gsl_math.h>
#include <gsl_monte_plain.h>
#include <gsl_rng.h>

int gsl_monte_plain(const gsl_monte_plain_state *state, 
		    const gsl_monte_f_T fun, 
		    const double* xl, const double* xu, const size_t num_dim, 
		    const size_t calls, double* res, double* err)
{
  int status = 0;
  double sum, sum2;
  double fval;
  double x[GSL_MONTE_MAX_DIM];
  double vol;
  size_t n, i;

  status = gsl_monte_plain_validate(state, xl, xu, num_dim, calls);

  vol = 1;
  for (i = 0; i < num_dim; i++) 
    vol *= xu[i]-xl[i];

  sum = sum2 = 0.0;
  
  for (n = 1; n <= calls; n++) {
    for (i = 0; i < num_dim; i++) 
      x[i] = xl[i] + gsl_rng_uniform(state->ranf)*(xu[i] - xl[i]);
    fval = (*fun)(x);
    sum += fval;
    sum2 += fval * fval;
  }
  *res = vol * sum/calls;
  *err = vol * sqrt(myMAX(TINY, (sum2-sum*sum/calls)/(calls*calls)));

  return status;
}

gsl_monte_plain_state* gsl_monte_plain_alloc(void)
{
  gsl_monte_plain_state *s =  
    (gsl_monte_plain_state *) malloc(sizeof (gsl_monte_plain_state));
  
  if ( s == (gsl_monte_plain_state*) NULL) {
    GSL_ERROR_RETURN ("failed to allocate space for miser state struct",
                        GSL_ENOMEM, 0);
  }

  return s;
}

int gsl_monte_plain_validate(gsl_monte_plain_state* state,
			     double xl[], double xu[], 
			     unsigned long num_dim, unsigned long calls)
{
  unsigned long i;
  char warning[100];

  if (state == (gsl_monte_plain_state*) NULL) {
    GSL_ERROR("Allocate state structure before calling!", GSL_EINVAL);

  }

  if (state->check_done) 
    return GSL_SUCCESS;
    
  if (num_dim <= 0) {
    sprintf(warning, "number of dimensions must be greater than zero, not %lu",
	    num_dim);
    GSL_ERROR(warning, GSL_EINVAL);
  }
  
  for (i=0; i < num_dim; i++ ) {
    if (xu[i] - xl[i] <= 0 ) {
      sprintf(warning, "xu[%lu] must be greater than xu[%lu]", i, i);
    GSL_ERROR(warning, GSL_EINVAL);
    }
    if (xu[i] - xl[i] > DBL_MAX) {
      sprintf(warning, 
	      "Range of integration is too large for cord %lu, please rescale", 
	      i);
      GSL_ERROR(warning, GSL_EINVAL);
    }
  }

  if ( calls <= 0 ) {
    sprintf(warning, "number of calls must be greater than zero, not %lu",
	    calls);
    GSL_ERROR(warning, GSL_EINVAL);
  }
  
  state->check_done = 1;

  return GSL_SUCCESS;
}  

/* Set some default values and whatever */
int gsl_monte_plain_init(gsl_monte_plain_state* state)
{

  if (state == (gsl_monte_plain_state*) NULL) {
    GSL_ERROR("Allocate state structure before calling!", GSL_EINVAL);
  }

  state->ranf = gsl_rng_alloc(gsl_rng_env_setup());

  state->init_done = 1;
  state->verbose = 1;
  return GSL_SUCCESS;
}

void gsl_monte_plain_free (gsl_monte_plain_state* s)
{
  free (s);
}
