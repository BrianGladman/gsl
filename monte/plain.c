/* Plain Monte-Carlo. */

/* Author: MJB */
/* RCS: $Id$ */

#define TINY DBL_MIN

#include <config.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_monte_plain.h>
#include <gsl_rng.h>
#include <utils.h>

int gsl_monte_plain_integrate(gsl_monte_plain_state *state, 
			      const gsl_monte_f_T fun, 
			      const double* xl, const double* xu, 
			      const size_t num_dim, 
			      const size_t calls, double* res, double* err)
{
  int status = 0;
  double sum, sum2;
  double fval;
  /*  double x[GSL_MONTE_MAX_DIM]; */
  /*   double* x = state->x; */ /* save some typing */
  double* x; 
  double vol;
  size_t n, i;

  x = gsl_monte_vector_alloc(num_dim);

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
  if ( calls > 1) {
    *err = vol * sqrt(GSL_MAX(TINY, (sum2-sum*sum/calls)/(calls*(calls-1))));
  }
  else {
    /* should't happen, validate should catch */
    *err = -1;
    status = 1;
  }
  gsl_monte_vector_free(x);

  return status;
}



gsl_monte_plain_state* gsl_monte_plain_alloc(size_t num_dim)
{
  gsl_monte_plain_state *s =  
    (gsl_monte_plain_state *) malloc(sizeof (gsl_monte_plain_state));
  
  if ( s == (gsl_monte_plain_state*) NULL) {
    GSL_ERROR_RETURN ("failed to allocate space for state struct",
                        GSL_ENOMEM, 0);
  }

  s->num_dim = num_dim;
  return s;
}

int gsl_monte_plain_validate(gsl_monte_plain_state* state,
			     const double xl[], const double xu[], 
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

  if (num_dim > state->num_dim) {
    sprintf(warning, 
	    "number of dimensions (%lu) greater than allocated size (%lu)",
	    num_dim, state->num_dim);
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

  if ( calls <= 1 ) {
    sprintf(warning, "number of calls must be greater than 1, not %lu",
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
    GSL_ERROR("Allocate state structure before calling!", GSL_EFAULT);
  }

  state->ranf = gsl_rng_alloc(gsl_rng_env_setup());

  state->init_done = 1;
  state->verbose = 1;
  return GSL_SUCCESS;
}

void gsl_monte_plain_free (gsl_monte_plain_state* s)
{
  if (s == (gsl_monte_plain_state*) NULL )
    GSL_ERROR_RETURN_NOTHING("Attempt to free null pointer", GSL_EFAULT);

  free (s);
}
