/* MISER.  Based on the algorithm described in Numerical recipes. */

/* Author: MJB */
/* RCS: $Id$ */

#include <config.h>
#include <math.h>
#include <stdlib.h>

/* gsl headers */
#include <gsl_math.h>
/* #include <gsl_vector.h> */
/* #include <gsl_vector_int.h> */
#include <gsl_errno.h>
#include <gsl_rng.h>

#include <gsl_monte.h>
#include <gsl_monte_miser.h>
#include <gsl_monte_plain.h>
#include <utils.h>

#define SQR(a) ((a)*(a))
#define myMAX(a,b) ((a) >= (b) ? (a) : (b))
#define myMIN(a,b) ((a) <= (b) ? (a) : (b))

int gsl_monte_miser_integrate(gsl_monte_miser_state* state,
		    gsl_monte_f_T func, double xl[], double xu[], 
		    unsigned long num_dim, unsigned long calls, 
		    double *res, double *err)
{
  int status = 0;

  double *xl_tmp, *xu_tmp;

  unsigned long n, estimate_calls, calls_l, calls_r;
  unsigned long i;
  int i_bisect;
  double res_l, err_l;
  double fraction_l;
  double f;
  double xbi_l, xbi_m, xbi_r, s;

  double weight_l, weight_r;
  double var, best_var;
  double vol;

  double beta; 

  double *sigma_l, *sigma_r;
  int *hits_l, *hits_r;
  double *x, *x_mid;

  gsl_rng *r;

  status = gsl_monte_miser_validate(state, xl, xu, num_dim, calls);

  r = state->ranf;

  x = gsl_monte_vector_alloc(num_dim);

  if (calls < state->min_calls_per_bisection) {
    status = gsl_monte_plain_integrate(state->plain_state, 
			     func, xl, xu, num_dim, calls, res, err);
  }
  else {
    /* FIXME: This is bad when the estimate_calls come out less than
       or near to num_dim, because then we will get subplanes.  This
       is also an issue for min_calls. 
    */
    estimate_calls = myMAX(state->min_calls, 
			   (unsigned long)(calls*(state->estimate_frac)) );
    if (estimate_calls <= num_dim) {
      GSL_WARNING("estimate calls is close to nun_dim!", GSL_ESANITY);
    }

    x_mid = gsl_monte_vector_alloc(num_dim);

    sigma_l = gsl_monte_vector_alloc(num_dim);
    sigma_r = gsl_monte_vector_alloc(num_dim);

    vol = 1;
    for (i = 0; i < num_dim; i++) {
      /* flip a coin to bisect the integration region with some fuzz */
      s = (0.5 - gsl_rng_uniform(r) >= 0.0 ? state->dither : -state->dither ); 
      x_mid[i] = (0.5 + s)*xl[i] + (0.5 - s)*xu[i];
      vol *= xu[i] - xl[i];
    }
    
    /* 
       The idea is to chose the direction to bisect based on which will
       give the smallest total variance.  We could (and may do so later)
       use MC to compute these variances.  But the NR guys simply estimate
       the variances by finding the min and max function values 
       for each half-region for each bisection.
    */
    if (state->estimate_style == ESTIMATE_STYLE_NR ) {
      /* NR way */
      double *fmax_l, *fmax_r, *fmin_l, *fmin_r;
      fmax_l = gsl_monte_vector_alloc(num_dim);
      fmax_r = gsl_monte_vector_alloc(num_dim);
      fmin_l = gsl_monte_vector_alloc(num_dim);
      fmin_r = gsl_monte_vector_alloc(num_dim);

      for (i = 0; i < num_dim; i++) {
	fmin_l[i] = fmin_r[i] = GSL_DBL_MAX;
	fmax_l[i] = fmax_r[i] = -GSL_DBL_MAX;
      }

      for (n = 1; n <= estimate_calls; n++) {
	for (i = 0; i < num_dim; i++)
	  x[i] = xl[i] + (xu[i] - xl[i])*gsl_rng_uniform(r);
	f = (*func)(x);
	for (i = 0; i < num_dim; i++) {
	  if (x[i] <= x_mid[i]) {
	    fmin_l[i] = myMIN(fmin_l[i], f);
	    fmax_l[i] = myMAX(fmax_l[i], f);
	  }
	  else {
	    fmin_r[i] = myMIN(fmin_r[i], f);
	    fmax_r[i] = myMAX(fmax_r[i], f);
	  }
	}
      }
      for (i = 0; i < num_dim; i++) {
	if (fmax_l[i] >= fmin_l[i] && 
	    fmax_r[i] >= fmin_r[i]) {
	  sigma_l[i] = 
	    myMAX(GSL_MACH_EPS, fmax_l[i] - fmin_l[i]);
	  sigma_r[i] = 
	    myMAX(GSL_MACH_EPS, fmax_r[i] - fmin_r[i]);
	}
	else {
	  /* must be that no points landed in one of the half-regions */
	  sigma_l[i] = sigma_r[i] = -1;
	}
      }
      gsl_monte_vector_free(fmin_r);
      gsl_monte_vector_free(fmin_l);
      gsl_monte_vector_free(fmax_r);
      gsl_monte_vector_free(fmax_l);
    }
    else if (state->estimate_style == ESTIMATE_STYLE_CORRELATED_MC) {
      /* assert estimate_style = 1 */
      double *sum_l, *sum_r, *sum2_l, *sum2_r;

      hits_l = gsl_monte_vector_int_alloc(num_dim);
      hits_r = gsl_monte_vector_int_alloc(num_dim);
      sum_l = gsl_monte_vector_alloc(num_dim);
      sum2_l = gsl_monte_vector_alloc(num_dim);
      sum_r = gsl_monte_vector_alloc(num_dim);
      sum2_r = gsl_monte_vector_alloc(num_dim);

      for (i = 0; i < num_dim; i++) {
	hits_l[i] = hits_r[i] = 0;
	sum_l[i] = sum_r[i] = 0.0;
	sum2_l[i] = sum2_r[i] = 0.0;
	sigma_l[i] = sigma_r[i] = -1;
      }

      for (n = 1; n <= estimate_calls; n++) {
	for (i = 0; i < num_dim; i++)
	  x[i] = xl[i] + (xu[i] - xl[i])*gsl_rng_uniform(r);
	f = (*func)(x);
	for (i = 0; i < num_dim; i++) {
	  if (x[i] <= x_mid[i]) {
	    sum_l[i] += f;
	    sum2_l[i] = f*f;
	    hits_l[i]++;
	  }
	  else {
	    sum_r[i] += f;
	    sum_r[i] = f*f;
	    hits_r[i]++;
	  }
	}
      }
      for (i = 0; i < num_dim; i++) {
	fraction_l = (x_mid[i] - xl[i])/(xu[i] - xl[i]);
	if (hits_l[i] > 0 ) {
	  sum_l[i] /= hits_l[i];
	  sigma_l[i] = 
	    sqrt(sum2_l[i] - SQR(sum_l[i])/hits_l[i]);
	  sigma_l[i] *= fraction_l*vol/hits_l[i];
	}
	if (hits_l[i] > 0 ) {
	  sum_r[i] /= hits_r[i];
	  sigma_r[i] = 
	    sqrt(sum2_r[i] - SQR(sum_r[i])/hits_r[i]);
	  sigma_r[i] *= (1 - fraction_l)*vol/hits_r[i];;
	}
      }

      gsl_monte_vector_int_free(hits_l);
      gsl_monte_vector_int_free(hits_r);
      gsl_monte_vector_free(sum_l);
      gsl_monte_vector_free(sum2_l);
      gsl_monte_vector_free(sum_r);
      gsl_monte_vector_free(sum2_r);
    }
    else if (state->estimate_style == ESTIMATE_STYLE_MC) {
      /* FIXME */
      /* Would do complete mc estimate for each half space.  Trick is
	 how to use these results - it would be a shame to waste them.
      */
      GSL_ERROR("ESTIMATE_STYLE_MC not yet implemented", GSL_EUNSUP);
    }
    else {
      GSL_ERROR("estimate_style has three possibilities!", GSL_ESANITY);
    }

    /* Now find direction with the smallest total "variance" */
    var = 0;
    best_var = GSL_DBL_MAX;
    i_bisect = -1;
    weight_l = weight_r = 1.0;
    
    beta = 2.0/(1.0 + state->alpha);
    for (i = 0; i < num_dim; i++) {
      if (sigma_l[i] >= 0 && sigma_r[i] >= 0) {
	/* estimates are okay */
	var = pow(sigma_l[i], beta) + pow(sigma_r[i], beta);
	if (var <= best_var) {
	  best_var = var;
	  i_bisect = i;
	  weight_l = pow(sigma_l[i], beta);
	  weight_r = pow(sigma_r[i], beta);
	}
      }
      else {
	char warning[100];
	if ((sigma_l[i] < 0) && (state->verbose > 0) ) {
	  /* FIXME: Get a proper error code here */
	  sprintf(warning, "no points in left-half space(%lu)!?", i);
	  GSL_WARNING(warning, GSL_ESANITY);
	}
	if ((sigma_r[i] < 0) && (state->verbose > 0) ) {
	  /* FIXME: Get a proper error code here */
	  sprintf(warning, "no points in right-half space(%lu)!?", i);
	  GSL_WARNING(warning, GSL_ESANITY);
	}
      }
    }

    if (i_bisect < 0) {
      /* All were same, so chose direction at random */
      
      i_bisect = num_dim*gsl_rng_uniform(r); 
    }

    xbi_l = xl[i_bisect];
    xbi_m = x_mid[i_bisect];
    xbi_r = xu[i_bisect];

    /* get the actual fractional sizes of the two "halves" */
    fraction_l = fabs((xbi_m - xbi_l)/(xbi_r - xbi_l));
    calls_l = state->min_calls;
    calls_l += (unsigned long)
      (calls - estimate_calls - 2*state->min_calls)*fraction_l*weight_l
       /(fraction_l*weight_l + (1.0 - fraction_l)*weight_r);
    calls_r = calls - estimate_calls - calls_l;

    xl_tmp = sigma_l; /* reuse vectors */
    xu_tmp = sigma_r;
    for (i = 0; i < num_dim; i++) {
      xl_tmp[i] = xl[i];
      xu_tmp[i] = xu[i];
    }

    xu_tmp[i_bisect] = x_mid[i_bisect];
    status = gsl_monte_miser_integrate(state, func, xl_tmp, xu_tmp, 
			     num_dim, calls_l, &res_l, &err_l);
    xl_tmp[i_bisect] = x_mid[i_bisect];
    xu_tmp[i_bisect] = xu[i_bisect];
    status = gsl_monte_miser_integrate(state, func, xl_tmp, xu_tmp, 
			     num_dim, calls_r, res, err);

    *res += res_l;
    *err = sqrt( SQR(err_l) + SQR(*err) );   

    gsl_monte_vector_free(xl_tmp);
    gsl_monte_vector_free(xu_tmp);
    gsl_monte_vector_free(x_mid);
  }

  gsl_monte_vector_free(x);

  return status;
}

gsl_monte_miser_state* gsl_monte_miser_alloc(size_t num_dim)
{
  gsl_monte_miser_state *s =  
    (gsl_monte_miser_state *) malloc(sizeof (gsl_monte_miser_state));
  
  if ( s == (gsl_monte_miser_state*) NULL) {
    GSL_ERROR_RETURN ("failed to allocate space for miser state struct",
                        GSL_ENOMEM, 0);
  }
  s->num_dim = num_dim;

  return s;
}

int gsl_monte_miser_validate(gsl_monte_miser_state* state,
			     double xl[], double xu[], 
			     unsigned long num_dim, unsigned long calls)
{
  unsigned long i;
  char warning[100];

  if (state == (gsl_monte_miser_state*) NULL) {
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
    if (xu[i] - xl[i] > GSL_DBL_MAX) {
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

  if ( state->alpha < 0 ) {
    sprintf(warning, "alpha less than zero (currently %f) does not make sense!",
	    state->alpha);
    GSL_ERROR(warning, GSL_EINVAL);
  }

  state->check_done = 1;

  return GSL_SUCCESS;
}  

/* Set some default values and whatever */
int gsl_monte_miser_init(gsl_monte_miser_state* state)
{

  if (state == (gsl_monte_miser_state*) NULL) {
    GSL_ERROR("Allocate state structure before calling!", GSL_EFAULT);
  }

  state->min_calls = 15;
  state->min_calls_per_bisection = 60;
  state->estimate_frac =  0.1;
  state->alpha = 2.0; 
  state->dither = 0.0;
  state->verbose = 0;
  state->estimate_style = ESTIMATE_STYLE_NR;
  state->ranf = gsl_rng_alloc(gsl_rng_env_setup());

  state->plain_state = gsl_monte_plain_alloc(state->num_dim);
  gsl_monte_plain_init(state->plain_state);

  state->init_done = 1;
  return GSL_SUCCESS;
}

void gsl_monte_miser_free (gsl_monte_miser_state* s)
{
  if (s == (gsl_monte_miser_state*) NULL )
    GSL_ERROR_RETURN_NOTHING("Attempt to free null pointer", GSL_EFAULT);

  free (s);
}
