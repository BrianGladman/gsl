/* directional_minimize.c -- wrapper for calling the one dimensional
   minimization algorithm for vector functions */

#include <gsl_multimin.h>
#include <gsl_blas_types.h>
#include <gsl_blas.h>

void
gsl_multimin_compute_evaluation_point(gsl_vector *evaluation_point,
				      const gsl_vector *starting_point,
				      double x,
				      const gsl_vector *direction) 
{
  /* evaluation_point = starting_point + x * direction */
  size_t i;
  for(i=0;i<starting_point->size;i++) 
    {
      gsl_vector_set(evaluation_point,i,
		     gsl_vector_get(starting_point,i)
		     + x * gsl_vector_get(direction,i));
    }
}

void
gsl_multimin_compute_ep(gsl_multimin_to_single_params *params,double x) 
{
  gsl_multimin_compute_evaluation_point(params->evaluation_point,
					params->starting_point,
					x,
					params->direction);
}

double gsl_multimin_to_single_eval(double x, void * oparams)
{
  gsl_multimin_to_single *params = (gsl_multimin_to_single *)oparams;
  
  gsl_multimin_compute_ep(params->params,x);
  return GSL_MULTIMIN_FN_EVAL(params->f,params->params->evaluation_point);
}

double gsl_multimin_to_single_eval_fdf(double x, void * oparams)
{
  gsl_multimin_to_single_fdf *params = (gsl_multimin_to_single_fdf *)oparams;
  double value;
  
  gsl_multimin_compute_ep(params->params,x);
  value = GSL_MULTIMIN_FN_EVAL_F(params->fdf,params->params->evaluation_point);
  /*  printf("x "); gsl_vector_fprintf (stdout,params->params->evaluation_point , "%g"); 
      printf("f(x) %g\n",value);*/
  return value;
}

gsl_multimin_to_single_params *
gsl_multimin_to_single_params_alloc(gsl_vector * starting_point,
				    gsl_vector * direction) {
  gsl_multimin_to_single_params *p;
  int status;

  p = (gsl_multimin_to_single_params *)malloc (sizeof (gsl_multimin_to_single_params));

  if (p == 0) {
    GSL_ERROR_RETURN ("failed to allocate space for multimin wrapper params",
		      GSL_ENOMEM, 0);
  }

  p->evaluation_point = gsl_vector_calloc (starting_point->size);

  if (p->evaluation_point == 0) {
    free(p);

    GSL_ERROR_RETURN ("failed to allocate space for evaluation_point", GSL_ENOMEM, 0);
  }

  status = gsl_multimin_to_single_params_set(p,starting_point,direction);

  if (status != GSL_SUCCESS) {
    gsl_vector_free (p->evaluation_point);
    free(p);
    
    GSL_ERROR_RETURN ("failed to set multimin wrapper params", status, 0);
  }

  return p;
}

int
gsl_multimin_to_single_params_set(gsl_multimin_to_single_params *p,
				  gsl_vector * starting_point,
				  gsl_vector * direction) {
  if (p->evaluation_point->size != starting_point->size 
      || p->evaluation_point->size != direction->size) {
    GSL_ERROR_RETURN ("vector length not compatible with function", 
		      GSL_EBADLEN, 0);
  }
  p->starting_point = starting_point;
  p->direction = direction;
  return GSL_SUCCESS;
}

void
gsl_multimin_to_single_params_free(gsl_multimin_to_single_params *p) {
   gsl_vector_free (p->evaluation_point);
   free(p);
}

gsl_multimin_to_single *
gsl_multimin_to_single_alloc(const gsl_multimin_function *f,
			     gsl_vector * starting_point,
			     gsl_vector * direction) {
  gsl_multimin_to_single *wrapper;
  const size_t n = f->n ;
  int status;

  if (starting_point->size != n || direction->size != n) {
    GSL_ERROR_RETURN ("vector length not compatible with function", 
		      GSL_EBADLEN, 0);
  }
  wrapper = (gsl_multimin_to_single *)malloc (sizeof (gsl_multimin_to_single));

  if (wrapper == 0) {
    GSL_ERROR_RETURN ("failed to allocate space for multimin wrapper struct",
		      GSL_ENOMEM, 0);
  }
  wrapper->params = gsl_multimin_to_single_params_alloc(starting_point,
							direction);
  if (wrapper->params == 0) {
    free(wrapper);

    GSL_ERROR_RETURN ("failed to allocate multimin wrapper params", GSL_ENOMEM, 0);
  }
  wrapper->f = f;
  return wrapper;
}

int
gsl_multimin_to_single_set(gsl_multimin_to_single *w,
			   const gsl_multimin_function *f,
			   gsl_vector * starting_point,
			   gsl_vector * direction) {
  w->f = f;
  return gsl_multimin_to_single_params_set(w->params,starting_point,direction);
}

void
gsl_multimin_to_single_free(gsl_multimin_to_single *w) {
  gsl_multimin_to_single_params_free(w->params);
  free(w);
}

gsl_function *
gsl_multimin_to_single_function_alloc(gsl_multimin_to_single *w) {
  gsl_function *f;
  
  f = (gsl_function*)malloc (sizeof (gsl_function));
  if (f == 0) {
     GSL_ERROR_RETURN ("failed to allocate space for function struct",
		       GSL_ENOMEM, 0);
  }
  f->params = (void *)w;
  f->function = gsl_multimin_to_single_eval;

  return f;
}

void
gsl_multimin_to_single_function_free(gsl_function *f) {
  free (f);
}
gsl_multimin_to_single_fdf *
gsl_multimin_to_single_fdf_alloc(const gsl_multimin_function_fdf *fdf,
				 gsl_vector * starting_point,
				 gsl_vector * direction) {
  gsl_multimin_to_single_fdf *wrapper;
  const size_t n = fdf->n ;
  int status;

  if (starting_point->size != n || direction->size != n) {
    GSL_ERROR_RETURN ("vector length not compatible with function", 
		      GSL_EBADLEN, 0);
  }
  wrapper = (gsl_multimin_to_single_fdf *)malloc (sizeof (gsl_multimin_to_single_fdf));

  if (wrapper == 0) {
    GSL_ERROR_RETURN ("failed to allocate space for multimin wrapper struct",
		      GSL_ENOMEM, 0);
  }
  wrapper->params = gsl_multimin_to_single_params_alloc(starting_point,
							direction);
  if (wrapper->params == 0) {
    free(wrapper);

    GSL_ERROR_RETURN ("failed to allocate multimin wrapper params", GSL_ENOMEM, 0);
  }
  wrapper->fdf = fdf;
  return wrapper;
}

int
gsl_multimin_to_single_fdf_set(gsl_multimin_to_single_fdf *w,
			       const gsl_multimin_function_fdf *fdf,
			       gsl_vector * starting_point,
			       gsl_vector * direction) {
  w->fdf = fdf;
  return gsl_multimin_to_single_params_set(w->params,starting_point,direction);
}

void
gsl_multimin_to_single_fdf_free(gsl_multimin_to_single_fdf *w) {
  gsl_multimin_to_single_params_free(w->params);
  free(w);
}

gsl_function *
gsl_multimin_to_single_function_fdf_alloc(gsl_multimin_to_single_fdf *w) {
  gsl_function *f;
  
  f = (gsl_function*)malloc (sizeof (gsl_function));
  if (f == 0) {
     GSL_ERROR_RETURN ("failed to allocate space for function struct",
		       GSL_ENOMEM, 0);
  }
  f->params = (void *)w;
  f->function = gsl_multimin_to_single_eval_fdf;

  return f;
}

void
gsl_multimin_to_single_function_fdf_free(gsl_function *f) {
  free (f);
}

