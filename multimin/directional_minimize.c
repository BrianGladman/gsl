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
gsl_multimin_compute_ep(gsl_multimin_to_single *params,double x) 
{
  gsl_multimin_compute_evaluation_point(params->evaluation_point,
					params->starting_point,
					x,
					params->direction);
}

double gsl_multimin_to_single_eval(double x, void * oparams)
{
  gsl_multimin_to_single *params = (gsl_multimin_to_single *)oparams;
  
  gsl_multimin_compute_ep(params,x);
  return GSL_MULTIMIN_FN_EVAL(&(params->f),params->evaluation_point);
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

  wrapper->evaluation_point = gsl_vector_calloc (starting_point->size); 

  if (wrapper->evaluation_point == 0) {
    free(wrapper);

    GSL_ERROR_RETURN ("failed to allocate space for evaluation_point", GSL_ENOMEM, 0);
  }

  status = gsl_multimin_to_single_set(wrapper,f,starting_point,direction);

  if (status != GSL_SUCCESS) {
    free(wrapper);
    gsl_vector_free (wrapper->evaluation_point);
    
    GSL_ERROR_RETURN ("failed to set multimin wrapper params", status, 0);
  }

  return wrapper;
}

int
gsl_multimin_to_single_set(gsl_multimin_to_single *w,
			   const gsl_multimin_function *f,
			   gsl_vector * starting_point,
			   gsl_vector * direction) {
  if (w->evaluation_point->size != starting_point->size 
      || w->evaluation_point->size != direction->size) {
      
    GSL_ERROR_RETURN ("vector length not compatible with function", 
		      GSL_EBADLEN, 0);
  }
  w->starting_point = starting_point;
  w->direction = direction;
  w->f.f = f->f;
  w->f.n = f->n;
  w->f.params = f->params;
  return GSL_SUCCESS;
}

void
gsl_multimin_to_single_free(gsl_multimin_to_single *w) {
  gsl_vector_free (w->evaluation_point);
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
gsl_multimin_to_single_function_free(gsl_function *f) 
{  
  gsl_multimin_to_single *params = (gsl_multimin_to_single *)(f->params);
  gsl_multimin_to_single_free(params);
  free (f);
}

gsl_multimin_to_single *
gsl_multimin_to_single_alloc_fdf(const gsl_multimin_function_fdf *fdf,
				 gsl_vector * starting_point,
				 gsl_vector * direction) {
  gsl_multimin_function f;

  f.f = fdf->f;
  f.n = fdf->n;
  f.params = fdf->params;

  return gsl_multimin_to_single_alloc(&f,starting_point,direction);
}

int
gsl_multimin_to_single_set_fdf(gsl_multimin_to_single *w,
			       const gsl_multimin_function_fdf *fdf,
			       gsl_vector * starting_point,
			       gsl_vector * direction) {
  gsl_multimin_function f;

  f.f = fdf->f;
  f.n = fdf->n;
  f.params = fdf->params;

  return gsl_multimin_to_single_set(w,&f,starting_point,direction);
}


