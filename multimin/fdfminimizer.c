#include <gsl_errno.h>
#include <gsl_multimin.h>

gsl_multimin_fdf_history *
gsl_multimin_fdf_history_alloc(gsl_multimin_function_fdf *fdf,
			       gsl_vector * x)
{
  gsl_multimin_fdf_history *h;
  const size_t n = fdf->n;  
  int status;

  if (x->size != n) {
      GSL_ERROR_RETURN ("vector length not compatible with function", 
                        GSL_EBADLEN, 0);
  }
  
  h = (gsl_multimin_fdf_history *) malloc(sizeof(gsl_multimin_fdf_history));

  if (h == 0) {
    GSL_ERROR_RETURN ("failed to allocate space for multimin history struct",
		      GSL_ENOMEM, 0);
  }

  h->x = gsl_vector_calloc (n);
  
  if (h->x == 0) {
    free (h);
    GSL_ERROR_RETURN ("failed to allocate space for x", GSL_ENOMEM, 0);
  }
  
  h->x1 = gsl_vector_calloc (n);
  
  if (h->x1 == 0) {
    free (h);
    gsl_vector_free(h->x);
    GSL_ERROR_RETURN ("failed to allocate space for x1", GSL_ENOMEM, 0);
  }
  
  h->g = gsl_vector_calloc (n);
  
  if (h->g == 0) {
    free (h);
    gsl_vector_free(h->x);
    gsl_vector_free(h->x1);
    GSL_ERROR_RETURN ("failed to allocate space for g", GSL_ENOMEM, 0);
  }
  
  h->g1 = gsl_vector_calloc (n);
  
  if (h->g1 == 0) {
    free (h);
    gsl_vector_free(h->x);
    gsl_vector_free(h->x1);
    gsl_vector_free(h->g);
    GSL_ERROR_RETURN ("failed to allocate space for g1", GSL_ENOMEM, 0);
  }
 
  status=gsl_multimin_fdf_history_set(h,fdf,x);
  if (status != GSL_SUCCESS) {
    free (h);
    gsl_vector_free(h->x);
    gsl_vector_free(h->x1);
    gsl_vector_free(h->g);
    gsl_vector_free(h->g1);
    GSL_ERROR_RETURN ("failed to set history", GSL_ENOMEM, 0);
  }
  return h;
}

int
gsl_multimin_fdf_history_set(gsl_multimin_fdf_history *h,
			     gsl_multimin_function_fdf *fdf,
			     gsl_vector * x)
{
  gsl_vector_copy(h->x,x);
  GSL_MULTIMIN_FN_EVAL_F_DF(fdf,h->x,&(h->f),h->g);
  if (!GSL_IS_REAL(h->f))
    GSL_ERROR("function not continuous", GSL_EBADFUNC);
  return GSL_SUCCESS;
}

int
gsl_multimin_fdf_history_set_with_value(gsl_multimin_fdf_history *h,
					gsl_multimin_function_fdf *fdf,
					gsl_vector * x,
					double fx)
{
  gsl_vector_copy(h->x,x);
  GSL_MULTIMIN_FN_EVAL_DF(fdf,h->x,h->g);
  h->f = fx;
  return GSL_SUCCESS;
}

int
gsl_multimin_fdf_history_step(gsl_multimin_fdf_history *h,
			      gsl_multimin_function_fdf *fdf,
			      const gsl_vector * direction,
			      double step)
{
  gsl_vector_copy(h->g1,h->g);
  gsl_vector_copy(h->x1,h->x);
  h->f1 = h->f;
  gsl_multimin_compute_evaluation_point(h->x,h->x1,step,direction);
  GSL_MULTIMIN_FN_EVAL_F_DF(fdf,h->x,&(h->f),h->g);
  if (!GSL_IS_REAL(h->f))
    GSL_ERROR("function not continuous", GSL_EBADFUNC);
  return GSL_SUCCESS;
}

int
gsl_multimin_fdf_history_step_with_value(gsl_multimin_fdf_history *h,
					 gsl_multimin_function_fdf *fdf,
					 const gsl_vector * direction,
					 double step,double fx) 
{
  gsl_vector_copy(h->g1,h->g);
  gsl_vector_copy(h->x1,h->x);
  h->f1 = h->f;
  gsl_multimin_compute_evaluation_point(h->x,h->x1,step,direction);
  GSL_MULTIMIN_FN_EVAL_DF(fdf,h->x,h->g);
  h->f = fx;
  return GSL_SUCCESS;
}

void
gsl_multimin_fdf_history_free(gsl_multimin_fdf_history *h)
{
  gsl_vector_free(h->x);
  gsl_vector_free(h->x1);
  gsl_vector_free(h->g);
  gsl_vector_free(h->g1);
  free(h);
}

gsl_multimin_fdf_minimizer *
gsl_multimin_fdf_minimizer_alloc(const gsl_multimin_fdf_minimizer_type *T,
				 gsl_multimin_function_fdf *fdf,
				 gsl_vector * x) 
{
  int status;
  gsl_vector *direction;
  gsl_multimin_to_single_fdf *w;
    
  gsl_multimin_fdf_minimizer *s;
  
  s = (gsl_multimin_fdf_minimizer *)malloc(sizeof(gsl_multimin_fdf_minimizer));
  if (s == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for minimizer struct",
			GSL_ENOMEM, 0);
    }

  s->state = malloc(T->size);

  if (s->state == 0)
    {
      free(s);
      GSL_ERROR_RETURN ("failed to allocate space for minimizer state",
			GSL_ENOMEM, 0);
    }

  status = (T->alloc)(s->state,fdf->n);

  if (status != GSL_SUCCESS)
    {
      free(s->state);
      free(s);
    
      GSL_ERROR_RETURN ("failed to initialize minimizer state",
			GSL_ENOMEM, 0);
    }
  
  s->type = T;
  s->fdf = fdf;
  s->history = gsl_multimin_fdf_history_alloc(fdf,x);

  if (s->history == 0)
    {
      (T->free)(s->state);
      free(s->state);
      free(s);
      GSL_ERROR_RETURN ("failed to allocate history",
			GSL_ENOMEM, 0);
    }

  direction = gsl_vector_calloc (x->size);

  if (direction == 0) 
    {
      gsl_multimin_fdf_history_free(s->history);
      (T->free)(s->state);
      free(s->state);
      free(s);
      
      GSL_ERROR_RETURN ("failed to allocate direction vector",
			GSL_ENOMEM, 0);
    }

  w = gsl_multimin_to_single_fdf_alloc(fdf,s->history->x,direction);

  if (w == 0) 
    {
      gsl_vector_free(direction);
      gsl_multimin_fdf_history_free(s->history);
      (T->free)(s->state);
      free(s->state);
      free(s);
      
      GSL_ERROR_RETURN ("failed to allocate gsl_multimin_to_single_fdf struct",
			GSL_ENOMEM, 0);
    }
  
  s->f_directional = gsl_multimin_to_single_function_fdf_alloc(w);
  
  if (s->f_directional == 0) 
    {
      gsl_multimin_to_single_fdf_free(w);
      gsl_vector_free(direction);
      gsl_multimin_fdf_history_free(s->history);
      (T->free)(s->state);
      free(s->state);
      free(s);

      GSL_ERROR_RETURN ("failed to allocate one dimensional representation",
			GSL_ENOMEM, 0);     
    }

  return s;
}

void
gsl_multimin_fdf_minimizer_free(gsl_multimin_fdf_minimizer *s) 
{
  gsl_multimin_fdf_history_free(s->history);
  (s->type->free)(s->state);
  free(s->state);
  free(s);
}

gsl_vector *
gsl_multimin_fdf_minimizer_direction_internal(gsl_multimin_fdf_minimizer *s)
{
  return ((gsl_multimin_to_single_fdf *)(s->f_directional->params))->params->direction;
}

const gsl_vector *
gsl_multimin_fdf_minimizer_direction(gsl_multimin_fdf_minimizer *s)
{
  return gsl_multimin_fdf_minimizer_direction_internal(s);
}

int
gsl_multimin_fdf_minimizer_next_direction(gsl_multimin_fdf_minimizer *s)
{
  return s->type->direction(s->state,s->history,gsl_multimin_fdf_minimizer_direction_internal(s));
}

int
gsl_multimin_fdf_minimizer_step_with_value(gsl_multimin_fdf_minimizer *s,
					   double step,double f_at_end)
{
  return 
    gsl_multimin_fdf_history_step_with_value(s->history,s->fdf,
					     gsl_multimin_fdf_minimizer_direction(s),
					     step,
					     f_at_end);
}

int
gsl_multimin_fdf_minimizer_step(gsl_multimin_fdf_minimizer *s,
				double step) 
{
  return 
    gsl_multimin_fdf_history_step(s->history,s->fdf,
				  gsl_multimin_fdf_minimizer_direction(s),
				  step);  
}

int
gsl_multimin_fdf_minimizer_restart(gsl_multimin_fdf_minimizer *s)
{
  return (s->type->restart)(s->state);
}
