/* multimin/fdfminimizer.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>

gsl_multimin_fdf_history *
gsl_multimin_fdf_history_alloc(gsl_multimin_function_fdf *fdf,
			       const gsl_vector * x)
{
  gsl_multimin_fdf_history *h;
  const size_t n = fdf->n;  
  int status;

  if (x->size != n) {
      GSL_ERROR_VAL ("vector length not compatible with function", 
                        GSL_EBADLEN, 0);
  }
  
  h = (gsl_multimin_fdf_history *) malloc(sizeof(gsl_multimin_fdf_history));

  if (h == 0) {
    GSL_ERROR_VAL ("failed to allocate space for multimin history struct",
		      GSL_ENOMEM, 0);
  }

  h->x = gsl_vector_calloc (n);
  
  if (h->x == 0) {
    free (h);
    GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
  }
  
  h->x1 = gsl_vector_calloc (n);
  
  if (h->x1 == 0) {
    free (h);
    gsl_vector_free(h->x);
    GSL_ERROR_VAL ("failed to allocate space for x1", GSL_ENOMEM, 0);
  }
  
  h->g = gsl_vector_calloc (n);
  
  if (h->g == 0) {
    free (h);
    gsl_vector_free(h->x);
    gsl_vector_free(h->x1);
    GSL_ERROR_VAL ("failed to allocate space for g", GSL_ENOMEM, 0);
  }
  
  h->g1 = gsl_vector_calloc (n);
  
  if (h->g1 == 0) {
    free (h);
    gsl_vector_free(h->x);
    gsl_vector_free(h->x1);
    gsl_vector_free(h->g);
    GSL_ERROR_VAL ("failed to allocate space for g1", GSL_ENOMEM, 0);
  }
 
  status=gsl_multimin_fdf_history_set(h,fdf,x);
  if (status != GSL_SUCCESS) {
    free (h);
    gsl_vector_free(h->x);
    gsl_vector_free(h->x1);
    gsl_vector_free(h->g);
    gsl_vector_free(h->g1);
    GSL_ERROR_VAL ("failed to set history", GSL_ENOMEM, 0);
  }
  return h;
}

int
gsl_multimin_fdf_history_set(gsl_multimin_fdf_history *h,
			     gsl_multimin_function_fdf *fdf,
			     const gsl_vector * x)
{
  gsl_vector_memcpy(h->x,x);
  GSL_MULTIMIN_FN_EVAL_F_DF(fdf,h->x,&(h->f),h->g);
  if (!finite(h->f))
    GSL_ERROR("function not continuous", GSL_EBADFUNC);
  return GSL_SUCCESS;
}

int
gsl_multimin_fdf_history_set_with_value(gsl_multimin_fdf_history *h,
					gsl_multimin_function_fdf *fdf,
					const gsl_vector * x,
					double fx)
{
  gsl_vector_memcpy(h->x,x);
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
  gsl_vector_memcpy(h->g1,h->g);
  gsl_vector_memcpy(h->x1,h->x);
  h->f1 = h->f;
  gsl_multimin_compute_evaluation_point(h->x,h->x1,step,direction);
  GSL_MULTIMIN_FN_EVAL_F_DF(fdf,h->x,&(h->f),h->g);
  if (!finite(h->f))
    GSL_ERROR("function not continuous", GSL_EBADFUNC);
  return GSL_SUCCESS;
}

int
gsl_multimin_fdf_history_step_with_value(gsl_multimin_fdf_history *h,
					 gsl_multimin_function_fdf *fdf,
					 const gsl_vector * direction,
					 double step,double fx) 
{
  gsl_vector_memcpy(h->g1,h->g);
  gsl_vector_memcpy(h->x1,h->x);
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

gsl_multimin_fdfminimizer *
gsl_multimin_fdfminimizer_alloc(const gsl_multimin_fdfminimizer_type *T,
				 gsl_multimin_function_fdf *fdf,
				 const gsl_vector * x,
				 gsl_min_bracketing_function bracket,
				 const gsl_min_fminimizer_type * T_line) 
{
  int status;
  gsl_vector *direction;
  gsl_multimin_to_single *w;
  double dummy_lower, dummy_upper;
    
  gsl_multimin_fdfminimizer *s;
  
  s = (gsl_multimin_fdfminimizer *)malloc(sizeof(gsl_multimin_fdfminimizer));
  if (s == 0) 
    {
      GSL_ERROR_VAL ("failed to allocate space for minimizer struct",
			GSL_ENOMEM, 0);
    }

  s->state = malloc(T->size);

  if (s->state == 0)
    {
      free(s);
      GSL_ERROR_VAL ("failed to allocate space for minimizer state",
			GSL_ENOMEM, 0);
    }

  status = (T->alloc)(s->state,fdf->n);

  if (status != GSL_SUCCESS)
    {
      free(s->state);
      free(s);
    
      GSL_ERROR_VAL ("failed to initialize minimizer state",
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
      GSL_ERROR_VAL ("failed to allocate history",
			GSL_ENOMEM, 0);
    }

  direction = gsl_vector_calloc (x->size);

  if (direction == 0) 
    {
      gsl_multimin_fdf_history_free(s->history);
      (T->free)(s->state);
      free(s->state);
      free(s);
      
      GSL_ERROR_VAL ("failed to allocate direction vector",
			GSL_ENOMEM, 0);
    }

  w = gsl_multimin_to_single_alloc_fdf(fdf,s->history->x,direction);

  if (w == 0) 
    {
      gsl_vector_free(direction);
      gsl_multimin_fdf_history_free(s->history);
      (T->free)(s->state);
      free(s->state);
      free(s);
      
      GSL_ERROR_VAL ("failed to allocate gsl_multimin_to_single_fdf struct",
			GSL_ENOMEM, 0);
    }
  
  s->f_directional = gsl_multimin_to_single_function_alloc(w);
  
  if (s->f_directional == 0) 
    {
      gsl_multimin_to_single_free(w);
      gsl_vector_free(direction);
      gsl_multimin_fdf_history_free(s->history);
      (T->free)(s->state);
      free(s->state);
      free(s);

      GSL_ERROR_VAL ("failed to allocate one dimensional representation",
			GSL_ENOMEM, 0);     
    }

  s->bracketing = bracket;
  s->line_search_type = T_line;
  /* Warning, we try to bypass checks */
  dummy_lower = 0;
  dummy_upper = 1;
  s->line_search = gsl_min_fminimizer_alloc(T_line);
  
  if (s->line_search == 0)
    {
      gsl_multimin_to_single_function_free(s->f_directional);
      gsl_multimin_fdf_history_free(s->history);
      (T->free)(s->state);
      free(s->state);
      free(s);

      GSL_ERROR_VAL ("failed to allocate one dimensional minimization algorithm",
			GSL_ENOMEM, 0);     
    }

  gsl_min_fminimizer_set_with_values (s->line_search, 
                                      s->f_directional, 
                                      0.5, 0.0, 
                                      dummy_lower, 1.0,
                                      dummy_upper, 1.0);

  return s;
}

void
gsl_multimin_fdfminimizer_free(gsl_multimin_fdfminimizer *s) 
{
  gsl_min_fminimizer_free(s->line_search);
  gsl_multimin_to_single_function_free(s->f_directional);
  gsl_multimin_fdf_history_free(s->history);
  (s->type->free)(s->state);
  free(s->state);
  free(s);
}

static gsl_vector * direction_internal(gsl_multimin_fdfminimizer *s);

static gsl_vector *
direction_internal(gsl_multimin_fdfminimizer *s)
{
  return ((gsl_multimin_to_single *)(s->f_directional->params))->direction;
}

const gsl_vector *
gsl_multimin_fdfminimizer_direction(gsl_multimin_fdfminimizer *s)
{
  return direction_internal(s);
}

int
gsl_multimin_fdfminimizer_next_direction(gsl_multimin_fdfminimizer *s)
{
  return s->type->direction(s->state,s->history,direction_internal(s));
}

int
gsl_multimin_fdfminimizer_bracket(gsl_multimin_fdfminimizer *s,
				   double first_step,size_t eval_max)
{
  int status;
  double bracket_lower, bracket_upper;
  double f_minimum,f_upper,f_lower;
  double minimum;

  
  bracket_upper = first_step;
  bracket_lower = 0.0;
  f_lower = s->history->f;
  f_upper = GSL_FN_EVAL(s->f_directional,first_step);
  if (!finite(f_upper))
    GSL_ERROR("function not continuous", GSL_EBADFUNC);
  status =  (s->bracketing)(s->f_directional,&minimum,
			    &f_minimum,
                            &bracket_lower, &f_lower, 
                            &bracket_upper, &f_upper,eval_max);
  if (status == GSL_SUCCESS)
    {
      return gsl_min_fminimizer_set_with_values(s->line_search,
						s->f_directional,
						minimum,f_minimum,
						bracket_lower, f_lower,
                                                bracket_upper, f_upper);
    }
  else
    {
      return status;
    }
}

int
gsl_multimin_fdfminimizer_iterate(gsl_multimin_fdfminimizer *s)
{
  return gsl_min_fminimizer_iterate(s->line_search);
}

int
gsl_multimin_fdfminimizer_step_with_value(gsl_multimin_fdfminimizer *s,
					   double step,double f_at_end)
{
  return 
    gsl_multimin_fdf_history_step_with_value(s->history,s->fdf,
					     gsl_multimin_fdfminimizer_direction(s),
					     step,
					     f_at_end);
}

int
gsl_multimin_fdfminimizer_step(gsl_multimin_fdfminimizer *s,
				double step) 
{
  return 
    gsl_multimin_fdf_history_step(s->history,s->fdf,
				  gsl_multimin_fdfminimizer_direction(s),
				  step);  
}

int
gsl_multimin_fdfminimizer_best_step(gsl_multimin_fdfminimizer *s)
{
  return 
    gsl_multimin_fdfminimizer_step_with_value(s,s->line_search->minimum,
					       s->line_search->f_minimum);
}

int
gsl_multimin_fdfminimizer_restart(gsl_multimin_fdfminimizer *s)
{
  return (s->type->restart)(s->state);
}
