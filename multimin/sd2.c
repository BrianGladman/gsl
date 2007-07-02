/* multimin/steepest_descent.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* steepest_descent.c -- the steepest descent algorithm */

#include <config.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>

typedef struct
{
  double step;
  double max_step;
  gsl_vector *x1;
  gsl_vector *g1;
}
steepest_descent_state_t;

static int
steepest_descent_alloc (void *vstate, size_t n)
{
  steepest_descent_state_t *state = (steepest_descent_state_t *) vstate;

  state->x1 = gsl_vector_alloc (n);

  if (state->x1 == NULL)
    {
      GSL_ERROR ("failed to allocate space for x1", GSL_ENOMEM);
    }

  state->g1 = gsl_vector_alloc (n);

  if (state->g1 == NULL)
    {
      gsl_vector_free (state->x1);
      GSL_ERROR ("failed to allocate space for g1", GSL_ENOMEM);
    }

  return GSL_SUCCESS;
}

static int
steepest_descent_set (void *vstate, gsl_multimin_function_fdf * fdf,
                      const gsl_vector * x, double *f,
                      gsl_vector * gradient, double step_size)
{
  steepest_descent_state_t *state = (steepest_descent_state_t *) vstate;

  GSL_MULTIMIN_FN_EVAL_F_DF (fdf, x, f, gradient);

  state->step = step_size;
  state->max_step = step_size;

  return GSL_SUCCESS;
}


static void
steepest_descent_free (void *vstate)
{
  steepest_descent_state_t *state = (steepest_descent_state_t *) vstate;

  gsl_vector_free (state->x1);
  gsl_vector_free (state->g1);
}

static int
steepest_descent_restart (void *vstate)
{
  steepest_descent_state_t *state = (steepest_descent_state_t *) vstate;

  state->step = state->max_step;

  return GSL_SUCCESS;
}

static int
steepest_descent_iterate (void *vstate, gsl_multimin_function_fdf * fdf,
                          gsl_vector * x, double *f,
                          gsl_vector * gradient, gsl_vector * dx)
{
  steepest_descent_state_t *state = (steepest_descent_state_t *) vstate;

  gsl_vector *x1 = state->x1;
  gsl_vector *g1 = state->g1;

  double f0 = *f;
  double f1, fm, f_trial;
  double step0 = 0.0, step1, stepm, step_trial, step = state->step;
  size_t iter = 0;

  /* compute new trial point at x1= x - step * dir, where dir is the
     normalized gradient */

  double gnorm = gsl_blas_dnrm2 (gradient);

  gsl_vector_set_zero (dx);
  gsl_blas_daxpy (-step / gnorm, gradient, dx);

  gsl_vector_memcpy (x1, x);
  gsl_blas_daxpy (1.0, dx, x1);

  /* evaluate function and gradient at new point x1 */

  GSL_MULTIMIN_FN_EVAL_F_DF (fdf, x1, &f1, g1);

  if (f1 < f0)
    {
      state->step = 2.0 * step;

      gsl_vector_memcpy (x, x1);
      gsl_vector_memcpy (gradient, g1);

      *f = f1;

      return GSL_SUCCESS;
    }

  step1 = step;

trial:

  stepm = 0.5 * step1;

  gsl_vector_set_zero (dx);
  gsl_blas_daxpy (-stepm / gnorm, gradient, dx);
  
  gsl_vector_memcpy (x1, x);
  gsl_blas_daxpy (1.0, dx, x1);

  fm = GSL_MULTIMIN_FN_EVAL_F (fdf, x1);
  printf("trying stepm = %.18e fm=%g\n", stepm, fm);
  if (fm >= f0)
    {
      /* downhill step failed, reduce step-size and try again */
      f1 = fm;
      step1 = stepm;
      goto trial;
    }

  /* We now have a triplet (0,f0) (stepm, fm) (step1, f1) */

minimize:
  iter++;

  if ((stepm - step0) > (step1 - stepm))
    {
      step_trial = stepm - 0.38 * (stepm - step0);
    }
  else
    {
      step_trial = stepm + 0.38 * (step1 - stepm);
    }
  
  gsl_vector_set_zero (dx);
  gsl_blas_daxpy (-step_trial / gnorm, gradient, dx);
  
  gsl_vector_memcpy (x1, x);
  gsl_blas_daxpy (1.0, dx, x1);
  
  f_trial = GSL_MULTIMIN_FN_EVAL_F (fdf, x1);
  
  if (f_trial > fm)
    {
      if (step_trial < stepm)
        {
          step0 = step_trial;
          f0 = f_trial;
        }
      else
        {
          step1 = step_trial;
          f1 = f_trial;
        }
    }
  else
    {
      if (step_trial < stepm)
        {
          step1 = stepm;
          f1 = fm;
        }
      else
        {
          step0 = stepm;
          f0 = fm;
        }

      stepm = step_trial;
      fm = f_trial;
    }

  printf("f = %.18e\n", fm);
  
  if (iter > 100) 
    {
      gsl_vector_set_zero (dx);
      gsl_blas_daxpy (-stepm / gnorm, gradient, dx);
      gsl_blas_daxpy (1.0, dx, x);
      
      GSL_MULTIMIN_FN_EVAL_DF (fdf, x, gradient);
      *f = fm;
      state->step = stepm;

      return GSL_SUCCESS;
    }

  goto minimize;
}

static const gsl_multimin_fdfminimizer_type steepest_descent_type =
  { "steepest_descent",         /* name */
  sizeof (steepest_descent_state_t),
  &steepest_descent_alloc,
  &steepest_descent_set,
  &steepest_descent_iterate,
  &steepest_descent_restart,
  &steepest_descent_free
};

const gsl_multimin_fdfminimizer_type
  *gsl_multimin_fdfminimizer_steepest_descent = &steepest_descent_type;
