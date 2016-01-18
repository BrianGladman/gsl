/* multilarge_nlinear/nlinear.c
 * 
 * Copyright (C) 2015 Patrick Alken
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

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge_nlinear.h>

gsl_multilarge_nlinear_workspace *
gsl_multilarge_nlinear_alloc (const gsl_multilarge_nlinear_type * T,
                              const gsl_multilarge_nlinear_parameters * params,
                              const size_t n, const size_t p)
{
  gsl_multilarge_nlinear_workspace * w;

  w = (gsl_multilarge_nlinear_workspace *) calloc (1, sizeof (gsl_multilarge_nlinear_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace",
                      GSL_ENOMEM);
    }

  w->type = T;
  w->fdf = NULL;
  w->niter = 0;
  w->n = n;
  w->p = p;

  w->x = gsl_vector_alloc (p);
  if (w->x == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_NULL ("failed to allocate space for x", GSL_ENOMEM);
    }

  w->f = gsl_vector_alloc (n);
  if (w->f == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_NULL ("failed to allocate space for f", GSL_ENOMEM);
    }

  w->g = gsl_vector_alloc (p);
  if (w->g == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_NULL ("failed to allocate space for g", GSL_ENOMEM);
    }

  w->dx = gsl_vector_alloc (p);
  if (w->dx == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_NULL ("failed to allocate space for dx", GSL_ENOMEM);
    }

  w->JTJ = gsl_matrix_alloc (p, p);
  if (w->JTJ == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_NULL ("failed to allocate space for JTJ", GSL_ENOMEM);
    }

  w->state = (w->type->alloc)(params, n, p);
  if (w->state == 0)
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_NULL ("failed to allocate space for solver state",
                      GSL_ENOMEM);
    }

  return w;
}

void
gsl_multilarge_nlinear_free (gsl_multilarge_nlinear_workspace * w)
{
  if (w->state)
    (w->type->free) (w->state);

  if (w->dx)
    gsl_vector_free (w->dx);

  if (w->x)
    gsl_vector_free (w->x);

  if (w->f)
    gsl_vector_free (w->f);

  if (w->g)
    gsl_vector_free (w->g);

  if (w->JTJ)
    gsl_matrix_free (w->JTJ);

  free (w);
}

gsl_multilarge_nlinear_parameters
gsl_multilarge_nlinear_default_parameters(void)
{
  gsl_multilarge_nlinear_parameters params;

  params.scale = gsl_multilarge_nlinear_scale_more;
  params.accel = 0;
  params.avmax = 0.75;
  params.h_fvv = 0.01;

  return params;
}

const char *
gsl_multilarge_nlinear_name (const gsl_multilarge_nlinear_workspace * w)
{
  return w->type->name;
}

size_t
gsl_multilarge_nlinear_niter (const gsl_multilarge_nlinear_workspace * w)
{
  return w->niter;
}

/*
gsl_multilarge_nlinear_init()
  Initialize nonlinear least squares solver

Inputs: x   - initial estimate of parameters
        fdf - user-supplied Jacobian and residual vector callback
        w   - workspace
*/

int
gsl_multilarge_nlinear_init (const gsl_vector * x, gsl_multilarge_nlinear_fdf * fdf,
                             gsl_multilarge_nlinear_workspace * w)
{
  if (x->size != w->p)
    {
      GSL_ERROR ("x vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status;

      fdf->nevalf = 0;
      fdf->nevaldf = 0;
      fdf->nevaldff = 0;
      fdf->nevalfvv = 0;

      w->fdf = fdf;
      gsl_vector_memcpy(w->x, x);
      w->niter = 0;

      status = (w->type->init)(w->fdf, w->x, w->f,
                               w->g, w->JTJ, w->state);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

int
gsl_multilarge_nlinear_iterate (gsl_multilarge_nlinear_workspace * w)
{
  int status;

  status = (w->type->iterate) (w->fdf, w->x, w->f, w->JTJ, w->g, w->dx, w->state);
  w->niter++;

  return status;
}

gsl_vector *
gsl_multilarge_nlinear_position (const gsl_multilarge_nlinear_workspace * w)
{
  return w->x;
}

gsl_vector *
gsl_multilarge_nlinear_residual (const gsl_multilarge_nlinear_workspace * w)
{
  return w->f;
}

gsl_matrix *
gsl_multilarge_nlinear_JTJ (const gsl_multilarge_nlinear_workspace * w)
{
  return w->JTJ;
}


int
gsl_multilarge_nlinear_rcond (double * rcond, const gsl_multilarge_nlinear_workspace * w)
{
  return (w->type->rcond) (w->JTJ, rcond, w->state);
}

/*
gsl_multilarge_nlinear_driver()
  Iterate the nonlinear least squares solver until completion

Inputs: maxiter  - maximum iterations to allow
        xtol     - tolerance in step x
        gtol     - tolerance in gradient
        ftol     - tolerance in ||f||
        callback - callback function to call each iteration
        callback_params - parameters to pass to callback function
        info     - (output) info flag on why iteration terminated
                   1 = stopped due to small step size ||dx|
                   2 = stopped due to small gradient
                   3 = stopped due to small change in f
                   GSL_ETOLX = ||dx|| has converged to within machine
                               precision (and xtol is too small)
                   GSL_ETOLG = ||g||_inf is smaller than machine
                               precision (gtol is too small)
                   GSL_ETOLF = change in ||f|| is smaller than machine
                               precision (ftol is too small)
        w        - workspace

Return:
GSL_SUCCESS if converged
GSL_MAXITER if maxiter exceeded without converging
GSL_ENOPROG if no accepted step found on first iteration
*/

int
gsl_multilarge_nlinear_driver (const size_t maxiter,
                               const double xtol,
                               const double gtol,
                               const double ftol,
                               void (*callback)(const size_t iter, void *params,
                                                const gsl_multilarge_nlinear_workspace *w),
                               void *callback_params,
                               int *info,
                               gsl_multilarge_nlinear_workspace * w)
{
  int status;
  size_t iter = 0;

  /* call user callback function prior to any iterations
   * with initial system state */
  if (callback)
    callback(iter, callback_params, w);

  do
    {
      status = gsl_multilarge_nlinear_iterate (w);

      /*
       * If the solver reports no progress on the first iteration,
       * then it didn't find a single step to reduce the
       * cost function and more iterations won't help so return.
       *
       * If we get a no progress flag on subsequent iterations,
       * it means we did find a good step in a previous iteration,
       * so continue iterating since the solver has now reset
       * mu to its initial value.
       */
      if (status == GSL_ENOPROG && iter == 0)
        {
          *info = status;
          return GSL_EMAXITER;
        }

      ++iter;

      if (callback)
        callback(iter, callback_params, w);

      /* test for convergence */
      status = gsl_multilarge_nlinear_test(xtol, gtol, ftol, info, w);
    }
  while (status == GSL_CONTINUE && iter < maxiter);

  /*
   * the following error codes mean that the solution has converged
   * to within machine precision, so record the error code in info
   * and return success
   */
  if (status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG)
    {
      *info = status;
      status = GSL_SUCCESS;
    }

  /* check if max iterations reached */
  if (iter >= maxiter && status != GSL_SUCCESS)
    status = GSL_EMAXITER;

  return status;
}

/*
gsl_multilarge_nlinear_applyW()
  Apply weight matrix to (X,y) LS system

Inputs: w    - weight vector n-by-1 or NULL for W = I
        J    - Jacobian matrix n-by-p (can be NULL)
        f    - residual vector n-by-1
*/

int
gsl_multilarge_nlinear_applyW(const gsl_vector * w,
                              gsl_matrix * J, gsl_vector * f)
{
  const size_t n = w->size;

  if (n != f->size)
    {
      GSL_ERROR("f vector does not match w", GSL_EBADLEN);
    }
  else if (J != NULL && n != J->size1)
    {
      GSL_ERROR("weight vector does not match J", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      if (J == NULL)
        {
          /* construct sqrt(W) f */
          for (i = 0; i < n; ++i)
            {
              double wi = gsl_vector_get(w, i);
              double swi;
              double *fi = gsl_vector_ptr(f, i);

              if (wi < 0.0)
                wi = 0.0;

              swi = sqrt(wi);
              *fi *= swi;
            }
        }
      else
        {
          /* construct sqrt(W) J and sqrt(W) f */
          for (i = 0; i < n; ++i)
            {
              double wi = gsl_vector_get(w, i);
              double swi;
              gsl_vector_view row = gsl_matrix_row(J, i);
              double *fi = gsl_vector_ptr(f, i);

              if (wi < 0.0)
                wi = 0.0;

              swi = sqrt(wi);
              gsl_vector_scale(&row.vector, swi);
              *fi *= swi;
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_multilarge_nlinear_eval_f()
  Compute residual vector f with user callback function, and apply
weighting transform if given:

Inputs: fdf  - callback function
        x    - model parameters
        f    - (output) residual vector f(x)
*/

int
gsl_multilarge_nlinear_eval_f(gsl_multilarge_nlinear_fdf *fdf,
                              const gsl_vector *x, gsl_vector *f)
{
  int s = ((*((fdf)->f)) (x, fdf->params, f));

  ++(fdf->nevalf);

  return s;
}

/*
gsl_multilarge_nlinear_eval_df()
  Compute residual vector f with user callback function, and apply
weighting transform if given:

Inputs: fdf  - callback function
        x    - model parameters
        f    - residual vector f(x)
        JTf  - (output) J^T f
        JTJ  - (output) J^T J
*/

int
gsl_multilarge_nlinear_eval_df(gsl_multilarge_nlinear_fdf *fdf,
                               const gsl_vector *x, const gsl_vector *f,
                               gsl_vector *JTf, gsl_matrix *JTJ)
{
  int status = ((*((fdf)->df)) (x, f, fdf->params, JTf, JTJ));

  if (JTJ)
    ++(fdf->nevaldf);

  if (JTf)
    ++(fdf->nevaldff);

  return status;
}

/*
gsl_multilarge_nlinear_eval_fvv()
  Compute second directional derivative fvv with user
callback function if given

Inputs: h     - step size for finite difference, if needed
        x     - model parameters, size p
        v     - velocity vector, size p
        g     - gradient J^T f, size p
        JTJ   - J^T J, p-by-p
        fdf   - callback function
        fvv   - (output) second directional derivative vector f_vv(x), size n
        JTfvv - (output) J^T fvv
        workp - workspace, size p
*/

int
gsl_multilarge_nlinear_eval_fvv(const double h, const gsl_vector *x, const gsl_vector *v,
                                const gsl_vector *g, const gsl_matrix *JTJ,
                                gsl_multilarge_nlinear_fdf *fdf,
                                gsl_vector *fvv, gsl_vector *JTfvv, gsl_vector *workp)
{
  int status;
  
  if (fdf->fvv != NULL)
    {
      /* analytic callback available, compute f_vv(x) */
      status = ((*((fdf)->fvv)) (x, v, fdf->params, fvv));
      if (status)
        return status;

      ++(fdf->nevalfvv);

      /* compute J^T fvv */
      status = gsl_multilarge_nlinear_eval_df(fdf, x, fvv, JTfvv, NULL);
      if (status)
        return status;
    }
  else
    {
      status = gsl_multilarge_nlinear_fdJTfvv(h, x, v, g, JTJ, fdf, JTfvv, fvv, workp);
      if (status)
        return status;
    }

  return GSL_SUCCESS;
}
