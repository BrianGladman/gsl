/* multifit_nlinear/fdf.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
#include <gsl/gsl_multifit_nlinear.h>

gsl_multifit_nlinear_workspace *
gsl_multifit_nlinear_alloc (const gsl_multifit_nlinear_type * T, 
                            const size_t n, const size_t p)
{
  int status;

  gsl_multifit_nlinear_workspace * w;

  if (n < p)
    {
      GSL_ERROR_VAL ("insufficient data points, n < p", GSL_EINVAL, 0);
    }

  w = calloc (1, sizeof (gsl_multifit_nlinear_workspace));
  if (w == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for multifit workspace",
                     GSL_ENOMEM, 0);
    }

  w->x = gsl_vector_calloc (p);

  if (w->x == 0) 
    {
      gsl_multifit_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
    }

  w->f = gsl_vector_calloc (n);

  if (w->f == 0) 
    {
      gsl_multifit_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for f", GSL_ENOMEM, 0);
    }

  w->dx = gsl_vector_calloc (p);

  if (w->dx == 0) 
    {
      gsl_multifit_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for dx", GSL_ENOMEM, 0);
    }

  w->g = gsl_vector_alloc (p);

  if (w->g == 0) 
    {
      gsl_multifit_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for g", GSL_ENOMEM, 0);
    }

  w->sqrt_wts = gsl_vector_calloc (n);

  if (w->sqrt_wts == 0) 
    {
      gsl_multifit_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for sqrt_wts", GSL_ENOMEM, 0);
    }

  w->state = calloc (1, T->size);

  if (w->state == 0)
    {
      gsl_multifit_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for multifit solver state",
                     GSL_ENOMEM, 0);
    }

  w->type = T ;

  status = (w->type->alloc)(w->state, n, p);

  if (status != GSL_SUCCESS)
    {
      gsl_multifit_nlinear_free (w);
      GSL_ERROR_VAL ("failed to set solver", status, 0);
    }

  w->fdf = NULL;
  
  w->niter = 0;

  return w;
}

int
gsl_multifit_nlinear_set (gsl_multifit_nlinear_workspace * w, 
                          gsl_multifit_nlinear_fdf * f,
                          const gsl_vector * x)
{
  return gsl_multifit_nlinear_wset(w, f, x, NULL);
}

int
gsl_multifit_nlinear_wset (gsl_multifit_nlinear_workspace * w,
                           gsl_multifit_nlinear_fdf * f, 
                           const gsl_vector * x,
                           const gsl_vector * wts)
{
  const size_t n = w->f->size;

  if (n != f->n)
    {
      GSL_ERROR ("function size does not match workspace", GSL_EBADLEN);
    }
  else if (w->x->size != x->size)
    {
      GSL_ERROR ("vector length does not match workspace", GSL_EBADLEN);
    }
  else if (wts != NULL && n != wts->size)
    {
      GSL_ERROR ("weight vector length does not match workspace", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      w->fdf = f;
      gsl_vector_memcpy(w->x, x);
      w->niter = 0;

      if (wts)
        {
          for (i = 0; i < n; ++i)
            {
              double wi = gsl_vector_get(wts, i);
              gsl_vector_set(w->sqrt_wts, i, sqrt(wi));
            }
        }
      else
        gsl_vector_set_all(w->sqrt_wts, 1.0);
  
      return (w->type->set) (w->state, w->sqrt_wts, w->fdf, w->x, w->f, w->dx);
    }
}

int
gsl_multifit_nlinear_iterate (gsl_multifit_nlinear_workspace * w)
{
  int status =
    (w->type->iterate) (w->state, w->sqrt_wts, w->fdf, w->x, w->f, w->dx);

  w->niter++;

  return status;
}

/*
gsl_multifit_nlinear_driver()
  Iterate the nonlinear least squares solver until completion

Inputs: w       - workspace
        maxiter - maximum iterations to allow
        xtol    - tolerance in step x
        gtol    - tolerance in gradient
        ftol    - tolerance in ||f||
        info    - (output) info flag on why iteration terminated
                  1 = stopped due to small step size ||dx|
                  2 = stopped due to small gradient
                  3 = stopped due to small change in f
                  GSL_ETOLX = ||dx|| has converged to within machine
                              precision (and xtol is too small)
                  GSL_ETOLG = ||g||_inf is smaller than machine
                              precision (gtol is too small)
                  GSL_ETOLF = change in ||f|| is smaller than machine
                              precision (ftol is too small)

Return: GSL_SUCCESS if converged, GSL_MAXITER if maxiter exceeded without
converging
*/

int
gsl_multifit_nlinear_driver (gsl_multifit_nlinear_workspace * w,
                             const size_t maxiter,
                             const double xtol,
                             const double gtol,
                             const double ftol,
                             int *info)
{
  int status;
  size_t iter = 0;

  do
    {
      status = gsl_multifit_nlinear_iterate (w);

      /*
       * if status is GSL_ENOPROG or GSL_SUCCESS, continue iterating,
       * otherwise the method has converged with a GSL_ETOLx flag
       */
      if (status != GSL_SUCCESS && status != GSL_ENOPROG)
        break;

      /* test for convergence */
      status = gsl_multifit_nlinear_test(w, xtol, gtol, ftol, info);
    }
  while (status == GSL_CONTINUE && ++iter < maxiter);

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
} /* gsl_multifit_nlinear_driver() */

int
gsl_multifit_nlinear_jac (gsl_multifit_nlinear_workspace * w, gsl_matrix * J)
{
  const size_t n = w->f->size;
  const size_t p = w->x->size;

  if (n != J->size1 || p != J->size2)
    {
      GSL_ERROR ("Jacobian dimensions do not match workspace", GSL_EBADLEN);
    }
  else
    {
      return (w->type->jac) (w->state, J);
    }
} /* gsl_multifit_nlinear_jac() */

void
gsl_multifit_nlinear_free (gsl_multifit_nlinear_workspace * w)
{
  RETURN_IF_NULL (w);

  if (w->state)
    {
      (w->type->free) (w->state);
      free (w->state);
    }

  if (w->dx)
    gsl_vector_free (w->dx);

  if (w->x)
    gsl_vector_free (w->x);

  if (w->f)
    gsl_vector_free (w->f);

  if (w->sqrt_wts)
    gsl_vector_free (w->sqrt_wts);

  if (w->g)
    gsl_vector_free (w->g);

  free (w);
}

const char *
gsl_multifit_nlinear_name (const gsl_multifit_nlinear_workspace * w)
{
  return w->type->name;
}

gsl_vector *
gsl_multifit_nlinear_position (const gsl_multifit_nlinear_workspace * w)
{
  return w->x;
}

gsl_vector *
gsl_multifit_nlinear_residual (const gsl_multifit_nlinear_workspace * w)
{
  return w->f;
}

size_t
gsl_multifit_nlinear_niter (const gsl_multifit_nlinear_workspace * w)
{
  return w->niter;
}

/*
gsl_multifit_eval_wf()
  Compute residual vector y with user callback function, and apply
weighting transform if given:

y~ = sqrt(W) y

Inputs: fdf  - callback function
        x    - model parameters
        swts - weight matrix sqrt(W) = sqrt(diag(w1,w2,...,wn))
               set to NULL for unweighted fit
        y    - (output) (weighted) residual vector
               y_i = sqrt(w_i) f_i where f_i is unweighted residual
*/

int
gsl_multifit_eval_wf(gsl_multifit_nlinear_fdf *fdf, const gsl_vector *x,
                     const gsl_vector *swts, gsl_vector *y)
{
  int s = ((*((fdf)->f)) (x, fdf->params, y));
  ++(fdf->nevalf);

  /* y <- sqrt(W) y */
  if (swts)
    gsl_vector_mul(y, swts);

  return s;
}

/*
gsl_multifit_eval_wdf()
  Compute Jacobian matrix J with user callback function, and apply
weighting transform if given:

J~ = sqrt(W) J

Inputs: fdf  - callback function
        x    - model parameters
        swts - weight matrix W = diag(w1,w2,...,wn)
               set to NULL for unweighted fit
        dy   - (output) (weighted) Jacobian matrix
               dy = sqrt(W) dy where dy is unweighted Jacobian
*/

int
gsl_multifit_eval_wdf(gsl_multifit_nlinear_fdf *fdf, const gsl_vector *x,
                      const gsl_vector *swts, gsl_matrix *dy)
{
  int status = ((*((fdf)->df)) (x, fdf->params, dy));

  ++(fdf->nevaldf);

  /* J <- sqrt(W) J */
  if (swts)
    {
      const size_t n = swts->size;
      size_t i;

      for (i = 0; i < n; ++i)
        {
          double swi = gsl_vector_get(swts, i);
          gsl_vector_view v = gsl_matrix_row(dy, i);

          gsl_vector_scale(&v.vector, swi);
        }
    }

  return status;
}
