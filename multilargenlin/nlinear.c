/* multilargenlin/nlinear.c
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
#include <gsl/gsl_multilarge_nlin.h>

#include "oct.c"

gsl_multilarge_nlinear_workspace *
gsl_multilarge_nlinear_alloc (const gsl_multilarge_nlinear_type * T, 
                              const size_t p)
{
  gsl_multilarge_nlinear_workspace * w;

  w = (gsl_multilarge_nlinear_workspace *) calloc (1, sizeof (gsl_multilarge_nlinear_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace",
                      GSL_ENOMEM);
    }

  w->type = T;
  w->callback = NULL;
  w->niter = 0;
  w->p = p;

  w->x = gsl_vector_alloc (p);
  if (w->x == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_NULL ("failed to allocate space for x", GSL_ENOMEM);
    }

  w->dx = gsl_vector_alloc (p);
  if (w->dx == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_NULL ("failed to allocate space for dx", GSL_ENOMEM);
    }

  w->state = (w->type->alloc)(p);
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

  free (w);
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
gsl_multilarge_nlinear_init (const gsl_vector * x, gsl_multilarge_function_fdf * fdf,
                             gsl_multilarge_nlinear_workspace * w)
{
  if (x->size != w->p)
    {
      GSL_ERROR ("x vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status;

      w->callback = fdf;
      gsl_vector_memcpy(w->x, x);
      w->niter = 0;

      status = (w->type->init)(x, fdf, w, w->state);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

/*
gsl_multilarge_nlinear_accumulate()
  Accumulate Jacobian and residual vector into linear LS system

Inputs: J - Jacobian matrix
        f - residual vector
        w - workspace
*/

int
gsl_multilarge_nlinear_accumulate (gsl_matrix * J, gsl_vector * f,
                                   gsl_multilarge_nlinear_workspace * w)
{
  const size_t n = J->size1;
  const size_t p = J->size2;

  if (f->size != n)
    {
      GSL_ERROR ("f vector does not match J", GSL_EBADLEN);
    }
  else if (p != w->p)
    {
      GSL_ERROR ("Jacobian has wrong number of columns", GSL_EBADLEN);
    }
  else
    {
      int status;

      status = (w->type->accum) (J, f, w->state);

      return status;
    }
}

int
gsl_multilarge_nlinear_iterate (gsl_multilarge_nlinear_workspace * w)
{
  int status;

  status = (w->type->iterate) (w->x, w->dx, w->callback, w, w->state);

  return status;
}

double
gsl_multilarge_nlinear_normf (const gsl_multilarge_nlinear_workspace * w)
{
  return (w->type->normf) (w->state);
}

gsl_vector *
gsl_multilarge_nlinear_position (const gsl_multilarge_nlinear_workspace * w)
{
  return w->x;
}

/*
gsl_multilarge_nlinear_driver()
  Iterate the nonlinear least squares solver until completion

Inputs: maxiter - maximum iterations to allow
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
        w       - workspace

Return: GSL_SUCCESS if converged, GSL_MAXITER if maxiter exceeded without
converging
*/

int
gsl_multilarge_nlinear_driver (const size_t maxiter,
                               const double xtol,
                               const double gtol,
                               const double ftol,
                               int *info,
                               gsl_multilarge_nlinear_workspace *w)
{
  int status;
  size_t iter = 0;

  do
    {
      status = gsl_multilarge_nlinear_iterate (w);

      /*
       * if status is GSL_ENOPROG or GSL_SUCCESS, continue iterating,
       * otherwise the method has converged with a GSL_ETOLx flag
       */
      if (status != GSL_SUCCESS && status != GSL_ENOPROG)
        break;

      /* test for convergence */
      status = gsl_multilarge_nlinear_test(xtol, gtol, ftol, info, w);
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
}
