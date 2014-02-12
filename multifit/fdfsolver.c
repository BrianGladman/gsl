/* multifit/fdfsolver.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit_nlin.h>

gsl_multifit_fdfsolver *
gsl_multifit_fdfsolver_alloc (const gsl_multifit_fdfsolver_type * T, 
                              size_t n, size_t p)
{
  int status;

  gsl_multifit_fdfsolver * s;

  if (n < p)
    {
      GSL_ERROR_VAL ("insufficient data points, n < p", GSL_EINVAL, 0);
    }

  s = (gsl_multifit_fdfsolver *) calloc (1, sizeof (gsl_multifit_fdfsolver));
  if (s == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for multifit solver struct",
                     GSL_ENOMEM, 0);
    }

  s->x = gsl_vector_calloc (p);

  if (s->x == 0) 
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
    }

  s->f = gsl_vector_calloc (n);

  if (s->f == 0) 
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for f", GSL_ENOMEM, 0);
    }

  s->J = gsl_matrix_calloc (n,p);

  if (s->J == 0) 
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for g", GSL_ENOMEM, 0);
    }

  s->dx = gsl_vector_calloc (p);

  if (s->dx == 0) 
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for dx", GSL_ENOMEM, 0);
    }

  s->g = gsl_vector_alloc (p);

  if (s->g == 0) 
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for g", GSL_ENOMEM, 0);
    }

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for multifit solver state",
                     GSL_ENOMEM, 0);
    }

  s->type = T ;

  status = (s->type->alloc)(s->state, n, p);

  if (status != GSL_SUCCESS)
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to set solver", status, 0);
    }

  s->fdf = NULL;
  
  s->niter = 0;

  return s;
}

int
gsl_multifit_fdfsolver_set (gsl_multifit_fdfsolver * s, 
                            gsl_multifit_function_fdf * f, 
                            const gsl_vector * x)
{
  if (s->f->size != f->n)
    {
      GSL_ERROR ("function size does not match solver", GSL_EBADLEN);
    }

  if (s->x->size != x->size)
    {
      GSL_ERROR ("vector length does not match solver", GSL_EBADLEN);
    }  

  s->fdf = f;
  gsl_vector_memcpy(s->x,x);
  s->niter = 0;
  
  return (s->type->set) (s->state, s->fdf, s->x, s->f, s->J, s->dx);
}

int
gsl_multifit_fdfsolver_iterate (gsl_multifit_fdfsolver * s)
{
  int status =
    (s->type->iterate) (s->state, s->fdf, s->x, s->f, s->J, s->dx);

  s->niter++;

  return status;
}

/*
gsl_multifit_fdfsolver_driver()
  Iterate the nonlinear least squares solver until completion

Inputs: s - fdfsolver
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
gsl_multifit_fdfsolver_driver (gsl_multifit_fdfsolver * s,
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
      status = gsl_multifit_fdfsolver_iterate (s);

      /*
       * if status is GSL_ENOPROG or GSL_SUCCESS, continue iterating,
       * otherwise the method has converged with a GSL_ETOLx flag
       */
      if (status != GSL_SUCCESS && status != GSL_ENOPROG)
        break;

      /* test for convergence */
      status = gsl_multifit_test_convergence(s, xtol, gtol, ftol, info);
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
} /* gsl_multifit_fdfsolver_driver() */

void
gsl_multifit_fdfsolver_free (gsl_multifit_fdfsolver * s)
{
  RETURN_IF_NULL (s);

  if (s->state)
    {
      (s->type->free) (s->state);
      free (s->state);
    }

  if (s->dx)
    gsl_vector_free (s->dx);

  if (s->x)
    gsl_vector_free (s->x);

  if (s->f)
    gsl_vector_free (s->f);

  if (s->J)
    gsl_matrix_free (s->J);

  if (s->g)
    gsl_vector_free (s->g);

  free (s);
}

const char *
gsl_multifit_fdfsolver_name (const gsl_multifit_fdfsolver * s)
{
  return s->type->name;
}

gsl_vector *
gsl_multifit_fdfsolver_position (const gsl_multifit_fdfsolver * s)
{
  return s->x;
}

