/* multilarge.c
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>

gsl_multilarge_linear_workspace *
gsl_multilarge_linear_alloc(const gsl_multilarge_linear_type *T,
                            const size_t nmax, const size_t p)
{
  gsl_multilarge_linear_workspace *w;

  w = calloc(1, sizeof(gsl_multilarge_linear_workspace));
  if (w == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for workspace",
                     GSL_ENOMEM);
    }

  w->type = T;

  w->state = w->type->alloc(nmax, p);
  if (w->state == NULL)
    {
      gsl_multilarge_linear_free(w);
      GSL_ERROR_NULL("failed to allocate space for multilarge state",
                     GSL_ENOMEM);
    }

  w->nmax = nmax;
  w->p = p;

  /* initialize newly allocated state */
  gsl_multilarge_linear_reset(w);

  return w;
}

void
gsl_multilarge_linear_free(gsl_multilarge_linear_workspace *w)
{
  RETURN_IF_NULL(w);

  if (w->state)
    w->type->free(w->state);

  free(w);
}

const char *
gsl_multilarge_linear_name(const gsl_multilarge_linear_workspace *w)
{
  return w->type->name;
}

int
gsl_multilarge_linear_reset(gsl_multilarge_linear_workspace *w)
{
  int status = w->type->reset(w->state);
  return status;
}

int
gsl_multilarge_linear_accumulate(const gsl_matrix * X, const gsl_vector * y,
                                 gsl_multilarge_linear_workspace * w)
{
  int status = w->type->accumulate(X, y, w->state);
  return status;
}

int
gsl_multilarge_linear_solve(const double lambda, gsl_vector * c,
                            double * rnorm, double * snorm,
                            gsl_multilarge_linear_workspace * w)
{
  int status = w->type->solve(lambda, c, rnorm, snorm, w->state);
  return status;
}

/*
gsl_multilarge_linear_wstdform1()
  Using regularization matrix
L = diag(l_1,l_2,...,l_p), transform to Tikhonov standard form:

X~ = sqrt(W) X L^{-1}
y~ = sqrt(W) y
c~ = L c

Inputs: L    - Tikhonov matrix as a vector of diagonal elements p-by-1;
               or NULL for L = I
        X    - least squares matrix n-by-p
        y    - right hand side vector n-by-1
        w    - weight vector n-by-1; or NULL for W = I
        Xs   - least squares matrix in standard form X~ n-by-p
        ys   - right hand side vector in standard form y~ n-by-1
        work - workspace

Return: success/error

Notes:
1) It is allowed for X = Xs and y = ys
*/

int
gsl_multilarge_linear_wstdform1 (const gsl_vector * L,
                                 const gsl_matrix * X,
                                 const gsl_vector * w,
                                 const gsl_vector * y,
                                 gsl_matrix * Xs,
                                 gsl_vector * ys,
                                 gsl_multilarge_linear_workspace * work)
{
  const size_t n = X->size1;
  const size_t p = X->size2;

  if (n > work->nmax)
    {
      GSL_ERROR("observation matrix larger than workspace", GSL_EBADLEN);
    }
  else if (p != work->p)
    {
      GSL_ERROR("observation matrix has different number of columns than workspace",
                GSL_EBADLEN);
    }
  else if (L != NULL && p != L->size)
    {
      GSL_ERROR("L vector does not match X", GSL_EBADLEN);
    }
  else if (n != y->size)
    {
      GSL_ERROR("y vector does not match X", GSL_EBADLEN);
    }
  else if (w != NULL && n != w->size)
    {
      GSL_ERROR("weight vector does not match X", GSL_EBADLEN);
    }
  else if (n != Xs->size1 || p != Xs->size2)
    {
      GSL_ERROR("Xs matrix dimensions do not match X", GSL_EBADLEN);
    }
  else if (n != ys->size)
    {
      GSL_ERROR("ys vector must be length n", GSL_EBADLEN);
    }
  else
    {
      int status = GSL_SUCCESS;

      /* compute Xs = sqrt(W) X and ys = sqrt(W) y */
      status = gsl_multifit_linear_applyW(X, w, y, Xs, ys);
      if (status)
        return status;

      if (L != NULL)
        {
          size_t j;

          /* construct X~ = sqrt(W) X * L^{-1} matrix */
          for (j = 0; j < p; ++j)
            {
              gsl_vector_view Xj = gsl_matrix_column(Xs, j);
              double lj = gsl_vector_get(L, j);

              if (lj == 0.0)
                {
                  GSL_ERROR("L matrix is singular", GSL_EDOM);
                }

              gsl_vector_scale(&Xj.vector, 1.0 / lj);
            }
        }

      return status;
    }
}

/*
gsl_multilarge_linear_genform1()
  Backtransform regularized solution vector using matrix
L = diag(L)
*/

int
gsl_multilarge_linear_genform1 (const gsl_vector * L,
                                const gsl_vector * cs,
                                gsl_vector * c,
                                gsl_multilarge_linear_workspace * work)
{
  if (L->size != work->p)
    {
      GSL_ERROR("L vector does not match workspace", GSL_EBADLEN);
    }
  else if (L->size != cs->size)
    {
      GSL_ERROR("cs vector does not match L", GSL_EBADLEN);
    }
  else if (L->size != c->size)
    {
      GSL_ERROR("c vector does not match L", GSL_EBADLEN);
    }
  else
    {
      /* compute true solution vector c = L^{-1} c~ */
      gsl_vector_memcpy(c, cs);
      gsl_vector_div(c, L);

      return GSL_SUCCESS;
    }
}
