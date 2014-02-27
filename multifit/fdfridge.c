/* multifit/fdfridge.c
 * 
 * Copyright (C) 2014 Patrick Alken
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

static int fdfridge_f(const gsl_vector * x, void * params, gsl_vector * f);
static int fdfridge_df(const gsl_vector * x, void * params, gsl_matrix * J);

gsl_multifit_fdfridge *
gsl_multifit_fdfridge_alloc (const gsl_multifit_fdfsolver_type * T,
                             const size_t n, const size_t p)
{
  gsl_multifit_fdfridge * work;

  work = calloc(1, sizeof(gsl_multifit_fdfridge));
  if (work == NULL)
    {
      GSL_ERROR_VAL("failed to allocate workspace",
                    GSL_ENOMEM, 0);
    }

  work->s = gsl_multifit_fdfsolver_alloc (T, n + p, p);
  if (work->s == NULL)
    {
      gsl_multifit_fdfridge_free(work);
      GSL_ERROR_VAL("failed to allocate space for fdfsolver",
                    GSL_ENOMEM, 0);
    }

  work->n = n;
  work->p = p;
  work->lambda = 0.0;

  return work;
} /* gsl_multifit_fdfridge_alloc() */

void
gsl_multifit_fdfridge_free(gsl_multifit_fdfridge *work)
{
  if (work->s)
    gsl_multifit_fdfsolver_free(work->s);

  free(work);
}

const char *
gsl_multifit_fdfridge_name(const gsl_multifit_fdfridge * w)
{
  return gsl_multifit_fdfsolver_name(w->s);
}

gsl_vector *
gsl_multifit_fdfridge_position (const gsl_multifit_fdfridge * w)
{
  return gsl_multifit_fdfsolver_position(w->s);
}

size_t
gsl_multifit_fdfridge_niter (const gsl_multifit_fdfridge * w)
{
  return w->s->niter;
}

int
gsl_multifit_fdfridge_set (gsl_multifit_fdfridge * w,
                           gsl_multifit_function_fdf * f,
                           const gsl_vector * x,
                           const double lambda)
{
  return gsl_multifit_fdfridge_wset(w, f, x, lambda, NULL);
} /* gsl_multifit_fdfridge_set() */

int
gsl_multifit_fdfridge_wset (gsl_multifit_fdfridge * w,
                            gsl_multifit_function_fdf * f,
                            const gsl_vector * x,
                            const double lambda,
                            const gsl_vector * wts)
{
  if (w->n != f->n || w->p != f->p)
    {
      GSL_ERROR ("function size does not match solver", GSL_EBADLEN);
    }
  else if (w->p != x->size)
    {
      GSL_ERROR ("vector length does not match solver", GSL_EBADLEN);
    }
  else if (wts != NULL && w->n != wts->size)
    {
      GSL_ERROR ("weight vector length does not match solver", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* save user defined fdf */
      w->fdf = f;

      /* build modified fdf for Tikhonov terms */
      w->fdftik.f = &fdfridge_f;
      w->fdftik.df = &fdfridge_df;
      w->fdftik.n = w->n + w->p; /* add p for Tikhonov terms */
      w->fdftik.p = w->p;
      w->fdftik.params = (void *) w;

      /* store damping parameter */
      w->lambda = lambda;
      w->L = NULL;

      status = gsl_multifit_fdfsolver_wset(w->s, &(w->fdftik), x, wts);

      return status;
    }
} /* gsl_multifit_fdfridge_wset() */

int
gsl_multifit_fdfridge_set2 (gsl_multifit_fdfridge * w,
                            gsl_multifit_function_fdf * f,
                            const gsl_vector * x,
                            const gsl_matrix * L)
{
  return gsl_multifit_fdfridge_wset2(w, f, x, L, NULL);
} /* gsl_multifit_fdfridge_set2() */

int
gsl_multifit_fdfridge_wset2 (gsl_multifit_fdfridge * w,
                             gsl_multifit_function_fdf * f,
                             const gsl_vector * x,
                             const gsl_matrix * L,
                             const gsl_vector * wts)
{
  if (w->n != f->n || w->p != f->p)
    {
      GSL_ERROR ("function size does not match solver", GSL_EBADLEN);
    }
  else if (w->p != x->size)
    {
      GSL_ERROR ("vector length does not match solver", GSL_EBADLEN);
    }
  else if (L->size2 != w->p)
    {
      GSL_ERROR ("L matrix size2 does not match solver", GSL_EBADLEN);
    }
  else if (wts != NULL && w->n != wts->size)
    {
      GSL_ERROR ("weight vector length does not match solver", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* save user defined fdf */
      w->fdf = f;
      w->fdf->nevalf = 0;
      w->fdf->nevaldf = 0;

      /* build modified fdf for Tikhonov terms */
      w->fdftik.f = &fdfridge_f;
      w->fdftik.df = &fdfridge_df;
      w->fdftik.n = w->n + w->p; /* add p for Tikhonov terms */
      w->fdftik.p = w->p;
      w->fdftik.params = (void *) w;

      /* store damping matrix */
      w->lambda = 0.0;
      w->L = L;

      status = gsl_multifit_fdfsolver_wset(w->s, &(w->fdftik), x, wts);

      return status;
    }
} /* gsl_multifit_fdfridge_wset2() */

int
gsl_multifit_fdfridge_iterate (gsl_multifit_fdfridge * w)
{
  int status = gsl_multifit_fdfsolver_iterate(w->s);
  return status;
}

int
gsl_multifit_fdfridge_driver (gsl_multifit_fdfridge * w,
                              const size_t maxiter,
                              const double xtol,
                              const double gtol,
                              const double ftol,
                              int *info)
{
  int status = gsl_multifit_fdfsolver_driver(w->s, maxiter, xtol,
                                             gtol, ftol, info);
  return status;
} /* gsl_multifit_fdfridge_driver() */

/*
fdfridge_f()
  Callback function to provide residuals, including extra p
Tikhonov terms. The residual vector will have the form:

f~ = [     f     ]
     [ \lambda x ]

where f is the user supplied residuals, x are the model
parameters, and \lambda is the Tikhonov damping parameter

Inputs: x      - model parameters (size p)
        params - pointer to fdfridge workspace
        f      - (output) (n+p) vector to store f~
*/

static int
fdfridge_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  int status;
  gsl_multifit_fdfridge *w = (gsl_multifit_fdfridge *) params;
  const size_t n = w->n;
  const size_t p = w->p;
  gsl_vector_view f_user = gsl_vector_subvector(f, 0, n);
  gsl_vector_view f_tik = gsl_vector_subvector(f, n, p);

  /* call user callback function to get residual vector f */
  status = gsl_multifit_eval_wf(w->fdf, x, NULL, &f_user.vector);
  if (status)
    return status;

  if (w->L)
    {
      /* store Lx in Tikhonov portion of f~ */
      gsl_blas_dgemv(CblasNoTrans, 1.0, w->L, x, 0.0, &f_tik.vector);
    }
  else
    {
      /* store \lambda x in Tikhonov portion of f~ */
      gsl_vector_memcpy(&f_tik.vector, x);
      gsl_vector_scale(&f_tik.vector, w->lambda);
    }

  return GSL_SUCCESS;
} /* fdfridge_f() */

static int
fdfridge_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  int status;
  gsl_multifit_fdfridge *w = (gsl_multifit_fdfridge *) params;
  const size_t n = w->n;
  const size_t p = w->p;
  gsl_matrix_view J_user = gsl_matrix_submatrix(J, 0, 0, n, p);
  gsl_matrix_view J_tik = gsl_matrix_submatrix(J, n, 0, p, p);
  gsl_vector_view diag = gsl_matrix_diagonal(&J_tik.matrix);

  /* compute user supplied Jacobian */
  status = gsl_multifit_eval_wdf(w->fdf, x, NULL, &J_user.matrix);
  if (status)
    return status;

  if (w->L)
    {
      /* store L in Tikhonov portion of J */
      gsl_matrix_memcpy(&J_tik.matrix, w->L);
    }
  else
    {
      /* store \lambda I_p in Tikhonov portion of J */
      gsl_matrix_set_zero(&J_tik.matrix);
      gsl_vector_set_all(&diag.vector, w->lambda);
    }

  return GSL_SUCCESS;
} /* fdfridge_df() */
