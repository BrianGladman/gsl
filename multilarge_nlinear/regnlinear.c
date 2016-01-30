/* multilarge_nlinear/regnlinear.c
 * 
 * Copyright (C) 2016 Patrick Alken
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
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_blas.h>

static int regnlinear_f(const gsl_vector * x, void * params, gsl_vector * f);
static int regnlinear_df(const gsl_vector * x, const gsl_vector * y, void * params,
                         gsl_vector * JTy, gsl_matrix * JTJ);

gsl_multilarge_regnlinear_workspace *
gsl_multilarge_regnlinear_alloc (const gsl_multilarge_nlinear_type * T,
                                 const gsl_multilarge_nlinear_parameters * params,
                                 const size_t n, const size_t p, const size_t m)
{
  gsl_multilarge_regnlinear_workspace * work;

  work = calloc(1, sizeof(gsl_multilarge_regnlinear_workspace));
  if (work == NULL)
    {
      GSL_ERROR_VAL("failed to allocate workspace",
                    GSL_ENOMEM, 0);
    }

  work->diag = gsl_vector_alloc(p);
  if (work->diag == NULL)
    {
      gsl_multilarge_regnlinear_free(work);
      GSL_ERROR_VAL("failed to allocate space for diag", GSL_ENOMEM, 0);
    }

  work->multilarge_nlinear_p = gsl_multilarge_nlinear_alloc (T, params, n + m, p);
  if (work->multilarge_nlinear_p == NULL)
    {
      gsl_multilarge_regnlinear_free(work);
      GSL_ERROR_VAL("failed to allocate space for multilarge workspace",
                    GSL_ENOMEM, 0);
    }

  work->n = n;
  work->p = p;
  work->m = m;
  work->lambda = 0.0;

  return work;
}

void
gsl_multilarge_regnlinear_free(gsl_multilarge_regnlinear_workspace *w)
{
  if (w->diag)
    gsl_vector_free(w->diag);

  if (w->multilarge_nlinear_p)
    gsl_multilarge_nlinear_free(w->multilarge_nlinear_p);

  free(w);
}

const char *
gsl_multilarge_regnlinear_name(const gsl_multilarge_regnlinear_workspace * w)
{
  return gsl_multilarge_nlinear_name(w->multilarge_nlinear_p);
}

gsl_vector *
gsl_multilarge_regnlinear_position (const gsl_multilarge_regnlinear_workspace * w)
{
  return gsl_multilarge_nlinear_position(w->multilarge_nlinear_p);
}

gsl_vector *
gsl_multilarge_regnlinear_residual (const gsl_multilarge_regnlinear_workspace * w)
{
  return gsl_multilarge_nlinear_residual(w->multilarge_nlinear_p);
}

gsl_matrix *
gsl_multilarge_regnlinear_JTJ (const gsl_multilarge_regnlinear_workspace * w)
{
  return gsl_multilarge_nlinear_JTJ(w->multilarge_nlinear_p);
}

size_t
gsl_multilarge_regnlinear_niter (const gsl_multilarge_regnlinear_workspace * w)
{
  return gsl_multilarge_nlinear_niter(w->multilarge_nlinear_p);
}

int
gsl_multilarge_regnlinear_init2 (const double lambda,
                                 const gsl_vector * L,
                                 const gsl_vector * x,
                                 gsl_multilarge_nlinear_fdf * f,
                                 gsl_multilarge_regnlinear_workspace * w)
{
  const size_t n = w->n;
  const size_t p = w->p;
  const size_t m = w->m;

  if (n != f->n || p != f->p)
    {
      GSL_ERROR ("function size does not match workspace", GSL_EBADLEN);
    }
  else if (p != x->size)
    {
      GSL_ERROR ("vector length does not match workspace", GSL_EBADLEN);
    }
  else if (L->size != p)
    {
      GSL_ERROR ("L vector size does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* save user defined fdf */
      w->fdf = f;
      w->fdf->nevalf = 0;
      w->fdf->nevaldf = 0;

      /* build modified fdf for Tikhonov terms */
      w->fdftik.f = regnlinear_f;
      w->fdftik.df = regnlinear_df;
      w->fdftik.fvv = NULL;
      w->fdftik.n = n + m; /* add m for Tikhonov terms */
      w->fdftik.p = p;
      w->fdftik.params = (void *) w;

      /* store damping matrix */
      w->lambda = lambda;
      gsl_vector_memcpy(w->diag, L);

      status = gsl_multilarge_nlinear_init(x, &(w->fdftik), w->multilarge_nlinear_p);

      /* update function/Jacobian evaluations */
      f->nevalf = w->fdftik.nevalf;
      f->nevaldf = w->fdftik.nevaldf;
      f->nevaldff = w->fdftik.nevaldff;

      return status;
    }
}

int
gsl_multilarge_regnlinear_iterate (gsl_multilarge_regnlinear_workspace * w)
{
  int status = gsl_multilarge_nlinear_iterate(w->multilarge_nlinear_p);

  /* update function/Jacobian evaluations */
  w->fdf->nevalf = w->fdftik.nevalf;
  w->fdf->nevaldf = w->fdftik.nevaldf;

  return status;
}

int
gsl_multilarge_regnlinear_driver (const size_t maxiter,
                                  const double xtol,
                                  const double gtol,
                                  const double ftol,
                                  void (*callback)(const size_t iter, void *params,
                                                   const gsl_multilarge_nlinear_workspace *w),
                                  void *callback_params,
                                  int *info,
                                  gsl_multilarge_regnlinear_workspace * w)
{
  int status = gsl_multilarge_nlinear_driver(maxiter, xtol, gtol, ftol,
                                             callback, callback_params, info,
                                             w->multilarge_nlinear_p);
  return status;
}

/*
regnlinear_f()
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
regnlinear_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  int status;
  gsl_multilarge_regnlinear_workspace *w = (gsl_multilarge_regnlinear_workspace *) params;
  const size_t n = w->n;
  const size_t p = w->p;
  const size_t m = w->m;
  gsl_vector_view f_user = gsl_vector_subvector(f, 0, n);
  gsl_vector_view f_tik = gsl_vector_subvector(f, n, m);

  /* call user callback function to get residual vector f */
  status = gsl_multilarge_nlinear_eval_f(w->fdf, x, &f_user.vector);
  if (status)
    return status;

  if (w->diag)
    {
      size_t i;

      /* store lambda*D*x in Tikhonov portion of f~ */

      for (i = 0; i < p; ++i)
        {
          double xi = gsl_vector_get(x, i);
          double di = gsl_vector_get(w->diag, i);

          gsl_vector_set(&f_tik.vector, i, w->lambda * di * xi);
        }
    }

  return GSL_SUCCESS;
}

static int
regnlinear_df(const gsl_vector * x, const gsl_vector * y, void * params,
              gsl_vector * JTy, gsl_matrix * JTJ)
{
  int status;
  gsl_multilarge_regnlinear_workspace *w = (gsl_multilarge_regnlinear_workspace *) params;
  const size_t n = w->n;
  const size_t p = w->p;
  const size_t m = w->m;
  gsl_vector_const_view y_user = gsl_vector_const_subvector(y, 0, n);
  gsl_vector_const_view y_tik = gsl_vector_const_subvector(y, n, m);

  /* compute Jacobian */
  status = gsl_multilarge_nlinear_eval_df(w->fdf, x, &y_user.vector, JTy, JTJ);
  if (status)
    return status;

  if (w->diag)
    {
      size_t i;

      /*
       * JTJ += lambda^2 D^T D 
       * JTy += lambda D^T y_tik
       */
      for (i = 0; i < p; ++i)
        {
          double *Aii = gsl_matrix_ptr(JTJ, i, i);
          double *JTyi = gsl_vector_ptr(JTy, i);
          double di = gsl_vector_get(w->diag, i);
          double yi = gsl_vector_get(&y_tik.vector, i);

          *Aii += w->lambda * w->lambda * di * di;
          *JTyi += w->lambda * di * yi;
        }
    }

  return GSL_SUCCESS;
}
