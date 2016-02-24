/* multifit_nlinear/qr.c
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

/*
 * This module handles the solution of the linear least squares
 * system:
 *
 * [     J      ] v = - [ f ]
 * [ sqrt(mu)*D ]       [ 0 ]
 *
 * using a QR approach. The solver is organized into 4 sequential steps:
 *
 * 1. init: initialize solver for a given f(x) and J(x), independent of mu
 * 2. presolve: further solver initialization once mu is selected
 * 3. solve_vel: solve above linear system for geodesic velocity (this is the
 *               standard LM step)
 * 4. solve_acc: solve a similar linear system for the geodesic acceleration:
 *
 * [     J      ] a = - [ fvv ]
 * [ sqrt(mu)*D ]       [  0  ]
 *
 * only the right hand side is different in this case (instead of f(x) we
 * use the second directional derivative in the velocity direction).
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

#include "qrsolv.c"

typedef struct
{
  gsl_matrix *R;             /* QR factorization of J */
  gsl_vector *tau;           /* Householder scalars */
  gsl_matrix *Q;             /* Householder reflectors for J */
  gsl_permutation *perm;     /* permutation matrix */
  gsl_vector *qtf;           /* Q^T f */
  gsl_vector *diag;          /* scaling matrix D */
  gsl_vector *workn;         /* workspace, length n */
  gsl_vector *workp;         /* workspace, length p */
  double mu;                 /* LM parameter */
} qr_state_t;

static int qr_init(const gsl_matrix * J, void * vstate);
static int qr_presolve(const double mu, const gsl_vector * diag,
                       void * vstate);
static int qr_solve(const gsl_vector * f, const gsl_vector * g,
                    const gsl_matrix *J, gsl_vector *x, void *vstate);

static void *
qr_alloc (const size_t n, const size_t p)
{
  qr_state_t *state;

  (void)n;
  
  state = calloc(1, sizeof(qr_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate qr state", GSL_ENOMEM);
    }

  state->R = gsl_matrix_alloc(n, p);
  if (state->R == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for R", GSL_ENOMEM);
    }

  state->Q = gsl_matrix_alloc(n, p);
  if (state->Q == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for Q", GSL_ENOMEM);
    }

  state->tau = gsl_vector_alloc(p);
  if (state->tau == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for tau",
                      GSL_ENOMEM);
    }

  state->qtf = gsl_vector_alloc(n);
  if (state->qtf == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for qtf",
                      GSL_ENOMEM);
    }

  state->perm = gsl_permutation_calloc(p);
  if (state->perm == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for perm",
                      GSL_ENOMEM);
    }

  state->workn = gsl_vector_alloc(n);
  if (state->workn == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workn",
                      GSL_ENOMEM);
    }

  state->diag = gsl_vector_alloc(p);
  if (state->diag == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for diag",
                      GSL_ENOMEM);
    }

  state->workp = gsl_vector_alloc(p);
  if (state->workp == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workp",
                      GSL_ENOMEM);
    }

  state->mu = 0.0;

  return state;
}

static void
qr_free(void *vstate)
{
  qr_state_t *state = (qr_state_t *) vstate;

  if (state->R)
    gsl_matrix_free(state->R);

  if (state->Q)
    gsl_matrix_free(state->Q);

  if (state->tau)
    gsl_vector_free(state->tau);

  if (state->qtf)
    gsl_vector_free(state->qtf);

  if (state->perm)
    gsl_permutation_free(state->perm);

  if (state->diag)
    gsl_vector_free(state->diag);

  if (state->workn)
    gsl_vector_free(state->workn);

  if (state->workp)
    gsl_vector_free(state->workp);

  free(state);
}

/* compute J = Q R PT */
static int
qr_init(const gsl_matrix * J, void * vstate)
{
  qr_state_t *state = (qr_state_t *) vstate;
  int signum;

  gsl_matrix_memcpy(state->R, J);
  gsl_linalg_QRPT_decomp(state->R, state->tau, state->perm,
                         &signum, state->workp);

  /* save Householder part of R matrix which is destroyed by qrsolv() */
  gsl_matrix_memcpy(state->Q, state->R);

  return GSL_SUCCESS;
}

static int
qr_presolve(const double mu, const gsl_vector * diag,
            void * vstate)
{
  qr_state_t *state = (qr_state_t *) vstate;

  state->mu = mu;
  gsl_vector_memcpy(state->diag, diag);

  return GSL_SUCCESS;
}

static int
qr_solve(const gsl_vector * f, const gsl_vector * g,
         const gsl_matrix *J, gsl_vector *x, void *vstate)
{
  qr_state_t *state = (qr_state_t *) vstate;
  const double sqrt_mu = sqrt(state->mu);
  int status;

  (void)g;

  /* compute qtf = Q^T f */
  gsl_vector_memcpy(state->qtf, f);
  gsl_linalg_QR_QTvec(state->Q, state->tau, state->qtf);

  if (state->mu == 0.0)
    {
      /* solve: J x = f; this system can arise from the dogleg method */

      const size_t p = state->R->size2;
      gsl_matrix_view R = gsl_matrix_submatrix(state->R, 0, 0, p, p);
      gsl_vector_view qtf = gsl_vector_subvector(state->qtf, 0, p);

      status = gsl_linalg_QRPT_Rsolve(&R.matrix, state->perm,
                                      &qtf.vector, x);
    }
  else
    {
      /*
       * solve:
       *
       * [     J      ] x = [ f ]
       * [ sqrt(mu) D ]     [ 0 ]
       *
       * using QRPT factorization of J
       */

      status = qrsolv(state->R, state->perm, sqrt_mu, state->diag,
                      state->qtf, x, state->workp, state->workn);
    }

  /* reverse step to go downhill */
  gsl_vector_scale(x, -1.0);

  (void)J; /* avoid unused parameter warning */

  return status;
}

static const gsl_multifit_nlinear_solver qr_type =
{
  "qr",
  qr_alloc,
  qr_init,
  qr_presolve,
  qr_solve,
  qr_free
};

const gsl_multifit_nlinear_solver *gsl_multifit_nlinear_solver_qr = &qr_type;
