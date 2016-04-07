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
 * using a QR approach.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permute_vector.h>

#include "common.c"
#include "qrsolv.c"

typedef struct
{
  gsl_matrix *R;             /* QR factorization of J */
  gsl_vector *tau;           /* Householder scalars */
  gsl_matrix *Q;             /* Householder reflectors for J */
  gsl_permutation *perm;     /* permutation matrix */
  gsl_vector *qtf;           /* Q^T f */
  gsl_vector *workn;         /* workspace, length n */
  gsl_vector *workp;         /* workspace, length p */
  double mu;                 /* LM parameter */
} qr_state_t;

static int qr_init(const void * vtrust_state, void * vstate);
static int qr_presolve(const double mu, const void * vtrust_state, void * vstate);
static int qr_solve(const gsl_vector * f, const gsl_vector * g,
                    gsl_vector *x, const void * vtrust_state, void *vstate);
static int qr_newton (const gsl_matrix * r, const gsl_permutation * perm,
                      const gsl_vector * qtf, gsl_vector * x);

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

  if (state->workn)
    gsl_vector_free(state->workn);

  if (state->workp)
    gsl_vector_free(state->workp);

  free(state);
}

/* compute J = Q R PT */
static int
qr_init(const void * vtrust_state, void * vstate)
{
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  qr_state_t *state = (qr_state_t *) vstate;
  int signum;

  gsl_matrix_memcpy(state->R, trust_state->J);
  gsl_linalg_QRPT_decomp(state->R, state->tau, state->perm,
                         &signum, state->workp);

  /* save Householder part of R matrix which is destroyed by qrsolv() */
  gsl_matrix_memcpy(state->Q, state->R);

  return GSL_SUCCESS;
}

static int
qr_presolve(const double mu, const void * vtrust_state, void * vstate)
{
  qr_state_t *state = (qr_state_t *) vstate;

  state->mu = mu;

  (void) vtrust_state;

  return GSL_SUCCESS;
}

static int
qr_solve(const gsl_vector * f, const gsl_vector * g,
         gsl_vector *x, const void * vtrust_state, void *vstate)
{
  qr_state_t *state = (qr_state_t *) vstate;
  int status;

  (void)g;

  /* compute qtf = Q^T f */
  gsl_vector_memcpy(state->qtf, f);
  gsl_linalg_QR_QTvec(state->Q, state->tau, state->qtf);

  if (state->mu == 0.0)
    {
      /*
       * compute Gauss-Newton direction by solving
       * J x = f
       */
      status = qr_newton(state->R, state->perm, state->qtf, x);
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

      const gsl_multifit_nlinear_trust_state *trust_state =
        (const gsl_multifit_nlinear_trust_state *) vtrust_state;
      double sqrt_mu = sqrt(state->mu);

      status = qrsolv(state->R, state->perm, sqrt_mu, trust_state->diag,
                      state->qtf, x, state->workp, state->workn);
    }

  /* reverse step to go downhill */
  gsl_vector_scale(x, -1.0);

  (void)vtrust_state; /* avoid unused parameter warning */

  return status;
}

/* compute Gauss-Newton direction */
static int
qr_newton (const gsl_matrix * r, const gsl_permutation * perm,
           const gsl_vector * qtf, gsl_vector * x)
{

  /* Compute and store in x the Gauss-Newton direction. If the
     Jacobian is rank-deficient then obtain a least squares
     solution. */

  const size_t n = r->size2;
  size_t i, j, nsing;

  for (i = 0; i < n; i++)
    {
      double qtfi = gsl_vector_get (qtf, i);
      gsl_vector_set (x, i,  qtfi);
    }

  nsing = qr_nonsing (r);

  for (i = nsing; i < n; i++)
    {
      gsl_vector_set (x, i, 0.0);
    }

  if (nsing > 0)
    {
      for (j = nsing; j > 0 && j--;)
        {
          double rjj = gsl_matrix_get (r, j, j);
          double temp = gsl_vector_get (x, j) / rjj;
          
          gsl_vector_set (x, j, temp);
          
          for (i = 0; i < j; i++)
            {
              double rij = gsl_matrix_get (r, i, j);
              double xi = gsl_vector_get (x, i);
              gsl_vector_set (x, i, xi - rij * temp);
            }
        }
    }

  gsl_permute_vector_inverse (perm, x);

  return GSL_SUCCESS;
}

static const gsl_multifit_nlinear_solver qr_type =
{
  "qr",
  qr_alloc,
  qr_init,
  qr_presolve,
  qr_solve,
  NULL, /* XXX */
  qr_free
};

const gsl_multifit_nlinear_solver *gsl_multifit_nlinear_solver_qr = &qr_type;
