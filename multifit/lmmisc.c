/* multifit/lmmisc.c
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

#include "qrsolv.c"

/* JTJ = J^T J */
static int
lm_calc_JTJ(const gsl_matrix *J, gsl_matrix *JTJ)
{
  int status;
  status = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, J, J, 0.0, JTJ);
  return status;
} /* lm_calc_JTJ() */

static int
lm_calc_dx2(const double mu, const gsl_matrix *A, const gsl_vector *rhs,
            gsl_vector *dx, lm_state_t *state)
{
  int status;
  gsl_matrix *A_copy = state->A_copy;
  gsl_vector_view diag = gsl_matrix_diagonal(A_copy);

  /* make a copy of matrix */
  gsl_matrix_memcpy(A_copy, A);

  /* augment normal equations with LM parameter: A -> A + mu*I */
  gsl_vector_add_constant(&diag.vector, mu);

  status = gsl_linalg_QR_decomp(A_copy, state->work);
  if (status)
    return status;

  status = gsl_linalg_QR_solve(A_copy, state->work, rhs, dx);
  if (status)
    return status;

  return GSL_SUCCESS;
} /* lm_calc_dx2() */

/*
lm_calc_dx()
  Compute the least squares solution dx to the system:

[ J(x)       ] dx =~ [ -f ]
[ sqrt(mu) I ]       [  0 ]

where mu is the LM parameter. This least squares solution
is needed in each iteration of LM. J is decomposed into
J = Q R P^T and then the system which is solved is

[ R P^T      ] dx =~ [ -Q^T f ]
[ sqrt(mu) I ]       [    0   ]
*/

static int
lm_calc_dx(const double mu, const gsl_matrix *J, const gsl_vector *f,
            gsl_vector *dx, lm_state_t *state)
{
  int status;
  const double sqrt_mu = sqrt(mu);
  gsl_matrix *R = state->R;
  gsl_vector *tau = state->qrtau;
  gsl_permutation *perm = state->perm;
  gsl_vector *qtf = state->qtf;
  gsl_vector *diag = state->diag;
  gsl_vector *sdiag = state->sdiag;
  int signum;

  /* make a copy of matrix */
  gsl_matrix_memcpy(R, J);

  status = gsl_linalg_QRPT_decomp(R, tau, perm, &signum, state->work);
  if (status)
    return status;

  /* compute -Q^T f */
  gsl_vector_memcpy(qtf, f);
  gsl_vector_scale(qtf, -1.0);
  gsl_linalg_QR_QTvec(R, tau, qtf);

  status = qrsolv(R, perm, sqrt_mu, diag, qtf, dx, sdiag, state->work);
  if (status)
    return status;

  return GSL_SUCCESS;
} /* lm_calc_dx() */

static void
lm_trial_step(const gsl_vector * x, const gsl_vector * dx,
              gsl_vector * x_trial)
{
  size_t i, N = x->size;

  for (i = 0; i < N; i++)
    {
      double dxi = gsl_vector_get (dx, i);
      double xi = gsl_vector_get (x, i);
      gsl_vector_set (x_trial, i, xi + dxi);
    }
} /* lm_trial_step() */

/*
lm_calc_dF()
  Compute dF = F(x) - F(x + dx) = 1/2 (f - f_new)^T (f + f_new)
*/
static double
lm_calc_dF(const gsl_vector *f, const gsl_vector *f_new)
{
  const size_t N = f->size;
  size_t i;
  double dF = 0.0;

  for (i = 0; i < N; ++i)
    {
      double fi = gsl_vector_get(f, i);
      double fnewi = gsl_vector_get(f_new, i);

      dF += (fi - fnewi) * (fi + fnewi);
    }

  dF *= 0.5;

  return dF;
} /* lm_calc_dF() */

/*
lm_calc_dL()
  Compute dL = L(0) - L(dx) = 1/2 dx^T (mu * D^T D dx - g)
Here, the mg input is -g
*/

static double
lm_calc_dL(const double mu, const gsl_vector *diag,
           const gsl_vector *dx, const gsl_vector *mg)
{
  const size_t p = dx->size;
  size_t i;
  double dL = 0.0;

  for (i = 0; i < p; ++i)
    {
      double dxi = gsl_vector_get(dx, i);
      double di = gsl_vector_get(diag, i);
      double mgi = gsl_vector_get(mg, i); /* -g_i */

      dL += dxi * (mu * di * di * dxi + mgi);
    }

  dL *= 0.5;

  return dL;
} /* lm_calc_dL() */

/*
lm_calc_dLmore()
  Compute dL = L(0) - L(dx) = 1/2 dx^T (mu * D^T D dx - g)
Here, the mg input is -g
*/

static double
lm_calc_dLmore(const double mu, const gsl_vector *diag,
               const gsl_matrix *J, const gsl_vector *dx,
               const gsl_vector *f)
{
  const size_t p = dx->size;
  const size_t n = J->size1;
  gsl_vector *work1 = gsl_vector_alloc(n);
  gsl_vector *work2 = gsl_vector_alloc(p);
  double nf = gsl_blas_dnrm2(f);
  double nJp, nDp;
  double ratio1, ratio2;
  double dL = 0.0;

  /* compute ||Jp|| */
  gsl_blas_dgemv(CblasNoTrans, 1.0, J, dx, 0.0, work1);
  nJp = gsl_blas_dnrm2(work1);

  /* compute ||Dp|| */
  gsl_vector_memcpy(work2, dx);
  gsl_vector_mul(work2, diag);
  nDp = gsl_blas_dnrm2(work2);

  ratio1 = nJp / nf;
  ratio2 = nDp / nf;

  dL = ratio1*ratio1 + 2*mu*ratio2*ratio2;

  gsl_vector_free(work1);
  gsl_vector_free(work2);

  return dL;
} /* lm_calc_dLmore() */

static void
lm_update_diag(const gsl_matrix *J, gsl_vector *diag)
{
  size_t j, p = diag->size;

  for (j = 0; j < p; j++)
    {
      double *dj = gsl_vector_ptr (diag, j);
      gsl_vector_const_view v = gsl_matrix_const_column(J, j);
      double norm = gsl_blas_dnrm2(&v.vector);

      if (norm == 0.0)
        norm = 1.0;

      if (norm > *dj)
        *dj = norm;
    }
}
