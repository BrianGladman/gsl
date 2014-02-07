/* multifit/lm.c
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

typedef struct
{
  gsl_matrix *A;         /* J^T J */
  gsl_matrix *A_copy;    /* copy of J^T J */
  gsl_vector *g;         /* -J^T f */
  gsl_vector *x_trial;   /* trial parameter vector */
  gsl_vector *f_trial;   /* trial function vector */
  gsl_vector *qrwork;    /* QR workspace */
  long nu;               /* nu */
  double mu;             /* LM damping parameter mu */
  double tau;            /* initial scale factor for mu */
  double gtol;           /* gradient tolerance */
  double xtol;           /* parameter tolerance */
  double ftol;           /* residual tolerance */
} lm_state_t;

static int lm_alloc (void *vstate, const size_t n, const size_t p);
static void lm_free(void *vstate);
static int lm_set(void *vstate, gsl_multifit_function_fdf *fdf,
                  gsl_vector *x, gsl_vector *f, gsl_matrix *J,
                  gsl_vector *dx);
static int lm_iterate(void *vstate, gsl_multifit_function_fdf *fdf,
                      gsl_vector *x, gsl_vector *f, gsl_matrix *J,
                      gsl_vector *dx);
static int lm_calc_JTJ(const gsl_matrix *J, gsl_matrix *JTJ);
static int lm_calc_JTf(const gsl_matrix *J, const gsl_vector *f,
                       gsl_vector *JTf);
static double lm_infnorm(const gsl_vector *v);
static void lm_trial_step(const gsl_vector * x, const gsl_vector * dx,
                          gsl_vector * x_trial);
static double lm_calc_dL(const double mu, const gsl_vector *dx,
                         const gsl_vector *g);
static double lm_calc_df_dot(const gsl_vector *f1, const gsl_vector *f2);

#define LM_ONE_THIRD         (0.333333333333333)

static int
lm_alloc (void *vstate, const size_t n, const size_t p)
{
  lm_state_t *state = (lm_state_t *) vstate;
  double eps;

  state->A = gsl_matrix_alloc(p, p);
  if (state->A == NULL)
    {
      GSL_ERROR ("failed to allocate space for A", GSL_ENOMEM);
    }

  state->g = gsl_vector_alloc(p);
  if (state->g == NULL)
    {
      GSL_ERROR ("failed to allocate space for g", GSL_ENOMEM);
    }

  state->qrwork = gsl_vector_alloc(p);
  if (state->qrwork == NULL)
    {
      GSL_ERROR ("failed to allocate space for qrwork", GSL_ENOMEM);
    }

  state->A_copy = gsl_matrix_alloc(p, p);
  if (state->A_copy == NULL)
    {
      GSL_ERROR ("failed to allocate space for A_copy", GSL_ENOMEM);
    }

  state->x_trial = gsl_vector_alloc(p);
  if (state->x_trial == NULL)
    {
      GSL_ERROR ("failed to allocate space for x_trial", GSL_ENOMEM);
    }

  state->f_trial = gsl_vector_alloc(n);
  if (state->f_trial == NULL)
    {
      GSL_ERROR ("failed to allocate space for f_trial", GSL_ENOMEM);
    }

  /* initialize convergence parameters */
  eps = pow(GSL_DBL_EPSILON, 0.9);
  state->xtol = eps;
  state->ftol = eps;
  state->gtol = eps;

  state->tau = 1.0e-3;

  return GSL_SUCCESS;
} /* lm_alloc() */

static void
lm_free(void *vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;

  if (state->A)
    gsl_matrix_free(state->A);

  if (state->g)
    gsl_vector_free(state->g);

  if (state->qrwork)
    gsl_vector_free(state->qrwork);

  if (state->A_copy)
    gsl_matrix_free(state->A_copy);

  if (state->x_trial)
    gsl_vector_free(state->x_trial);

  if (state->f_trial)
    gsl_vector_free(state->f_trial);
} /* lm_free() */

static int
lm_set(void *vstate, gsl_multifit_function_fdf *fdf, gsl_vector *x,
       gsl_vector *f, gsl_matrix *J, gsl_vector *dx)
{
  lm_state_t *state = (lm_state_t *) vstate;

  /* set default parameters */
  state->nu = 2;
  state->mu = -1.0;

  return GSL_SUCCESS;
} /* lm_set() */

static int
lm_iterate(void *vstate, gsl_multifit_function_fdf *fdf, gsl_vector *x,
           gsl_vector *f, gsl_matrix *J, gsl_vector *dx)
{
  int status;
  lm_state_t *state = (lm_state_t *) vstate;
  gsl_matrix *A = state->A;                      /* J^T J */
  gsl_matrix *A_copy = state->A_copy;            /* copy of J^T J */
  gsl_vector_view diag = gsl_matrix_diagonal(A); /* diag(J^T J) */
  gsl_vector *g = state->g;                      /* -J^T f */
  gsl_vector *qrwork = state->qrwork;            /* QR workspace */
  gsl_vector *x_trial = state->x_trial;          /* trial parameters */
  gsl_vector *f_trial = state->f_trial;          /* trial function */
  double g_inf;                                  /* ||g||_inf */
  double x_2;                                    /* ||x||_2 */
  double dx_2;                                   /* ||dx||_2 */
  double f_sq;                                   /* ||f||^2 */
  double f_trial_sq;                             /* ||f_trial||^2 */
  double dL;                                     /* L(0) - L(dx) */
  double rho;

  /* evaluate function and Jacobian at x */
  status = GSL_MULTIFIT_FN_EVAL_F_DF(fdf, x, f, J);
  if (status)
    return status;

  /* compute A = J^T J */
  status = lm_calc_JTJ(J, A);
  if (status)
    return status;

  /* compute g = -J^T f */
  status = lm_calc_JTf(J, f, g);
  if (status)
    return status;

  /* compute ||f||^2 */
  gsl_blas_ddot(f, f, &f_sq);

  if (state->mu < 0.0)
    {
      /* first iteration, initialize mu */
      state->mu = state->tau * gsl_vector_max(&diag.vector);
    }

  /* check ||J^T f||_inf <= gtol */
  g_inf = lm_infnorm(g);
  if (g_inf <= state->gtol)
    {
      fprintf(stderr, "GTOL FLAG\n");
      return GSL_SUCCESS; /* converged */
    }

  /* compute ||x||_2 */
  x_2 = gsl_blas_dnrm2(x);

  /* save J^T J in case we need multiple decomps of augmented matrix */
  gsl_matrix_memcpy(A_copy, A);

  /* inner loop starts here */
  do
    {
      /* augment normal equations with LM parameter: A -> A + mu*I */
      gsl_vector_add_constant(&diag.vector, state->mu);

      /* solve (A + mu*I) dx = g */
      status = gsl_linalg_QR_decomp(A, qrwork);
      if (status)
        return status;

      status = gsl_linalg_QR_solve(A, qrwork, g, dx);
      if (status)
        return status;

      /* compute ||dx||_2 and check for convergence */
      dx_2 = gsl_blas_dnrm2(dx);
      if (dx_2 <= state->xtol * x_2)
        {
          fprintf(stderr, "XTOL FLAG\n");
          return GSL_SUCCESS; /* converged */
        }

      /* compute x_trial = x + dx */
      lm_trial_step(x, dx, x_trial);

      /* evaluate function at x + dx */
      status = GSL_MULTIFIT_FN_EVAL_F (fdf, x_trial, f_trial);
      if (status)
       return status;

      /* compute ||f_trial||^2 */
      gsl_blas_ddot(f_trial, f_trial, &f_trial_sq);

      /* compute dL = L(0) - L(dx) = dx^T (mu*dx - g) */
      dL = lm_calc_dL(state->mu, dx, g);

      /* compute rho = (||f||^2 - ||f_trial||^2) / (L(0) - L(dx)) */
      rho = (f_sq - f_trial_sq) / dL;
      if (rho > 0.0)
        {
          /* reduction in error, step acceptable */

          double tmp, df_sq;

          /* compute || f(x+dx) - f(x) ||^2 */
          df_sq = lm_calc_df_dot(f, f_trial);

          /* update x = x_trial */
          gsl_vector_memcpy(x, x_trial);

          /* update f = f_trial */
          gsl_vector_memcpy(f, f_trial);

          if (df_sq <= state->ftol * GSL_MAX(GSL_MAX(f_sq, f_trial_sq), 1.0))
            {
              fprintf(stderr, "FTOL FLAG\n");
              return GSL_SUCCESS; /* converged */
            }

          f_sq = f_trial_sq;

          /* update LM parameter mu */
          tmp = 2.0 * rho - 1.0;
          tmp = 1.0 - tmp*tmp*tmp;
          state->mu *= GSL_MAX(LM_ONE_THIRD, tmp);
          state->nu = 2;
        }
      else
        {
          long nu2;

          /* step did not reduce error, reject step */
          state->mu *= state->nu;
          nu2 = state->nu << 1; /* 2*nu */
          if (nu2 <= state->nu)
            {
              /* nu has wrapped around / overflown */
              GSL_ERROR("nu parameter has overflown", GSL_EOVRFLW);
            }
          state->nu = nu2;

          /* restore J^T J for next inner loop iteration */
          gsl_matrix_memcpy(A, A_copy);
        }
    }
  while (rho <= 0.0);

  /* not yet converged */
  return GSL_CONTINUE;
} /* lm_iterate() */

/* JTJ = J^T J */
static int
lm_calc_JTJ(const gsl_matrix *J, gsl_matrix *JTJ)
{
  int status;
  status = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, J, J, 0.0, JTJ);
  return status;
} /* lm_calc_JTJ() */

/* JTf = -J^T f (note minus sign) */
static int
lm_calc_JTf(const gsl_matrix *J, const gsl_vector *f, gsl_vector *JTf)
{
  int status;
  status = gsl_blas_dgemv(CblasTrans, -1.0, J, f, 0.0, JTf);
  return status;
} /* lm_calc_JTf() */

static double
lm_infnorm(const gsl_vector *v)
{
  CBLAS_INDEX_t idx = gsl_blas_idamax(v);
  return fabs(gsl_vector_get(v, idx));
}

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

static double
lm_calc_dL(const double mu, const gsl_vector *dx, const gsl_vector *g)
{
  const size_t N = dx->size;
  size_t i;
  double dL = 0.0;

  for (i = 0; i < N; ++i)
    {
      double dxi = gsl_vector_get(dx, i);
      double gi = gsl_vector_get(g, i);

      dL += dxi * (mu * dxi + gi);
    }

  return dL;
} /* lm_calc_dL() */

/* compute || f1 - f2 ||^2 */
static double
lm_calc_df_dot(const gsl_vector *f1, const gsl_vector *f2)
{
  const size_t N = f1->size;
  size_t i;
  double sum = 0.0;

  for (i = 0; i < N; ++i)
    {
      double f1i = gsl_vector_get(f1, i);
      double f2i = gsl_vector_get(f2, i);
      sum += (f2i - f1i) * (f2i - f1i);
    }

  return sum;
}

static const gsl_multifit_fdfsolver_type lm_type =
{
  "lmsder",
  sizeof(lm_state_t),
  &lm_alloc,
  &lm_set,
  &lm_iterate,
  &lm_free
};

const gsl_multifit_fdfsolver_type *gsl_multifit_fdfsolver_lmsder = &lm_type;
