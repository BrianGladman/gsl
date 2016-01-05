/* multifit_nlinear/lmmisc.c
 * 
 * Copyright (C) 2014, 2015 Patrick Alken
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

static int lm_init_lambda(const gsl_matrix * J, lm_state_t * state);
static double lm_calc_rho(const double lambda, const gsl_vector * v,
                          const gsl_vector * g,
                          const gsl_vector * f,
                          const gsl_vector * f_trial,
                          lm_state_t * state);
static int lm_check_step(const gsl_vector * v, const gsl_vector * g,
                         const gsl_vector * f, const gsl_vector * f_trial,
                         double * rho, lm_state_t * state);
static double lm_scaled_norm(const gsl_vector *a, const gsl_vector *b,
                             gsl_vector *work);

/* initialize damping parameter lambda; state->diag must first be
 * initialized */
static int
lm_init_lambda(const gsl_matrix * J, lm_state_t * state)
{
  state->lambda = state->lambda0;

  if (state->init_diag == init_diag_levenberg)
    {
      /* when D = I, set lambda = lambda0 * max(diag(J^T J)) */

      const size_t p = J->size2;
      size_t j;
      double max = -1.0;

      for (j = 0; j < p; ++j)
        {
          gsl_vector_const_view v = gsl_matrix_const_column(J, j);
          double norm = gsl_blas_dnrm2(&v.vector);
          max = GSL_MAX(max, norm);
        }

      state->lambda *= max * max;
    }

  return GSL_SUCCESS;
}

/* compute x_trial = x + dx */
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
lm_calc_rho()
  Calculate ratio of actual reduction to predicted
reduction, given by Eq 4.4 of More, 1978.

Inputs: lambda  - LM parameter
        v       - velocity vector (p in Eq 4.4)
        g       - gradient J^T f
        f       - f(x)
        f_trial - f(x + dx)
        state   - workspace
*/

static double
lm_calc_rho(const double lambda, const gsl_vector * v,
            const gsl_vector * g, const gsl_vector * f,
            const gsl_vector * f_trial, lm_state_t * state)
{
  const double normf = gsl_blas_dnrm2(f);
  const double normf_trial = gsl_blas_dnrm2(f_trial);
  double rho;
  double actual_reduction;
  double pred_reduction;
  double u;
  double norm_Dp; /* || D p || */

  /* if ||f(x+dx)|| > ||f(x)|| reject step immediately */
  if (normf_trial >= normf)
    return -1.0;

  /* compute numerator of rho */
  u = normf_trial / normf;
  actual_reduction = 1.0 - u*u;

  /* compute || D p || */
  norm_Dp = lm_scaled_norm(state->diag, v, state->workp);

  /*
   * compute denominator of rho; instead of computing J*v,
   * we note that:
   *
   * ||Jv||^2 + 2*lambda*||Dv||^2 = lambda*||Dv||^2 - v^T g
   * and g = J^T f
   */
  u = norm_Dp / normf;
  pred_reduction = lambda * u * u;

  gsl_blas_ddot(v, g, &u);
  pred_reduction -= u / (normf * normf);

  if (pred_reduction > 0.0)
    rho = actual_reduction / pred_reduction;
  else
    rho = -1.0;

  return rho;
}

/*
lm_check_step()
  Check if a proposed step should be accepted or
rejected

Inputs: v       - proposed step velocity
        g       - gradient J^T f
        f       - f(x)
        f_trial - f(x + dx)
        rho     - (output) ratio of actual to predicted reduction
        state   - workspace

Return: GSL_SUCCESS to accept
        GSL_FAILURE to reject
*/

static int
lm_check_step(const gsl_vector * v, const gsl_vector * g,
              const gsl_vector * f, const gsl_vector * f_trial,
              double * rho, lm_state_t * state)
{
  /* if using geodesic acceleration, check that |a|/|v| < alpha */
  if (state->accel)
    {
      double anorm = lm_scaled_norm(state->diag, state->acc, state->workp);
      double vnorm = lm_scaled_norm(state->diag, state->vel, state->workp);
      double ratio = anorm / vnorm;

      /* reject step if acceleration is too large compared to velocity */
      if (ratio > state->accel_alpha)
        return GSL_FAILURE;
    }

  *rho = lm_calc_rho(state->lambda, v, g, f, f_trial, state);

  /* if rho <= 0, the step does not reduce the cost function, reject */
  if (*rho <= 0.0)
    return GSL_FAILURE;

  return GSL_SUCCESS;
}

/* compute || diag(D) a || */
static double
lm_scaled_norm(const gsl_vector *D, const gsl_vector *a,
               gsl_vector *work)
{
  const size_t n = a->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double Di = gsl_vector_get(D, i);
      double ai = gsl_vector_get(a, i);

      gsl_vector_set(work, i, Di * ai);
    }

  return gsl_blas_dnrm2(work);
}
