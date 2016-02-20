/* multifit_nlinear/common.c
 * 
 * Copyright (C) 2014, 2015, 2016 Patrick Alken
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

static void trial_step(const gsl_vector * x, const gsl_vector * dx,
                       gsl_vector * x_trial);
static double scaled_norm(const gsl_vector *a, const gsl_vector *b);

/* compute x_trial = x + dx */
static void
trial_step(const gsl_vector * x, const gsl_vector * dx,
           gsl_vector * x_trial)
{
  size_t i, N = x->size;

  for (i = 0; i < N; i++)
    {
      double dxi = gsl_vector_get (dx, i);
      double xi = gsl_vector_get (x, i);
      gsl_vector_set (x_trial, i, xi + dxi);
    }
}

#if 0
/*
lm_calc_rho()
  Calculate ratio of actual reduction to predicted
reduction, given by Eq 4.4 of More, 1978.

Inputs: mu      - LM parameter
        v       - velocity vector (p in Eq 4.4)
        g       - gradient J^T f
        f       - f(x)
        f_trial - f(x + dx)
        state   - workspace
*/

static double
lm_calc_rho(const double mu, const gsl_vector * v,
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
  norm_Dp = lm_scaled_norm(state->diag, v);

  /*
   * compute denominator of rho; instead of computing J*v,
   * we note that:
   *
   * ||Jv||^2 + 2*mu*||Dv||^2 = mu*||Dv||^2 - v^T g
   * and g = J^T f
   */
  u = norm_Dp / normf;
  pred_reduction = mu * u * u;

  gsl_blas_ddot(v, g, &u);
  pred_reduction -= u / (normf * normf);

  if (pred_reduction > 0.0)
    rho = actual_reduction / pred_reduction;
  else
    rho = -1.0;

  return rho;
}
#endif

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

Notes:
1) If using geodesic acceleration, state->avratio
   is updated with |a| / |v|
*/

#if 0
static int
lm_check_step(const gsl_vector * v, const gsl_vector * g,
              const gsl_vector * f, const gsl_vector * f_trial,
              double * rho, lm_state_t * state)
{
  const gsl_multifit_nlinear_parameters *params = &(state->params);

  /* if using geodesic acceleration, check that |a|/|v| < alpha */
  if (params->accel)
    {
      double anorm = lm_scaled_norm(state->diag, state->acc);
      double vnorm = lm_scaled_norm(state->diag, state->vel);

      /* store |a| / |v| */
      state->avratio = anorm / vnorm;

      /* reject step if acceleration is too large compared to velocity */
      if (state->avratio > params->avmax)
        return GSL_FAILURE;
    }

  *rho = lm_calc_rho(state->mu, v, g, f, f_trial, state);

  /* if rho <= 0, the step does not reduce the cost function, reject */
  if (*rho <= 0.0)
    return GSL_FAILURE;

  return GSL_SUCCESS;
}
#endif

/* compute || diag(D) a || */
static double
scaled_norm(const gsl_vector *D, const gsl_vector *a)
{
  const size_t n = a->size;
  double e2 = 0.0;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double Di = gsl_vector_get(D, i);
      double ai = gsl_vector_get(a, i);
      double u = Di * ai;

      e2 += u * u;
    }

  return sqrt (e2);
}
