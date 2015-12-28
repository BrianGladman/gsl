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
static double lm_calc_rho(const double lambda, const gsl_vector * dx,
                          const gsl_vector * minus_g,
                          const gsl_vector * f,
                          const gsl_vector * f_trial,
                          lm_state_t * state);

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
*/

static double
lm_calc_rho(const double lambda, const gsl_vector * dx,
            const gsl_vector * minus_g, const gsl_vector * f,
            const gsl_vector * f_trial, lm_state_t * state)
{
  const double normf = gsl_blas_dnrm2(f);
  const double normf_trial = gsl_blas_dnrm2(f_trial);
  double rho;
  double actual_reduction;
  double pred_reduction;
  double u;
  double norm_Ddx; /* || D dx || */
  size_t i;

  /* if ||f(x+dx)|| > ||f(x)|| reject step immediately */
  if (normf_trial >= normf)
    return -1.0;

  /* compute numerator of rho */
  u = normf_trial / normf;
  actual_reduction = 1.0 - u*u;

  /* compute || D dx || */
  for (i = 0; i < dx->size; ++i)
    {
      double dxi = gsl_vector_get(dx, i);
      double di = gsl_vector_get(state->diag, i);

      gsl_vector_set(state->workp, i, dxi * di);
    }

  norm_Ddx = gsl_blas_dnrm2(state->workp);

  /* compute denominator of rho */
  u = norm_Ddx / normf;
  pred_reduction = lambda * u * u;

  for (i = 0; i < dx->size; ++i)
    {
      double dxi = gsl_vector_get(dx, i);
      double mgi = gsl_vector_get(minus_g, i);

      pred_reduction += (dxi / normf) * (mgi / normf);
    }

  if (pred_reduction > 0.0)
    rho = actual_reduction / pred_reduction;
  else
    rho = -1.0;

  return rho;
}
