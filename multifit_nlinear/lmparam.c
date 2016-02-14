/* multifit_nlinear/lmparam.c
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

/*
 * This module contains routines for updating the Levenberg-Marquardt
 * damping parameter on each iteration. Each update method needs 5
 * routines:
 *
 * 1. alloc  - allocate workspace state
 * 2. init   - initialize parameter prior to iteration
 * 3. accept - update parameter after a step has been accepted
 * 4. reject - update parameter after a step has been rejected
 * 5. free   - free workspace state
 *
 * There are currently two methods implemented, following the
 * publications:
 *
 * [1] H. B. Nielsen, K. Madsen, Introduction to Optimization and
 *     Data Fitting, Informatics and Mathematical Modeling,
 *     Technical University of Denmark (DTU), 2010.
 *
 * [2] J. J. More, The Levenberg-Marquardt Algorithm: Implementation
 *     and Theory, Lecture Notes in Mathematics, v630, 1978.
 */

#define LM_ONE_THIRD         (0.333333333333333)

typedef struct
{
  double mu0;  /* initial scale factor for mu */
  long nu;
} nielsen_state_t;

static void *nielsen_alloc(void);
static void nielsen_free(void *vstate);
static int nielsen_init(const gsl_matrix * J, const int diag,
                        double * mu, void * vstate);
static int nielsen_accept(const double rho, double * mu, void * vstate);
static int nielsen_reject(double * mu, void * vstate);

static void *
nielsen_alloc(void)
{
  nielsen_state_t *state;

  state = calloc(1, sizeof(nielsen_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate nielsen state", GSL_ENOMEM);
    }

  state->mu0 = 1.0e-3;

  return state;
}

static void
nielsen_free(void *vstate)
{
  nielsen_state_t *state = (nielsen_state_t *) vstate;
  free(state);
}

static int
nielsen_init(const gsl_matrix * J, const int diag,
             double * mu, void * vstate)
{
  nielsen_state_t *state = (nielsen_state_t *) vstate;

  state->nu = 2;
  *mu = state->mu0;

  if (!diag)
    {
      /* when D = I, set mu = mu0 * max(diag(J^T J)) */

      const size_t p = J->size2;
      size_t j;
      double max = -1.0;

      for (j = 0; j < p; ++j)
        {
          gsl_vector_const_view v = gsl_matrix_const_column(J, j);
          double norm = gsl_blas_dnrm2(&v.vector);
          max = GSL_MAX(max, norm);
        }

      *mu *= max * max;
    }

  return GSL_SUCCESS;
}

static int
nielsen_accept(const double rho, double * mu, void * vstate)
{
  nielsen_state_t *state = (nielsen_state_t *) vstate;
  double  b;
  
  /* reset nu */
  state->nu = 2;

  b = 2.0 * rho - 1.0;
  b = 1.0 - b*b*b;
  *mu *= GSL_MAX(LM_ONE_THIRD, b);

  return GSL_SUCCESS;
}

static int
nielsen_reject(double * mu, void * vstate)
{
  nielsen_state_t *state = (nielsen_state_t *) vstate;

  *mu *= (double) state->nu;

  /* nu := 2*nu */
  state->nu <<= 1;

  return GSL_SUCCESS;
}

static const gsl_multifit_nlinear_update nielsen_type =
{
  "nielsen",
  nielsen_alloc,
  nielsen_init,
  nielsen_accept,
  nielsen_reject,
  nielsen_free
};

const gsl_multifit_nlinear_update *gsl_multifit_nlinear_update_nielsen = &nielsen_type;
