/* multifit_nlinear/lmtrust.c
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
 * damping parameter on each iteration using the trust region
 * method described in:
 *
 * [1] J. J. More, The Levenberg-Marquardt Algorithm: Implementation
 *     and Theory, Lecture Notes in Mathematics, v630, 1978.
 *
 * 5 routines are needed to implement the update procedure:
 *
 * 1. alloc  - allocate workspace state
 * 2. init   - initialize parameter prior to iteration
 * 3. accept - update parameter after a step has been accepted
 * 4. reject - update parameter after a step has been rejected
 * 5. free   - free workspace state
 */

typedef struct
{
  double delta; /* trust region size */
} trust_state_t;

static void *trust_alloc(void);
static void trust_free(void *vstate);
static int trust_init(const gsl_matrix * J, const gsl_vector * diag,
                      const gsl_vector * x, double * mu, void * vstate);
static int trust_accept(const double rho, double * mu, void * vstate);
static int trust_reject(double * mu, void * vstate);

static void *
trust_alloc(void)
{
  trust_state_t *state;

  state = calloc(1, sizeof(trust_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate trust state", GSL_ENOMEM);
    }

  return state;
}

static void
trust_free(void *vstate)
{
  trust_state_t *state = (trust_state_t *) vstate;
  free(state);
}

static int
trust_init(const gsl_matrix * J, const gsl_vector * diag,
           const gsl_vector * x, double * mu, void * vstate)
{
  trust_state_t *state = (trust_state_t *) vstate;
  const double factor = 100.0; /* recommended value from MINPACK */
  double Dx = scaled_norm(diag, x);

  /* initialize delta = factor * || D x || */
  state->delta = (Dx > 0.0) ? factor * Dx : factor;

  return GSL_SUCCESS;
}

static int
trust_accept(const double rho, double * mu, void * vstate)
{
  trust_state_t *state = (trust_state_t *) vstate;

  return GSL_SUCCESS;
}

static int
trust_reject(double * mu, void * vstate)
{
  trust_state_t *state = (trust_state_t *) vstate;

  return GSL_SUCCESS;
}

static const gsl_multifit_nlinear_update trust_type =
{
  "trust",
  trust_alloc,
  trust_init,
  trust_accept,
  trust_reject,
  trust_free
};

const gsl_multifit_nlinear_update *gsl_multifit_nlinear_update_trust = &trust_type;
