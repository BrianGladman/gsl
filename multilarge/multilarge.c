/* multilarge.c
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multilarge.h>

gsl_multilarge_linear_workspace *
gsl_multilarge_linear_alloc(const gsl_multilarge_linear_type *T,
                            const size_t nmax, const size_t p)
{
  gsl_multilarge_linear_workspace *w;

  w = calloc(1, sizeof(gsl_multilarge_linear_workspace));
  if (w == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for workspace",
                     GSL_ENOMEM);
    }

  w->type = T;

  w->state = w->type->alloc(nmax, p);
  if (w->state == NULL)
    {
      gsl_multilarge_linear_free(w);
      GSL_ERROR_NULL("failed to allocate space for multilarge state",
                     GSL_ENOMEM);
    }

  /* initialize newly allocated state */
  gsl_multilarge_linear_reset(w);

  return w;
}

void
gsl_multilarge_linear_free(gsl_multilarge_linear_workspace *w)
{
  RETURN_IF_NULL(w);

  if (w->state)
    w->type->free(w->state);

  free(w);
}

const char *
gsl_multilarge_linear_name(const gsl_multilarge_linear_workspace *w)
{
  return w->type->name;
}

int
gsl_multilarge_linear_reset(gsl_multilarge_linear_workspace *w)
{
  int status = w->type->reset(w->state);
  return status;
}

int
gsl_multilarge_linear_accumulate(const gsl_matrix * X, const gsl_vector * y,
                                 gsl_multilarge_linear_workspace * w)
{
  int status = w->type->accumulate(X, y, w->state);
  return status;
}

int
gsl_multilarge_linear_solve(const double lambda, gsl_vector * c,
                            double * rnorm, double * snorm,
                            gsl_multilarge_linear_workspace * w)
{
  int status = w->type->solve(lambda, c, rnorm, snorm, w->state);
  return status;
}
