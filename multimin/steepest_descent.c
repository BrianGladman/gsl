/* multimin/steepest_descent.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* steepest_descent.c -- the dumb steepest descent algorithm */

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>

typedef struct
  {
    double dummy;
  }
steepest_descent_state_t;

static int 
steepest_descent_alloc(void *vstate, size_t n)
{
  return GSL_SUCCESS;
}

static int 
steepest_descent_restart(void *vstate)
{
  return GSL_SUCCESS;
}

static void
steepest_descent_free(void *vstate)
{
  /* nothing */
}

static int 
steepest_descent_direction(void *state,gsl_multimin_fdf_history *h ,gsl_vector * dir) 
{
  size_t i;
  
  for(i = 0; i<dir->size; i++) 
    {
      gsl_vector_set(dir,i,-gsl_vector_get(h->g,i));
    }
  return GSL_SUCCESS;
}

static const gsl_multimin_fdf_minimizer_type steepest_descent_type =
{"steepest_descent",			/* name */
 sizeof (steepest_descent_state_t),
 &steepest_descent_alloc,
 &steepest_descent_restart,
 &steepest_descent_direction,
 &steepest_descent_free};

const gsl_multimin_fdf_minimizer_type *gsl_multimin_fdf_minimizer_steepest_descent = &steepest_descent_type;
