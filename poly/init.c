/* poly/init.c
 * 
 * Copyright (C) 2002 Gert Van den Eynde
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

/*-*-*-*-*-*-*-*-*-*-*-* Allocators *-*-*-*-*-*-*-*-*-*-*-*/

gsl_poly *
gsl_poly_alloc (const size_t size)
{
  gsl_poly *p;

  if (size == 0)
    {
      GSL_ERROR_VAL ("polynomial size must be positive integer",
		     GSL_EINVAL, 0);
    }

  p = (gsl_poly *) malloc (sizeof (gsl_poly));

  if (p == 0)
    {
      GSL_ERROR_VAL ("failed to allocate gsl_poly struct", GSL_ENOMEM, 0);
    }

  p->size = size;

  p->c = (double *) malloc (size * sizeof (double));

  if (p->c == 0)
    {
      free (p);
      GSL_ERROR_VAL ("failed to allocate poly coefficients", GSL_ENOMEM, 0);
    }

  return p;
}

gsl_poly *
gsl_poly_calloc (const size_t size)
{
  size_t k;

  gsl_poly *p = gsl_poly_alloc (size);

  if (p == 0)
    return 0;

  for (k = 0; k < size; k++)
    p->c[k] = 0.0;

  return p;
}

void
gsl_poly_free (gsl_poly * p)
{
  free (p->c);
  free (p);
}
