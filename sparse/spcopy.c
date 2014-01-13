/* spcopy.c
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gsl_spmatrix.h"

gsl_spmatrix *
gsl_spmatrix_memcpy(const gsl_spmatrix *src)
{
  gsl_spmatrix *dest;
  size_t n;

  dest = gsl_spmatrix_alloc_nzmax(src->size1, src->size2, src->nz, src->flags);
  if (!dest)
    return NULL;

  /* copy indices and data to dest */
  if (GSLSP_ISTRIPLET(src))
    {
      for (n = 0; n < src->nz; ++n)
        {
          dest->i[n] = src->i[n];
          dest->p[n] = src->p[n];
          dest->data[n] = src->data[n];
        }
    }
  else if (GSLSP_ISCCS(src))
    {
      for (n = 0; n < src->nz; ++n)
        {
          dest->i[n] = src->i[n];
          dest->data[n] = src->data[n];
        }

      for (n = 0; n < src->size2 + 1; ++n)
        {
          dest->p[n] = src->p[n];
        }
    }
  else
    {
      GSL_ERROR_NULL("invalid matrix type for src", GSL_EINVAL);
    }

  dest->nz = src->nz;

  return dest;
} /* gsl_spmatrix_memcpy() */
