/* vector/view_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001 Gerard Jungman, Brian Gough
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

int
FUNCTION(gsl_vector, view_from_vector) (TYPE(gsl_vector) * v,
                                        TYPE(gsl_vector) * base,
                                        size_t offset, size_t n, size_t stride)
{
  if (n == 0)
    {
      GSL_ERROR ("vector length n must be positive integer", GSL_EINVAL);
    }

  if (stride == 0)
    {
      GSL_ERROR ("stride must be positive integer", GSL_EINVAL);
    }

  if (base->size <= offset + (n - 1) * stride)
    {
      GSL_ERROR ("vector would extend past end of vector", GSL_EINVAL);
    }

  if (v->block != 0)
    {
      GSL_ERROR ("vector already has memory allocated to it", GSL_ENOMEM);
    }

  v->data = base->data + MULTIPLICITY * base->stride * offset ;
  v->size = n;
  v->stride = base->stride * stride;
  v->block = base->block;
  v->owner = 0;

  return GSL_SUCCESS;
}


int
FUNCTION(gsl_vector, view_from_array) (TYPE(gsl_vector) * v,
                                       ATOMIC * base,
                                       size_t offset, size_t n, size_t stride)
{
  if (n == 0)
    {
      GSL_ERROR ("vector length n must be positive integer", GSL_EINVAL);
    }

  if (stride == 0)
    {
      GSL_ERROR ("stride must be positive integer", GSL_EINVAL);
    }

  v->data = base + offset ;
  v->size = n;
  v->stride = stride;
  v->block = 0;
  v->owner = 0;

  return GSL_SUCCESS;
}


TYPE(gsl_vector)
FUNCTION(gsl_vector, view) (ATOMIC * base, size_t n)
{
  TYPE(gsl_vector) v = {0, 0, 0, 0, 0};

  if (n == 0)
    {
      GSL_ERROR_VAL ("vector length n must be positive integer", 
                     GSL_EINVAL, v);
    }

  v.data = base  ;
  v.size = n;
  v.stride = 1;
  v.block = 0;
  v.owner = 0;

  return v;
}
