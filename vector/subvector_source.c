/* vector/subvector_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
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

QUALIFIED_TYPE(gsl_vector)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
  TYPE(gsl_vector) s = {0, 0, 0, 0};

  if (n == 0)
    {
      GSL_ERROR_VAL ("vector length n must be positive integer", GSL_EINVAL, s);
    }

  if (offset + (n - 1) >= v->size)
    {
      GSL_ERROR_VAL ("vector would extend past end of vector", GSL_EINVAL, s);
    }

  s.data = v->data +  MULTIPLICITY * v->stride * offset ;
  s.size = n;
  s.stride = v->stride;

  return s;
}

QUALIFIED_TYPE(gsl_vector)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n, size_t stride)
{
  TYPE(gsl_vector) s = {0, 0, 0, 0};

  if (n == 0)
    {
      GSL_ERROR_VAL ("vector length n must be positive integer", GSL_EINVAL, s);
    }

  if (stride == 0)
    {
      GSL_ERROR_VAL ("stride must be positive integer", GSL_EINVAL, s);
    }

  if (offset + (n - 1) * stride >= v->size)
    {
      GSL_ERROR_VAL ("vector would extend past end of vector", GSL_EINVAL, s);
    }

  s.data = v->data + MULTIPLICITY * v->stride * offset ;
  s.size = n;
  s.stride = v->stride * stride;

  return s;
}
