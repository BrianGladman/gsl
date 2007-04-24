/* eigen/subrowcol_source.c
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

static inline
QUALIFIED_VIEW(_gsl_vector,view)
FUNCTION (gsl_matrix, subrow) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i, const size_t offset, const size_t n)
{
  QUALIFIED_VIEW(_gsl_vector,view) view = NULL_VECTOR_VIEW;
  
  if (i >= m->size1)
    {
      GSL_ERROR_VAL ("row index is out of range", GSL_EINVAL, view);
    }
  else if (n == 0)
    {
      GSL_ERROR_VAL ("vector length n must be positive integer",
                     GSL_EINVAL, view);
    }
  else if (offset + n > m->size1)
    {
      GSL_ERROR_VAL ("dimension n overflows matrix", GSL_EINVAL, view);
    }
  
  {
    TYPE(gsl_vector) v = NULL_VECTOR;
    
    v.data = m->data + MULTIPLICITY * (i * m->tda + offset);
    v.size = n;
    v.stride = 1;
    v.block = m->block;
    v.owner = 0;
    
    view.vector = v;
    return view;
  }
}

static inline
QUALIFIED_VIEW(_gsl_vector,view)
FUNCTION (gsl_matrix, subcolumn) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j, const size_t offset, const size_t n)
{
  QUALIFIED_VIEW(_gsl_vector,view) view = NULL_VECTOR_VIEW;
  
  if (j >= m->size2)
    {
      GSL_ERROR_VAL ("column index is out of range", GSL_EINVAL, view);
    }
  else if (n == 0)
    {
      GSL_ERROR_VAL ("vector length n must be positive integer",
                     GSL_EINVAL, view);
    }
  else if (offset + n > m->size2)
    {
      GSL_ERROR_VAL ("dimension n overflows matrix", GSL_EINVAL, view);
    }

  {
    TYPE(gsl_vector) v = NULL_VECTOR;
    
    v.data = m->data + MULTIPLICITY * (offset * m->tda + j);
    v.size = n;
    v.stride = m->tda;
    v.block = m->block;
    v.owner = 0;

    view.vector = v;
    return view;
  }
}
