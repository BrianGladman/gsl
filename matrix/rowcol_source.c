/* matrix/rowcol_source.c
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

QUALIFIED_TYPE(gsl_matrix)
FUNCTION (gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * matrix, 
                                  const size_t i, const size_t j,
                                  const size_t m, const size_t n)
{
  TYPE(gsl_matrix) s = {0, 0, 0, 0};

  if (i >= matrix->size1)
    {
      GSL_ERROR_VAL ("row index is out of range", GSL_EINVAL, s);
    }
  else if (j >= matrix->size2)
    {
      GSL_ERROR_VAL ("column index is out of range", GSL_EINVAL, s);
    }
  else if (m == 0)
    {
      GSL_ERROR_VAL ("first dimension must be non-zero", GSL_EINVAL, s);
    }
  else if (n == 0)
    {
      GSL_ERROR_VAL ("second dimension must be non-zero", GSL_EINVAL, s);
    }
  else if (i + m > matrix->size1)
    {
      GSL_ERROR_VAL ("first dimension overflows matrix", GSL_EINVAL, s);
    }
  else if (j + n > matrix->size2)
    {
      GSL_ERROR_VAL ("second dimension overflows matrix", GSL_EINVAL, s);
    }

  s.data = matrix->data + i * matrix->tda + j;
  s.size1 = m;
  s.size2 = n;
  s.tda = matrix->tda;

  return s;
}

QUALIFIED_TYPE(gsl_vector)
FUNCTION (gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
  TYPE(gsl_vector) v = {0, 0, 0, 0};

  if (i >= m->size1)
    {
      GSL_ERROR_VAL ("row index is out of range", GSL_EINVAL, v);
    }

  v.data = m->data + i * m->tda;
  v.size = m->size2;
  v.stride = 1;

  return v;
}

QUALIFIED_TYPE(gsl_vector)
FUNCTION (gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
  TYPE(gsl_vector) v = {0, 0, 0, 0};

  if (j >= m->size2)
    {
      GSL_ERROR_VAL ("column index is out of range", GSL_EINVAL, v);
    }

  v.data = m->data + j;
  v.size = m->size1;
  v.stride = m->tda;

  return v;
}

QUALIFIED_TYPE(gsl_vector)
FUNCTION (gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
  TYPE(gsl_vector) v = {0, 0, 0, 0};

  v.data = m->data;
  v.size = GSL_MIN(m->size1,m->size2);
  v.stride = m->tda + 1;

  return v;
}
