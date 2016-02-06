/* spio.c
 * 
 * Copyright (C) 2016 Patrick Alken
 * Copyright (C) 2016 Alexis Tantet
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
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>

int
gsl_spmatrix_fprintf(FILE *stream, const gsl_spmatrix *m,
                     const char *format)
{
  int status;

  /* print header */
  status = fprintf(stream, "%zu %zu %zu\n",
                   m->size1, m->size2, m->nz);
  if (status < 0)
    {
      GSL_ERROR("fprintf failed for header", GSL_EFAILED);
    }

  if (GSL_SPMATRIX_ISTRIPLET(m))
    {
      size_t n;

      for (n = 0; n < m->nz; ++n)
        {
          status = fprintf(stream, "%zu %zu ", m->i[n], m->p[n]);
          if (status < 0)
            {
              GSL_ERROR("fprintf failed", GSL_EFAILED);
            }

          status = fprintf(stream, format, m->data[n]);
          if (status < 0)
            {
              GSL_ERROR("fprintf failed", GSL_EFAILED);
            }

          status = putc('\n', stream);
          if (status == EOF)
            {
              GSL_ERROR("putc failed", GSL_EFAILED);
            }
        }
    }
  else
    {
      GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
    }

  return GSL_SUCCESS;
}

int
gsl_spmatrix_fscanf(FILE *stream, gsl_spmatrix *m)
{
  int status;
  size_t size1, size2, nz;

  /* read header */
  status = fscanf (stream, "%zu %zu %zu",
                   &size1, &size2, &nz);
  if (status < 3)
    {
      GSL_ERROR ("fscanf failed reading header", GSL_EFAILED);
    }
  else if (nz > m->nzmax)
    {
      GSL_ERROR ("nzmax not large enough", GSL_EBADLEN);
    }
  else if (size1 != m->size1)
    {
      GSL_ERROR ("matrix has wrong size1", GSL_EBADLEN);
    }
  else if (size2 != m->size2)
    {
      GSL_ERROR ("matrix has wrong size2", GSL_EBADLEN);
    }
  else if (GSL_SPMATRIX_ISTRIPLET(m))
    {
      size_t i, j;
      double val;

      while ((status = fscanf (stream, "%zu %zu %lg", &i, &j, &val)) != EOF)
        {
          if (status < 3)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }

          gsl_spmatrix_set(m, i, j, val);
        }
    }
  else
    {
      GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
    }

  return GSL_SUCCESS;
}
