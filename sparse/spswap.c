/* spswap.c
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sparse.h>

gsl_spmatrix *
gsl_spmatrix_transpose_memcpy(const gsl_spmatrix *src)
{
  const size_t M = src->size1;
  const size_t N = src->size2;
  const size_t nz = src->nz;
  gsl_spmatrix *dest;

  /* allocate space for transposed matrix */
  dest = gsl_spmatrix_alloc_nzmax(N, M, nz, src->flags);

  if (GSLSP_ISTRIPLET(src))
    {
      size_t n;

      for (n = 0; n < nz; ++n)
        {
          dest->i[n] = src->p[n];
          dest->p[n] = src->i[n];
          dest->data[n] = src->data[n];
        }
    }
  else if (GSLSP_ISCCS(src))
    {
      size_t *Ai = src->i;
      size_t *Ap = src->p;
      double *Ad = src->data;
      size_t *ATi = dest->i;
      size_t *ATp = dest->p;
      double *ATd = dest->data;
      size_t *w = dest->work;
      size_t p, j;

      /* initialize to 0 */
      for (p = 0; p < M + 1; ++p)
        ATp[p] = 0;

      /* compute row counts of A (= column counts for A^T) */
      for (p = 0; p < nz; ++p)
        ATp[Ai[p]]++;

      /* compute row pointers for A (= column pointers for A^T) */
      gsl_spmatrix_cumsum(M, ATp);

      /* make copy of row pointers */
      for (j = 0; j < M; ++j)
        w[j] = ATp[j];

      for (j = 0; j < N; ++j)
        {
          for (p = Ap[j]; p < Ap[j + 1]; ++p)
            {
              size_t k = w[Ai[p]]++;
              ATi[k] = j;
              ATd[k] = Ad[p];
            }
        }
    }
  else
    {
      gsl_spmatrix_free(dest);
      GSL_ERROR_NULL("unknown sparse matrix type", GSL_EINVAL);
    }

  dest->nz = nz;

  return dest;
} /* gsl_spmatrix_transpose_memcpy() */
