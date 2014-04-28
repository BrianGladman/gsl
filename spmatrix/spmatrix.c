/* spmatrix.c
 * 
 * Copyright (C) 2012 Patrick Alken
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

#include "avl.c"

static int compare_triplet(const void *pa, const void *pb, void *param);

/*
gsl_spmatrix_alloc()
  Allocate a sparse matrix in triplet representation

Inputs: n1 - number of rows
        n2 - number of columns

Notes: if (n1,n2) are not known at allocation time, they can each be
set to 1, and they will be expanded as elements are added to the matrix
*/

gsl_spmatrix *
gsl_spmatrix_alloc(const size_t n1, const size_t n2)
{
  const double density = 0.1; /* estimate */
  size_t nzmax = (size_t) floor(n1 * n2 * density);

  if (nzmax == 0)
    nzmax = 10;

  return gsl_spmatrix_alloc_nzmax(n1, n2, nzmax, GSL_SPMATRIX_TRIPLET);
} /* gsl_spmatrix_alloc() */

/*
gsl_spmatrix_alloc_nzmax()
  Allocate a sparse matrix with given nzmax

Inputs: n1     - number of rows
        n2     - number of columns
        nzmax  - maximum number of matrix elements
        sptype - type of matrix (triplet, compressed column)

Notes: if (n1,n2) are not known at allocation time, they can each be
set to 1, and they will be expanded as elements are added to the matrix
*/

gsl_spmatrix *
gsl_spmatrix_alloc_nzmax(const size_t n1, const size_t n2,
                         const size_t nzmax, const size_t sptype)
{
  gsl_spmatrix *m;

  if (n1 == 0)
    {
      GSL_ERROR_VAL ("matrix dimension n1 must be positive integer",
                     GSL_EINVAL, 0);
    }
  else if (n2 == 0)
    {
      GSL_ERROR_VAL ("matrix dimension n2 must be positive integer",
                     GSL_EINVAL, 0);
    }

  m = calloc(1, sizeof(gsl_spmatrix));
  if (!m)
    {
      GSL_ERROR_VAL("failed to allocate space for spmatrix struct",
                    GSL_ENOMEM, 0);
    }

  m->size1 = n1;
  m->size2 = n2;
  m->nz = 0;
  m->nzmax = GSL_MAX(nzmax, 1);
  m->sptype = sptype;

  m->i = malloc(m->nzmax * sizeof(size_t));
  if (!m->i)
    {
      gsl_spmatrix_free(m);
      GSL_ERROR_VAL("failed to allocate space for row indices",
                    GSL_ENOMEM, 0);
    }

  if (sptype == GSL_SPMATRIX_TRIPLET)
    {
      m->p = malloc(m->nzmax * sizeof(size_t));
      if (!m->p)
        {
          gsl_spmatrix_free(m);
          GSL_ERROR_VAL("failed to allocate space for column indices",
                        GSL_ENOMEM, 0);
        }

      m->btree = avl_create(compare_triplet, (void *) m, NULL);
      if (!m->btree)
        {
          gsl_spmatrix_free(m);
          GSL_ERROR_VAL("failed to allocate space for AVL tree",
                        GSL_ENOMEM, 0);
        }
    }
  else if (sptype == GSL_SPMATRIX_CCS)
    {
      m->p = malloc((n2 + 1) * sizeof(size_t));
      m->work = malloc(GSL_MAX(n1, n2) *
                       GSL_MAX(sizeof(size_t), sizeof(double)));
      if (!m->p || !m->work)
        {
          gsl_spmatrix_free(m);
          GSL_ERROR_VAL("failed to allocate space for column pointers",
                        GSL_ENOMEM, 0);
        }
    }

  m->data = malloc(m->nzmax * sizeof(double));
  if (!m->data)
    {
      gsl_spmatrix_free(m);
      GSL_ERROR_VAL("failed to allocate space for data",
                    GSL_ENOMEM, 0);
    }

  return m;
} /* gsl_spmatrix_alloc_nzmax() */

/*
gsl_spmatrix_free()
  Free sparse matrix object
*/

void
gsl_spmatrix_free(gsl_spmatrix *m)
{
  if (m->i)
    free(m->i);

  if (m->p)
    free(m->p);

  if (m->data)
    free(m->data);

  if (m->work)
    free(m->work);

  if (m->btree)
    avl_destroy(m->btree, NULL);

  free(m);
} /* gsl_spmatrix_free() */

/*
gsl_spmatrix_realloc()
  As elements are added to the sparse matrix, its possible that they
will exceed the previously specified nzmax - reallocate the matrix
with a new nzmax
*/

int
gsl_spmatrix_realloc(const size_t nzmax, gsl_spmatrix *m)
{
  int s = GSL_SUCCESS;
  void *ptr;

  if (nzmax < m->nz)
    {
      GSL_ERROR("new nzmax is less than current nz", GSL_EINVAL);
    }

  ptr = realloc(m->i, nzmax * sizeof(size_t));
  if (!ptr)
    {
      GSL_ERROR("failed to allocate space for row indices", GSL_ENOMEM);
    }

  m->i = (size_t *) ptr;

  if (GSL_SPMATRIX_ISTRIPLET(m))
    {
      ptr = realloc(m->p, nzmax * sizeof(size_t));
      if (!ptr)
        {
          GSL_ERROR("failed to allocate space for column indices", GSL_ENOMEM);
        }

      m->p = (size_t *) ptr;
    }

  ptr = realloc(m->data, nzmax * sizeof(double));
  if (!ptr)
    {
      GSL_ERROR("failed to allocate space for data", GSL_ENOMEM);
    }

  m->data = (double *) ptr;

  m->nzmax = nzmax;

  /* rebuild binary tree */
  if (GSL_SPMATRIX_ISTRIPLET(m))
    {
      size_t n;

      avl_destroy(m->btree, NULL);
      m->btree = avl_create(compare_triplet, (void *) m, NULL);
      if (!m->btree)
        {
          gsl_spmatrix_free(m);
          GSL_ERROR("failed to allocate space for AVL tree", GSL_ENOMEM);
        }

      for (n = 0; n < m->nz; ++n)
        {
          void *ptr = avl_insert(m->btree, &m->data[n]);
          if (ptr != NULL)
            {
              GSL_ERROR("detected duplicate entry", GSL_EINVAL);
            }
        }
    }

  return s;
} /* gsl_spmatrix_realloc() */

int
gsl_spmatrix_set_zero(gsl_spmatrix *m)
{
  m->nz = 0;
  return GSL_SUCCESS;
} /* gsl_spmatrix_set_zero() */

size_t
gsl_spmatrix_nnz(const gsl_spmatrix *m)
{
  return m->nz;
} /* gsl_spmatrix_nnz() */

int
gsl_spmatrix_fprintf(FILE *stream, const gsl_spmatrix *m,
                     const char *format)
{
  int s = 0;

  return s;
} /* gsl_spmatrix_fprintf() */

/*
gsl_spmatrix_compare_idx()
  Comparison function for searching binary tree in triplet
representation.

To detect duplicate elements in the tree, we want to determine
if there already exists an entry for (i,j) in the tree. Since
the actual tree node stores only the data elements data[n],
we will do pointer arithmetic to get from the given data[n]
to the row/column indices i[n] and j[n].

This compare function will sort the tree first by row i,
and for equal rows, it will then sort by column j

Inputs: ia - row index of element a
        ja - column index of element a
        ib - row index of element b
        jb - column index of element b

Return:
  -1 if pa < pb: (ia,ja) < (ib,jb)
  +1 if pa > pb: (ia,ja) > (ib,jb)
   0 if pa = pb: (ia,ja) == (ib,jb)
*/

int
gsl_spmatrix_compare_idx(const size_t ia, const size_t ja,
                         const size_t ib, const size_t jb)
{
  if (ia < ib)
    return -1;
  else if (ia > ib)
    return 1;
  else
    {
      /* row indices are equal, sort by column index */
      if (ja < jb)
        return -1;
      else if (ja > jb)
        return 1;
      else
        return 0; /* row and column indices are equal */
    }
}

/*
compare_triplet()
  Comparison function for searching binary tree in triplet
representation.

To detect duplicate elements in the tree, we want to determine
if there already exists an entry for (i,j) in the tree. Since
the actual tree node stores only the data elements data[n],
we will do pointer arithmetic to get from the given data[n]
to the row/column indices i[n] and j[n].

This compare function will sort the tree first by row i,
and for equal rows, it will then sort by column j

Inputs: pa    - element 1 for comparison (double *) 
        pb    - element 2 for comparison (double *)
        param - parameter (gsl_spmatrix)

Return:
  -1 if pa < pb: (ia,ja) < (ib,jb)
  +1 if pa > pb: (ia,ja) > (ib,jb)
   0 if pa = pb: (ia,ja) == (ib,jb)
*/

static int
compare_triplet(const void *pa, const void *pb, void *param)
{
  gsl_spmatrix *m = (gsl_spmatrix *) param;

  /* pointer arithmetic to find indices in data array */
  const size_t idxa = (double *) pa - m->data;
  const size_t idxb = (double *) pb - m->data;

  return gsl_spmatrix_compare_idx(m->i[idxa], m->p[idxa], m->i[idxb], m->p[idxb]);
} /* compare_triplet() */
