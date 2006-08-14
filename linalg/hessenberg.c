/* linalg/hessenberg.c
 * 
 * Copyright (C) 2006 Patrick Alken
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <gsl/gsl_linalg.h>


/* Compute the Householder reduction to Hessenberg form of a
 * square matrix A.
 * 
 * See Golub & Van Loan, "Matrix Computations" (3rd ed), algorithm
 * 7.4.2
 */

int
gsl_linalg_hessenberg (gsl_matrix * A, gsl_vector * w)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR ("Hessenberg reduction requires square matrix", GSL_ENOTSQR);
    }
  else
    {
      const size_t N = A->size1;
      size_t i, j;
      gsl_vector_view v;
      gsl_matrix_view m;
      double tau_i;             /* beta in algorithm 7.4.2 */

      for (i = 0; i < N - 2; ++i)
        {
          /* make a copy of A(i + 1:n, i) */
          {
            v = gsl_matrix_column (A, i);
            gsl_vector_memcpy (w, &v.vector);
          }
          v = gsl_vector_subvector (w, i + 1, N - (i + 1));

          /* compute householder transformation of A(i+1:n,i) */
          tau_i = gsl_linalg_householder_transform (&v.vector);

          /* apply left householder matrix (I - tau_i v v') to A */
          m = gsl_matrix_submatrix (A, i + 1, i, N - (i + 1), N - i);
          gsl_linalg_householder_hm (tau_i, &v.vector, &m.matrix);

          /* apply right householder matrix (I - tau_i v v') to A */
          m = gsl_matrix_submatrix (A, 0, i + 1, N, N - (i + 1));
          gsl_linalg_householder_mh (tau_i, &v.vector, &m.matrix);
        }

      return GSL_SUCCESS;
    }
}
