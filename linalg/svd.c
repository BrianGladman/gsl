/* linalg/svd.c
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

/* Author:  G. Jungman */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "gsl_linalg.h"

#define REAL double

/* Factorise a general M x N matrix A into,
 *
 *   A = U D Q^T
 *
 * where U is a column-orthogonal M x N matrix (U^T U = I), 
 * D is a diagonal N x N matrix, 
 * and Q is an N x N orthogonal matrix (Q^T Q = Q Q^T = I)
 *
 * U is stored in the original matrix A, which has the same size
 *
 * Q is stored as a separate matrix (not Q^T). You must take the
 * transpose to form the product above.
 *
 * The diagonal matrix D is stored in the vector S,  D_ii = S_i
 *
 * Algorithm due to J.C. Nash, Compact Numerical Methods for
 * Computers (New York: Wiley and Sons, 1979), chapter 3.
 * See also Algorithm 4.1 in
 * James Demmel, Kresimir Veselic, "Jacobi's Method is more
 * accurate than QR", Lapack Working Note 15 (LAWN15), October 1989.
 * Available from netlib.
 *
 * Based on code by Arthur Kosowsky, Rutgers University
 *                  kosowsky@physics.rutgers.edu  
 */

int
gsl_linalg_SV_decomp (gsl_matrix * A,
                      gsl_matrix * Q,
                      gsl_vector * S)
{
  if (Q->size1 < Q->size2)
    {
      /* FIXME: only implemented  M>=N case so far */

      GSL_ERROR ("svd of MxN matrix, M<N, is not implemented", GSL_EUNIMPL);
    }
  else if (Q->size1 != A->size2)
    {
      GSL_ERROR ("square matrix Q must match second dimension of matrix A", GSL_EBADLEN);
    }
  else if (Q->size1 != Q->size2)
    {
      GSL_ERROR ("matrix Q must be square", GSL_ENOTSQR);
    }
  else if (S->size != A->size2)
    {
      GSL_ERROR ("length of vector S must match second dimension of matrix A", GSL_EBADLEN);
    }
  else
    {
      const int M = A->size1;
      const int N = A->size2;
      int i, j, k;

      /* Initialize the rotation counter and the sweep counter. */
      int count = 1;
      int sweep = 0;
      int sweepmax = N;

      double tolerance = 10 * GSL_DBL_EPSILON ;

      /* Always do at least 12 sweeps. */
      sweepmax = GSL_MAX (sweepmax, 12);

      /* Set Q to the identity matrix. */
      gsl_matrix_set_identity (Q);

      /* Orthogonalize A by plane rotations. */

      while (count > 0 && sweep <= sweepmax)
	{
	  /* Initialize rotation counter. */
	  count = N * (N - 1) / 2;

	  for (j = 0; j < N - 1; j++)
	    {
	      for (k = j + 1; k < N; k++)
		{
                  double a = 0.0;
                  double b = 0.0;
		  double p = 0.0;
		  double q = 0.0;
		  double r = 0.0;
		  double cosine, sine;
		  double v;

                  gsl_vector cj = gsl_matrix_column (A, j);
                  gsl_vector ck = gsl_matrix_column (A, k);

                  gsl_blas_ddot (&cj, &ck, &p);
                  
                  a = gsl_blas_dnrm2 (&cj);
                  b = gsl_blas_dnrm2 (&ck);

                  q = a * a ;
                  r = b * b ;

		  /* NOTE: this could be handled better by scaling
		   * the calculation of the inner products above.
		   * But I'm too lazy. This will have to do. [GJ]
		   */

                  /* FIXME: This routine is a hack. We need to get the
                     state of the art in Jacobi SVD's here ! */
                  
                  /* This is an adhoc method of testing for a "zero"
                     singular value. We consider it to be zero if it
                     is sufficiently small compared with the currently
                     leading column. Note that b <= a is guaranteed by
                     the sweeping algorithm. BJG */
                  
                  if (b <= tolerance * a) 
                    {
                      /* probably |b| = 0 */
                      count--;
                      continue;
                    }

		  if (fabs(p) <= tolerance * a * b)
		    {
		      /* columns j,k orthogonal
		       * note that p*p/(q*r) is automatically <= 1.0
		       */
		      count--;
		      continue;
		    }

		  /* calculate rotation angles */
		  if (q < r)
		    {
		      cosine = 0.0;
		      sine = 1.0;
		    }
		  else
		    {
		      q -= r;
		      v = hypot (2.0 * p, q);
		      cosine = sqrt ((v + q) / (2.0 * v));
		      sine = p / (v * cosine);
		    }

		  /* apply rotation to A */
		  for (i = 0; i < M; i++)
		    {
		      const REAL Aik = gsl_matrix_get (A, i, k);
		      const REAL Aij = gsl_matrix_get (A, i, j);
		      gsl_matrix_set (A, i, j, Aij * cosine + Aik * sine);
		      gsl_matrix_set (A, i, k, -Aij * sine + Aik * cosine);
		    }

		  /* apply rotation to Q */
		  for (i = 0; i < N; i++)
		    {
		      const REAL Qij = gsl_matrix_get (Q, i, j);
		      const REAL Qik = gsl_matrix_get (Q, i, k);
		      gsl_matrix_set (Q, i, j, Qij * cosine + Qik * sine);
		      gsl_matrix_set (Q, i, k, -Qij * sine + Qik * cosine);
		    }
		}
	    }

	  /* Sweep completed. */
	  sweep++;
	}

      /* 
       * Orthogonalization complete. Compute singular values.
       */

      {
        double prev_norm = -1.0;

        for (j = 0; j < N; j++)
          {
            gsl_vector column = gsl_matrix_column (A, j);
            double norm =  gsl_blas_dnrm2 (&column);
            
            /* Determine if singular value is zero, according to the
               criteria used in the main loop above (i.e. comparison
               with norm of previous column). */

            if (norm == 0 || (j > 0 && norm <= tolerance * prev_norm))
              {
                gsl_vector_set (S, j, 0.0);       /* singular */
                gsl_vector_set_zero (&column);    /* annihilate column */

                prev_norm = 0;        
              }
            else
              {
                gsl_vector_set (S, j, norm);            /* non-singular */
                gsl_vector_scale (&column, 1.0/norm);   /* normalize column */

                prev_norm = norm;   
              }
          }
      }

      if (count > 0)
	{
	  /* reached sweep limit */
          GSL_ERROR ("Jacobi iterations did not reach desired tolerance", GSL_ETOL);
	}

      return GSL_SUCCESS;
    }
}

/*  Solves the system A x = b using the SVD factorization
 *
 *  A = U S Q^T
 *
 *  to obtain x. For M x N systems it finds the solution in the least
 *  squares sense.  
 */

int
gsl_linalg_SV_solve (const gsl_matrix * U,
                     const gsl_matrix * Q,
                     const gsl_vector * S,
                     const gsl_vector * b,
                     gsl_vector * x)
{
  if (U->size1 != b->size)
    {
      GSL_ERROR ("first dimension of matrix Q must size of vector b", GSL_EBADLEN);
    }
  else if (U->size2 != S->size)
    {
      GSL_ERROR ("length of vector S must match second dimension of matrix U", GSL_EBADLEN);
    }
  else if (Q->size1 != Q->size2)
    {
      GSL_ERROR ("matrix Q must be square", GSL_ENOTSQR);
    }
  else if (S->size != Q->size1)
    {
      GSL_ERROR ("length of vector S must match size of matrix Q", GSL_EBADLEN);
    }
  else if (Q->size2 != x->size)
    {
      GSL_ERROR ("size of matrix Q must match size of vector x", GSL_EBADLEN);
    }
  else
    {
      const size_t N = U->size2;
      size_t i;

      gsl_vector * w = gsl_vector_calloc (N);
      
      gsl_blas_dgemv (CblasTrans, 1.0, U, b, 0.0, w) ;
      
      for (i = 0 ; i < N ; i++)
        {
          double wi = gsl_vector_get (w, i);
          double alpha = gsl_vector_get (S, i);
          if (alpha != 0) alpha = 1.0 / alpha;
          gsl_vector_set (w, i, alpha * wi);
        }

      gsl_blas_dgemv (CblasNoTrans, 1.0, Q, w, 0.0, x);
        
      gsl_vector_free (w);

      return GSL_SUCCESS;
    }
}
