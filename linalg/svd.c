/* Author:  G. Jungman */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "gsl_linalg.h"

#define REAL double

/* Factorise a general M x N matrix A into,
 *
 *   A = U D Q^T
 *
 * where U is a an M x N matrix (U^T U = I), D is a diagonal N x N
 * matrix and Q is an N x N orthogonal matrix (Q^T Q = I)
 *
 * U is stored in the original matrix A, which has the same size
 *
 * Q is stored as a separate matrix (not Q^T). You must take the
 * transpose to form the product above.
 *
 * The diagonal matrix D is stored in the vector S,  D_ii = S_i
 *
 * The one-sided Jacobi algorithm is used (see Algorithm 4.1 in the
 * following paper).
 *
 * See James Demmel, Kresimir Veselic, "Jacobi's Method is more
 * accurate than QR", Lapack Working Note 15 (LAWN15), October 1989.
 * Available from netlib.  
 */

int
gsl_linalg_SV_decomp (gsl_matrix * A,
                      gsl_matrix * Q,
                      gsl_vector * S,
                      double tolerance)
{
  if (Q->size1 != A->size2)
    {
      GSL_ERROR ("", GSL_EBADLEN);
    }
  else if (Q->size1 != Q->size2)
    {
      GSL_ERROR ("", GSL_ENOTSQR);
    }
  else if (S->size != A->size2)
    {
      GSL_ERROR ("", GSL_EBADLEN);
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
		  double p = 0.0;
		  double q = 0.0;
		  double r = 0.0;
		  double cosine, sine;
		  double v;

		  for (i = 0; i < M; i++)
		    {
		      /* quantities in rotation angles */
		      const REAL Aij = gsl_matrix_get (A, i, j);
		      const REAL Aik = gsl_matrix_get (A, i, k);
		      p += Aij * Aik;
		      q += Aij * Aij;
		      r += Aik * Aik;
		    }

		  /* NOTE: this could be handled better by scaling
		   * the calculation of the inner products above.
		   * But I'm too lazy. This will have to do. [GJ]
		   */
		  if (!(q * r < GSL_DBL_MAX))
		    {
		      /* overflow occured or will occur */
		      return GSL_EOVRFLW;
		    }
		  if (!(q * r > GSL_DBL_MIN))
		    {
		      /* underflow occured or will occur */
		      return GSL_EUNDRFLW;
		    }

		  if (q * r == 0.0)
		    {
		      /* column elements of A are vanishingly small */
		      count--;
		      continue;
		    }

		  if ((double) (p * p) / (double) (q * r) < tolerance)
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
		      v = sqrt (4.0 * p * p + q * q);
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

      if (count > 0)
	{
	  /* reached sweep limit */
	}

      for (j = 0; j < N; j++)
	{
	  double q = 0.0;

	  /* Calculate singular values. */

	  for (i = 0; i < M; i++)
	    {
	      const REAL Aij = gsl_matrix_get (A, i, j);
	      q += Aij * Aij;
	    }
	  gsl_vector_set (S, j, sqrt (q));

	  /* Normalize vectors. */

	  for (i = 0; i < M; i++)
	    {
	      const REAL Aij = gsl_matrix_get (A, i, j);
	      const REAL Sj = gsl_vector_get (S, j);
	      gsl_matrix_set (A, i, j, Aij / Sj);
	    }
	}

      return GSL_SUCCESS;
    }
}
