/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Simple linear algebra operations, operating directly
 * on the gsl_vector and gsl_matrix objects. These are
 * meant for "generic" and "small" systems. Anyone
 * interested in large systems will want to use more
 * sophisticated methods, presumably involving native
 * BLAS operations, specialized data representations,
 * or other optimizations.
 */
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "gsl_linalg.h"

#define REAL double

/* [Engeln-Mullges + Uhlig, Alg. 4.42]
 */
int
gsl_la_solve_HH (gsl_matrix * matrix, gsl_vector * vec)
{
  if (matrix->size1 > matrix->size2)
    {
      /* System is underdetermined. */

      GSL_ERROR ("System is underdetermined", GSL_EINVAL);
    }
  else if (matrix->size2 != vec->size)
    {
      GSL_ERROR ("matrix and vector sizes must be equal", GSL_EBADLEN);
    }
  else
    {
      const int N = matrix->size1;
      const int M = matrix->size2;
      int i, j, k;
      REAL *d = (REAL *) malloc (N * sizeof (REAL));

      if (d == 0)
	{
	  GSL_ERROR ("could not allocate memory for workspace", GSL_ENOMEM);
	}

      /* Perform Householder transformation. */

      for (i = 0; i < N; i++)
	{
	  const REAL elem_ii = gsl_matrix_get (matrix, i, i);
	  REAL alpha;
	  REAL f;
	  REAL ak;
	  REAL max_norm = 0.0;
	  REAL r = 0.0;

	  for (k = i; k < M; k++)
	    {
	      REAL elem_ki = gsl_matrix_get (matrix, k, i);
	      r += elem_ki * elem_ki;
	    }

	  if (r == 0.0)
	    {
	      /* Rank of matrix is less than size1. */
              free (d);
              GSL_ERROR ("matrix is rank deficient", GSL_ESING);
	    }

	  alpha = sqrt (r) * GSL_SIGN (elem_ii);

	  ak = 1.0 / (r + alpha * elem_ii);
	  gsl_matrix_set (matrix, i, i, elem_ii + alpha);

	  d[i] = -alpha;

	  for (k = i + 1; k < N; k++)
	    {
	      REAL norm = 0.0;
	      f = 0.0;
	      for (j = i; j < M; j++)
		{
		  REAL elem_jk = gsl_matrix_get (matrix, j, k);
		  REAL elem_ji = gsl_matrix_get (matrix, j, i);
		  norm += elem_jk * elem_jk;
		  f += elem_jk * elem_ji;
		}
	      max_norm = GSL_MAX (max_norm, norm);

	      f *= ak;

	      for (j = i; j < M; j++)
		{
		  REAL elem_jk = gsl_matrix_get (matrix, j, k);
		  REAL elem_ji = gsl_matrix_get (matrix, j, i);
		  gsl_matrix_set (matrix, j, k, elem_jk - f * elem_ji);
		}
	    }

	  if (fabs (alpha) < 2.0 * GSL_DBL_EPSILON * sqrt (max_norm))
	    {
	      /* Apparent singularity. */
	      free (d);
	      GSL_ERROR("apparent singularity detected", GSL_ESING);
	    }

	  /* Perform update of RHS. */

	  f = 0.0;
	  for (j = i; j < M; j++)
	    {
	      f += gsl_vector_get (vec, j) * gsl_matrix_get (matrix, j, i);
	    }
	  f *= ak;
	  for (j = i; j < M; j++)
	    {
	      REAL vec_j = gsl_vector_get (vec, j);
              REAL mat_ji = gsl_matrix_get (matrix, j, i);
	      gsl_vector_set (vec, j, vec_j - f * mat_ji);
	    }
	}

      /* Perform back-substitution. */

      for (i = N - 1; i >= 0; i--)
	{
	  REAL vec_i = gsl_vector_get (vec, i);
	  REAL sum = 0.0;
	  for (k = i + 1; k < N; k++)
	    {
	      sum += gsl_matrix_get (matrix, i, k) * gsl_vector_get (vec, k);
	    }

	  gsl_vector_set (vec, i, (vec_i - sum) / d[i]);
	}

      free (d);
      return GSL_SUCCESS;
    }
}

