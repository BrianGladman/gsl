/* Author:  G. Jungman */

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
gsl_linalg_HH_solve (gsl_matrix * a, gsl_vector * x)
{
  if (a->size1 > a->size2)
    {
      /* System is underdetermined. */

      GSL_ERROR ("System is underdetermined", GSL_EINVAL);
    }
  else if (a->size2 != x->size)
    {
      GSL_ERROR ("matrix and vector sizes must be equal", GSL_EBADLEN);
    }
  else
    {
      const int N = a->size1;
      const int M = a->size2;
      int i, j, k;
      REAL *d = (REAL *) malloc (N * sizeof (REAL));

      if (d == 0)
	{
	  GSL_ERROR ("could not allocate memory for workspace", GSL_ENOMEM);
	}

      /* Perform Householder transformation. */

      for (i = 0; i < N; i++)
	{
	  const REAL aii = gsl_matrix_get (a, i, i);
	  REAL alpha;
	  REAL f;
	  REAL ak;
	  REAL max_norm = 0.0;
	  REAL r = 0.0;

	  for (k = i; k < M; k++)
	    {
	      REAL aki = gsl_matrix_get (a, k, i);
	      r += aki * aki;
	    }

	  if (r == 0.0)
	    {
	      /* Rank of matrix is less than size1. */
              free (d);
              GSL_ERROR ("matrix is rank deficient", GSL_ESING);
	    }

	  alpha = sqrt (r) * GSL_SIGN (aii);

	  ak = 1.0 / (r + alpha * aii);
	  gsl_matrix_set (a, i, i, aii + alpha);

	  d[i] = -alpha;

	  for (k = i + 1; k < N; k++)
	    {
	      REAL norm = 0.0;
	      f = 0.0;
	      for (j = i; j < M; j++)
		{
		  REAL ajk = gsl_matrix_get (a, j, k);
		  REAL aji = gsl_matrix_get (a, j, i);
		  norm += ajk * ajk;
		  f += ajk * aji;
		}
	      max_norm = GSL_MAX (max_norm, norm);

	      f *= ak;

	      for (j = i; j < M; j++)
		{
		  REAL ajk = gsl_matrix_get (a, j, k);
		  REAL aji = gsl_matrix_get (a, j, i);
		  gsl_matrix_set (a, j, k, ajk - f * aji);
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
	      f += gsl_vector_get (x, j) * gsl_matrix_get (a, j, i);
	    }
	  f *= ak;
	  for (j = i; j < M; j++)
	    {
	      REAL xj = gsl_vector_get (x, j);
              REAL aji = gsl_matrix_get (a, j, i);
	      gsl_vector_set (x, j, xj - f * aji);
	    }
	}

      /* Perform back-substitution. */

      for (i = N - 1; i >= 0; i--)
	{
	  REAL xi = gsl_vector_get (x, i);
	  REAL sum = 0.0;
	  for (k = i + 1; k < N; k++)
	    {
	      sum += gsl_matrix_get (a, i, k) * gsl_vector_get (x, k);
	    }

	  gsl_vector_set (x, i, (xi - sum) / d[i]);
	}

      free (d);
      return GSL_SUCCESS;
    }
}

