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
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "gsl_linalg.h"

#define REAL double

int
gsl_la_decomp_LU_impl (gsl_matrix * matrix, gsl_permutation * permutation, int *signum)
{
  if (matrix->size1 != matrix->size2)
    {
      GSL_ERROR ("LU decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (permutation->size != matrix->size1)
    {
      GSL_ERROR ("permuation length must match matrix size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = matrix->size1;
      size_t i, j, k;
      size_t i_pivot = 0;
      REAL *scale = (REAL *) malloc (N * sizeof (REAL));

      if (scale == 0)
	{
	  GSL_ERROR ("failed to allocate memory for scale", GSL_ENOMEM);
	}

      /* Prepare permutation and scaling information. */

      *signum = 1;

      for (i = 0; i < N; i++)
	{
	  REAL max_row_element = 0.0;

	  for (j = 0; j < N; j++)
	    {
	      REAL aij = fabs (gsl_matrix_get (matrix, i, j));
	      max_row_element = GSL_MAX (max_row_element, aij);
	    }

	  if (max_row_element == 0.0)
	    {
	      /* Trap exact singularity. */
	      *signum = 0;
	      free (scale);
	      GSL_ERROR ("exact singularity", GSL_ESING);
	    }

	  scale[i] = 1.0 / max_row_element;
	}


      /* Perform Crout method. */

      for (j = 0; j < N; j++)
	{
	  REAL max_row_element = 0.0;

	  for (i = 0; i < j; i++)
	    {			/* equation (2.3.12) except for i = j */
	      REAL sum = gsl_matrix_get (matrix, i, j);
	      for (k = 0; k < i; k++)
		{
		  REAL aik = gsl_matrix_get (matrix, i, k);
		  REAL akj = gsl_matrix_get (matrix, k, j);
		  sum -= aik * akj;
		}
	      gsl_matrix_set (matrix, i, j, sum);
	    }

	  for (i = j; i < N; i++)
	    {			/* equation (2.3.13) */
	      REAL dum;
	      REAL sum = gsl_matrix_get (matrix, i, j);
	      for (k = 0; k < j; k++)
		{
		  REAL aik = gsl_matrix_get (matrix, i, k);
		  REAL akj = gsl_matrix_get (matrix, k, j);
		  sum -= aik * akj;
		}
	      gsl_matrix_set (matrix, i, j, sum);

	      dum = scale[i] * fabs (sum);

              /* Is the figure of merit for the pivot the best so far? */

	      if (dum >= max_row_element)
		{
		  max_row_element = dum;
		  i_pivot = i;
		}
	    }

	  /* Perform pivot if non-null. */

	  if (j != i_pivot)
	    {
	      gsl_matrix_swap_rows (matrix, j, i_pivot);
	      *signum = -(*signum);
	      scale[i_pivot] = scale[j];
	    }

	  permutation->data[j] = i_pivot;

	  /* Trap apparent singularity. */

	  if (gsl_matrix_get (matrix, j, j) == 0.0)
	    {
	      *signum = 0;
	      free (scale);
	      GSL_ERROR ("apparent singularity", GSL_ESING);
	    }

	  if (j != N - 1)
	    {
	      REAL ajj = gsl_matrix_get (matrix, j, j);
	      for (i = j + 1; i < N; i++)
		{
		  REAL aij = gsl_matrix_get (matrix, i, j);
		  gsl_matrix_set (matrix, i, j, aij / ajj);
		}
	    }
	}

      free (scale);
      return GSL_SUCCESS;
    }
}

int
gsl_la_solve_LU_impl (const gsl_matrix * lu_matrix, const gsl_permutation * permutation, const gsl_vector * rhs, gsl_vector * solution)
{
  if (lu_matrix->size1 != lu_matrix->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (lu_matrix->size1 != permutation->size)
    {
      GSL_ERROR ("permuation length must match matrix size", GSL_EBADLEN);
    }
  else if (lu_matrix->size1 != rhs->size)
    {
      GSL_ERROR ("matrix size must match rhs size", GSL_EBADLEN);
    }
  else if (lu_matrix->size1 != solution->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = lu_matrix->size1;
      int kk = -1;
      size_t j;
      int k;

      for (k = 0; k < N; k++)
	{
	  gsl_vector_set (solution, k, gsl_vector_get (rhs, k));
	}

      /* Apply permutation to RHS and perform update. */

      for (k = 0; k < N; k++)
	{
	  int perm_index_k = gsl_permutation_get (permutation, k);
	  REAL sum = gsl_vector_get (solution, perm_index_k);
          REAL s_k = gsl_vector_get (solution, k);
	  gsl_vector_set (solution, perm_index_k, s_k);

	  if (kk >= 0)
	    {
	      for (j = kk; j < k; j++)
		{
		  REAL sol_j = gsl_vector_get (solution, j);
		  REAL lum_kj = gsl_matrix_get (lu_matrix, k, j);
		  sum -= lum_kj * sol_j;
		}
	    }
	  else if (sum != 0.0)
	    {
	      kk = k;
	    }
	  gsl_vector_set (solution, k, sum);
	}

      /* Perform back-substitution. */

      for (k = N - 1; k >= 0; k--)
	{
	  REAL sum = gsl_vector_get (solution, k);
	  REAL lum_kk = gsl_matrix_get (lu_matrix, k, k);
	  for (j = k + 1; j < N; j++)
	    {
	      REAL lum_kj = gsl_matrix_get (lu_matrix, k, j);
	      REAL sol_j = gsl_vector_get (solution, j);
	      sum -= lum_kj * sol_j;
	    }

	  if (lum_kk == 0.0)
	    {
	      GSL_ERROR ("singularity in solution", GSL_EINVAL);
	    }
	  else
	    {
	      gsl_vector_set (solution, k, sum / lum_kk);
	    }
	}

      return GSL_SUCCESS;
    }
}


int
gsl_la_invert_LU_impl (const gsl_matrix * lu_matrix,
		       const gsl_permutation * permutation,
		       gsl_matrix * inverse)
{
  size_t i, n = lu_matrix->size1;

  gsl_matrix_set_identity (inverse);

  for (i = 0; i < n; i++)
    {
      int status;

      gsl_vector w = {0, 0, 0, 0};
      gsl_vector_view_col_from_matrix (&w, inverse, i);
      status = gsl_la_solve_LU_impl (lu_matrix, permutation, &w, &w);

      if (status)
	{
	  return status;
	}
    }

  return GSL_SUCCESS;
}

double
gsl_la_det_LU (gsl_matrix * lu_matrix, int signum)
{
  size_t i, n = lu_matrix->size1;

  double det = (double) signum;

  for (i = 0; i < n; i++)
    {
      det *= gsl_matrix_get (lu_matrix, i, i);
    }

  return det;
}


double
gsl_la_lndet_LU (gsl_matrix * lu_matrix)
{
  size_t i, n = lu_matrix->size1;

  double lndet = 0.0;

  for (i = 0; i < n; i++)
    {
      lndet += log (fabs (gsl_matrix_get (lu_matrix, i, i)));
    }

  return lndet;
}


int
gsl_la_sgndet_LU (gsl_matrix * lu_matrix, int signum)
{
  size_t i, n = lu_matrix->size1;

  int s = signum;

  for (i = 0; i < n; i++)
    {
      double u = gsl_matrix_get (lu_matrix, i, i);

      if (u < 0)
	{
	  s *= -1;
	}
      else if (u == 0)
	{
	  s = 0;
	  break;
	}
    }

  return s;
}
