/* Author:  G. Jungman */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>

#include "gsl_linalg.h"

#define REAL double

/* Factorise a general N x N matrix A into,
 *
 *   P A = L U
 *
 * where P is a permutation matrix, L is unit lower triangular and R
 * is upper triangular.
 *
 * L is stored in the strict lower triangular part of the input
 * matrix. The diagonal elements of L are unity and are not stored.
 *
 * R is stored in the diagonal and upper triangular part of the
 * input matrix.  
 * 
 * P is stored in the permutation p. Column j of P is column k of the
 * identity matrix, where k = permutation->data[j]
 *
 * signum gives the sign of the permutation, (-1)^n, where n is the
 * number of interchanges in the permutation. 
 *
 * See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss
 * Elimination with Partial Pivoting).
 */

int
gsl_linalg_LU_decomp (gsl_matrix * A, gsl_permutation * p, int *signum)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR ("LU decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (p->size != A->size1)
    {
      GSL_ERROR ("permuation length must match matrix size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = A->size1;
      size_t i, j, k;

      *signum = 1;
      gsl_permutation_init (p);

      for (j = 0; j < N - 1; j++)
	{
	  /* Find maximum in the j-th column */

	  REAL ajj, max = fabs (gsl_matrix_get (A, j, j));
	  size_t i_pivot = j;

	  for (i = j + 1; i < N; i++)
	    {
	      REAL aij = fabs (gsl_matrix_get (A, i, j));

	      if (aij > max)
		{
		  max = aij;
		  i_pivot = i;
		}
	    }

	  if (i_pivot != j)
	    {
	      gsl_matrix_swap_rows (A, j, i_pivot);
	      gsl_permutation_swap (p, j, i_pivot);
	      *signum = -(*signum);
	    }

	  ajj = gsl_matrix_get (A, j, j);

	  if (ajj != 0.0)
	    {
	      for (i = j + 1; i < N; i++)
		{
		  REAL aij = gsl_matrix_get (A, i, j) / ajj;
		  gsl_matrix_set (A, i, j, aij);

		  for (k = j + 1; k < N; k++)
		    {
		      REAL aik = gsl_matrix_get (A, i, k);
		      REAL ajk = gsl_matrix_get (A, j, k);
		      gsl_matrix_set (A, i, k, aik - aij * ajk);
		    }
		}
	    }
	}

#ifdef USE_CROUT
      REAL *scale = (REAL *) malloc (N * sizeof (REAL));

      if (scale == 0)
	{
	  GSL_ERROR ("failed to allocate memory for scale", GSL_ENOMEM);
	}

      /* Prepare permutation and scaling information. */

      *signum = 1;

      gsl_permutation_init (p);

      for (i = 0; i < N; i++)
	{
	  REAL max_row_element = 0.0;

	  for (j = 0; j < N; j++)
	    {
	      REAL aij = fabs (gsl_matrix_get (A, i, j));

	      if (Aij > max_row_element)
		{
		  max_row_element = aij;
		}
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
	      REAL sum = gsl_matrix_get (A, i, j);
	      for (k = 0; k < i; k++)
		{
		  REAL aik = gsl_matrix_get (A, i, k);
		  REAL akj = gsl_matrix_get (A, k, j);
		  sum -= aik * akj;
		}
	      gsl_matrix_set (A, i, j, sum);
	    }

	  for (i = j; i < N; i++)
	    {			/* equation (2.3.13) */
	      REAL dum;
	      REAL sum = gsl_matrix_get (A, i, j);
	      for (k = 0; k < j; k++)
		{
		  REAL aik = gsl_matrix_get (A, i, k);
		  REAL akj = gsl_matrix_get (A, k, j);
		  sum -= aik * akj;
		}
	      gsl_matrix_set (A, i, j, sum);

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
	      gsl_matrix_swap_rows (A, j, i_pivot);
	      gsl_permutation_swap (p, j, i_pivot);
	      *signum = -(*signum);
	      scale[i_pivot] = scale[j];
	    }

	  /* Trap apparent singularity. */

	  if (gsl_matrix_get (A, j, j) == 0.0)
	    {
	      *signum = 0;
	      free (scale);
	      GSL_ERROR ("apparent singularity", GSL_ESING);
	    }

	  if (j != N - 1)
	    {
	      REAL ajj = gsl_matrix_get (A, j, j);
	      for (i = j + 1; i < N; i++)
		{
		  REAL aij = gsl_matrix_get (A, i, j);
		  gsl_matrix_set (A, i, j, aij / ajj);
		}
	    }
	}

      free (scale);
#endif
      return GSL_SUCCESS;
    }
}

int
gsl_linalg_LU_solve (const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)
{
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permuation length must match matrix size", GSL_EBADLEN);
    }
  else if (LU->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LU->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Copy x <- b */

      gsl_vector_memcpy (x, b);

      /* Solve for x */

      gsl_linalg_LU_svx (LU, p, x);

      return GSL_SUCCESS;
    }
}


int
gsl_linalg_LU_svx (const gsl_matrix * LU, const gsl_permutation * p, gsl_vector * x)
{
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permuation length must match matrix size", GSL_EBADLEN);
    }
  else if (LU->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution/rhs size", GSL_EBADLEN);
    }
  else
    {
      /* Apply permutation to RHS */

      gsl_permute_vector (p, x);

      /* Solve for c using forward-substitution, L c = P b */

      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasUnit, LU, x);

      /* Perform back-substitution, U x = c */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LU, x);

      return GSL_SUCCESS;
    }
}


int
gsl_linalg_LU_refine (const gsl_matrix * A, const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x, gsl_vector * residual)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR ("matrix a must be square", GSL_ENOTSQR);
    }
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (A->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be decomposition of a", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permuation length must match matrix size", GSL_EBADLEN);
    }
  else if (LU->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LU->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Compute residual, residual = (A * x  - b) */

      gsl_vector_memcpy (residual, b);
      gsl_blas_dgemv (CblasNoTrans, 1.0, A, x, -1.0, residual);

      /* Find correction, delta = - (A^-1) * residual, and apply it */

      gsl_linalg_LU_svx (LU, p, residual);
      gsl_blas_daxpy (-1.0, residual, x);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_LU_invert (const gsl_matrix * LU, const gsl_permutation * p, gsl_matrix * inverse)
{
  size_t i, n = LU->size1;

  int status = GSL_SUCCESS;

  gsl_matrix_set_identity (inverse);

  for (i = 0; i < n; i++)
    {
      gsl_vector c = gsl_matrix_column (inverse, i);
      int status_i = gsl_linalg_LU_svx (LU, p, &c);

      if (status_i)
	status = status_i;
    }

  return status;
}

double
gsl_linalg_LU_det (gsl_matrix * LU, int signum)
{
  size_t i, n = LU->size1;

  double det = (double) signum;

  for (i = 0; i < n; i++)
    {
      det *= gsl_matrix_get (LU, i, i);
    }

  return det;
}


double
gsl_linalg_LU_lndet (gsl_matrix * LU)
{
  size_t i, n = LU->size1;

  double lndet = 0.0;

  for (i = 0; i < n; i++)
    {
      lndet += log (fabs (gsl_matrix_get (LU, i, i)));
    }

  return lndet;
}


int
gsl_linalg_LU_sgndet (gsl_matrix * LU, int signum)
{
  size_t i, n = LU->size1;

  int s = signum;

  for (i = 0; i < n; i++)
    {
      double u = gsl_matrix_get (LU, i, i);

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
