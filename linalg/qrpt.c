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
#include <gsl_math.h>
#include <gsl_vector.h>
#include <gsl_matrix.h>
#include "gsl_linalg.h"

#define REAL double

#include "norm.c"
#include "permut.c"
#include "givens.c"

/* Factorise a general NxN matrix A into 

   A P = Q R

   where P is a permutation matrix, Q is orthogonal and R is upper triangular.

   Q: diagonal and lower triangle of matrix contains a packed set of
   Householder transformations (see unpack_QR subroutine for unpacking 
   procedure)

   R: strict upper triangle of matrix, with diagonal elements rdiag

   P: column j of P is column k of the identity matrix, 
   where k = permutation->data[j]

   Uses pivoting. From SLATEC, qrfac.f */

int
gsl_la_decomp_QRPT_impl (gsl_matrix * matrix,
                         gsl_vector * rdiag,
                         gsl_vector_int * permutation,
                         int *signum)
{
  if (matrix == 0 || permutation == 0 || signum == 0)
    {
      return GSL_EFAULT;
    }
  else if (matrix->size1 != matrix->size2)
    {
      return GSL_ENOTSQR;
    }
  else if (permutation->size != matrix->size2)
    {
      return GSL_EBADLEN;
    }
  else if (rdiag->size != matrix->size2)
    {
      return GSL_EBADLEN;
    }
  else if (matrix->size1 == 0)
    {
      return GSL_SUCCESS;	/* FIXME: what to do with this idiot case ? */
    }
  else
    {
      const int N = matrix->size1;
      const int M = matrix->size2;
      int i, j, k;
      REAL *wa = (REAL *) malloc (N * sizeof (REAL));

      if (wa == 0)
	{
	  return GSL_ENOMEM;
	}

      *signum = 1;

      for (j = 0; j < N; j++)
	{
	  REAL norm = column_norm (matrix, 0, M, j);
	  gsl_vector_set (rdiag, j, norm);
	  wa[j] = norm;
	  gsl_vector_int_set (permutation, j, j);
	}

      for (j = 0; j < N; j++)
	{
	  REAL ajnorm, rk = 0, rkmax = 0;

	  /* Bring the column of largest norm into the pivot position */

	  int kmax = j;

	  for (k = j; k < N; k++)
	    {
	      rk = gsl_vector_get (rdiag, k);

	      if (rk > rkmax)
		{
		  kmax = k;
		  rkmax = rk;
		}
	    }

	  if (kmax != k)
	    {
	      gsl_matrix_swap_cols (matrix, j, kmax);
	      gsl_vector_int_swap (permutation, j, kmax);
	      gsl_vector_set (rdiag, kmax, gsl_vector_get (rdiag, j));
	      wa[kmax] = wa[j];
	      (*signum) = -(*signum);
	    }

	  /* Compute the Householder transformation to reduce the j-th
	     column of the matrix to a multiple of the j-th unit vector */

	  ajnorm = column_norm (matrix, j, M, j);

	  if (ajnorm == 0)
	    {
	      gsl_vector_set (rdiag, j, 0.0);
	      continue;
	    }

	  if (gsl_matrix_get (matrix, j, j) < 0)
	    ajnorm *= -1;

	  for (i = j; i < M; i++)
	    {
	      REAL aij = gsl_matrix_get (matrix, i, j);
	      gsl_matrix_set (matrix, i, j, aij / ajnorm);
	    }

	  gsl_matrix_set (matrix, j, j, 1.0 + gsl_matrix_get (matrix, j, j));

	  /* Apply the transformation to the remaining columns and
	     update the norms */

	  for (k = j + 1; k < N; k++)
	    {
	      REAL temp, sum = 0.0;

	      for (i = j; i < M; i++)
		{
		  REAL aij = gsl_matrix_get (matrix, i, j);
		  REAL aik = gsl_matrix_get (matrix, i, k);
		  sum += aij * aik;
		}

	      temp = sum / gsl_matrix_get (matrix, j, j);

	      for (i = j; i < M; i++)
		{
		  REAL aij = gsl_matrix_get (matrix, i, j);
		  REAL aik = gsl_matrix_get (matrix, i, k);
		  gsl_matrix_set (matrix, i, k, aik - temp * aij);
		}

	      rk = gsl_vector_get (rdiag, k);

	      if (rk == 0)
		continue;

	      temp = gsl_matrix_get (matrix, j, k) / rk;

	      if (fabs (temp) >= 1)
		rk = 0.0;
	      else
		rk = rk * sqrt (1 - temp * temp);

	      if (fabs (rk / wa[k]) < sqrt (20.0) * GSL_SQRT_DBL_EPSILON)
		{
		  rk = column_norm (matrix, j + 1, M, k);
		  wa[k] = rk;
		}

	      gsl_vector_set (rdiag, k, rk);
	    }

	  gsl_vector_set (rdiag, j, -ajnorm);
	}
      free (wa);
      return GSL_SUCCESS;
    }
}

/* Solves the system A x = rhs using the Q R P^T factorisation,

   R z = Q^T rhs

   x = P z;

   to obtain x. Based on SLATEC code. */

int
gsl_la_solve_QRPT_impl (const gsl_matrix * qr_matrix,
                        const gsl_vector * rdiag,
                        const gsl_vector_int * permutation,
                        const gsl_vector * rhs,
                        gsl_vector * solution)
{
  if (qr_matrix == 0 || permutation == 0 || rhs == 0 || solution == 0)
    {
      return GSL_EFAULT;
    }
  else if (solution->data == 0 || rhs->data == 0)
    {
      return GSL_EFAULT;
    }
  else if (qr_matrix->size1 != qr_matrix->size2)
    {
      return GSL_ENOTSQR;
    }
  else if (qr_matrix->size1 != permutation->size
	   || qr_matrix->size1 != rhs->size
	   || qr_matrix->size1 != solution->size)
    {
      return GSL_EBADLEN;
    }
  else if (qr_matrix->size1 == 0)
    {
      return GSL_SUCCESS;	/* FIXME: dumb case */
    }
  else
    {
      /* compute sol = Q^T b */

      gsl_la_QTvec_QR_impl (qr_matrix, rhs, solution);

      /* Solve R z = sol, storing z inplace in sol */

      gsl_la_Rsolve_QRPT_impl (qr_matrix, rdiag, permutation, solution);

      return GSL_SUCCESS;
    }
}

int
gsl_la_qrsolve_QRPT_impl (const gsl_matrix * q, const gsl_matrix * r,
                          const gsl_vector_int * permutation,
                          const gsl_vector * rhs,
                          gsl_vector * solution)
{
  if (q == 0 || r == 0 || permutation == 0 || rhs == 0 || solution == 0)
    {
      return GSL_EFAULT;
    }
  else if (solution->data == 0 || rhs->data == 0)
    {
      return GSL_EFAULT;
    }
  else if (q->size1 != q->size2 || r->size1 != r-> size2)
    {
      return GSL_ENOTSQR;
    }
  else if (q->size1 != permutation->size || q->size1 != r->size1
	   || q->size1 != rhs->size
	   || q->size1 != solution->size)
    {
      return GSL_EBADLEN;
    }
  else if (q->size1 == 0)
    {
      return GSL_SUCCESS;	/* FIXME: dumb case */
    }
  else
    {
      const size_t M = q->size1;
      const size_t N = q->size2;
      size_t i,j,n;

      /* compute sol = Q^T b */

      for (j = 0; j < M; j++)
        {
          double sum = 0;
          for (i = 0; i < N; i++)
            {
              sum += gsl_matrix_get (q, i, j) * gsl_vector_get (rhs, i);
            }
          gsl_vector_set (solution, j, sum);
        }

      /* Solve R z = sol, storing z inplace in sol */

      for (n = N; n > 0; n--)
	{
	  size_t j, i = n - 1;

	  REAL bi, di, sum = 0.0;

	  for (j = n; j < N; j++)
	    {
	      REAL aij = gsl_matrix_get (r, i, j);
	      REAL bj = gsl_vector_get (solution, j);
	      sum += aij * bj;
	    }

	  bi = gsl_vector_get (solution, i);
	  di = gsl_matrix_get (r, i, i);

	  gsl_vector_set (solution, i, (bi - sum) / di);
	}

      /* Apply permutation to solution in place */

      invpermute (permutation, solution);

      return GSL_SUCCESS;
    }
}

int
gsl_la_Rsolve_QRPT_impl (const gsl_matrix * qr_matrix,
                         const gsl_vector * rdiag,
                         const gsl_vector_int * permutation,
                         gsl_vector * solution)
{
  if (qr_matrix == 0 || permutation == 0 || solution == 0)
    {
      return GSL_EFAULT;
    }
  else if (solution->data == 0)
    {
      return GSL_EFAULT;
    }
  else if (qr_matrix->size1 != qr_matrix->size2)
    {
      return GSL_ENOTSQR;
    }
  else if (qr_matrix->size1 != permutation->size
	   || qr_matrix->size1 != solution->size)
    {
      return GSL_EBADLEN;
    }
  else if (qr_matrix->size1 == 0)
    {
      return GSL_SUCCESS;	/* FIXME: dumb case */
    }
  else
    {
      const size_t N = qr_matrix->size1;
      size_t n;

      /* Solve R z = sol, storing z inplace in sol */

      for (n = N; n > 0; n--)
	{
	  size_t j, i = n - 1;

	  REAL bi, di, sum = 0.0;

	  for (j = n; j < N; j++)
	    {
	      REAL aij = gsl_matrix_get (qr_matrix, i, j);
	      REAL bj = gsl_vector_get (solution, j);
	      sum += aij * bj;
	    }

	  bi = gsl_vector_get (solution, i);
	  di = gsl_vector_get (rdiag, i);

	  gsl_vector_set (solution, i, (bi - sum) / di);
	}

      /* Apply permutation to solution in place */

      invpermute (permutation, solution);

      return GSL_SUCCESS;
    }
}


/* Update a Q R P^T factorisation for A P= Q R ,  A' = A + u v^T,

   Q' R' P^-1 = QR P^-1 + u v^T
              = Q (R + Q^T u v^T P ) P^-1
              = Q (R + w v^T P) P^-1

   where w = Q^T u.

   Algorithm from Golub and Van Loan, "Matrix Computations", Section
   12.5 (Updating Matrix Factorizations, Rank-One Changes)  */

int
gsl_la_update_QRPT_impl (gsl_matrix * q, gsl_matrix * r,
                         const gsl_vector_int * permutation,
                         const gsl_vector * u, const gsl_vector * v,
                         gsl_vector * w)
{
  if (q->size1 != q->size2 || r->size1 != r-> size2)
    {
      return GSL_ENOTSQR;
    }
  else if (r->size1 != q->size2 || u->size != q->size2 || v->size != q->size2
           || w->size != q->size2 )
    {
      return GSL_EBADLEN;
    }
  else if (q->size1 == 0)
    {
      return GSL_SUCCESS;	/* FIXME: what to do with this idiot case ? */
    }
  else
    {
      size_t i, j, k;
      const size_t M = q->size1;
      const size_t N = q->size2;
      double w0;
      
      /* First compute w = Q^T u, (Equation 12.5.1) */
      
      for (j = 0; j < M; j++)
        {
          double sum = 0;
          
          for (i = 0; i < N; i++)
            {
              sum += gsl_matrix_get (q, i, j) * gsl_vector_get (u, i);
            }
          
          gsl_vector_set (w, j, sum);
        }
      
      /* Apply Given's rotations to reduce w to (|w|, 0, 0, ... , 0) 
         
         J_1^T .... J_(n-1)^T w = +/- |w| e_1
         
         simultaneously applied to R,  H = J_1^T ... J^T_(n-1) R
         so that H is upper Hessenberg.  (12.5.2) */
      
      for (k = N - 1; k > 0; k--)
        {
          double c, s;
          double wk = gsl_vector_get (w, k);
          double wkm1 = gsl_vector_get (w, k - 1);
          
          create_givens (wkm1, wk, &c, &s);
          apply_givens_vec (w, k - 1, k, c, s);
          apply_givens_qr (M, N, q, r, k - 1, k, c, s);
        }
      
      w0 = gsl_vector_get (w, 0);
      
      /* Add in w v^T  (Equation 12.5.3) */
      
      for (j = 0; j < N; j++)
        {
          double r0j = gsl_matrix_get (r, 0, j);
          int perm_j = gsl_vector_int_get (permutation, j);
          double vj = gsl_vector_get (v, perm_j);
          gsl_matrix_set (r, 0, j, r0j + w0 * vj);
        }
      
      /* Apply Givens transformations R' = G_(n-1)^T ... G_1^T H  
         Equation 12.5.4 */
      
      for (k = 1; k < N; k++)
        {
          double c, s;
          double diag = gsl_matrix_get (r, k - 1, k - 1);
          double offdiag = gsl_matrix_get (r, k, k - 1);
          
          create_givens (diag, offdiag, &c, &s);
          apply_givens_qr (M, N, q, r, k - 1, k, c, s);
        }
      
      return GSL_SUCCESS;
    }
}
