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

#include "givens.c"

/* Factorise a general M x N matrix A into
 *
 *   A P = Q R
 *
 * where Q is orthogonal (M x M) and R is upper triangular (M x N).
 *
 * Q is stored as a packed set of Householder transformations in the
 * strict lower triangular part of the input matrix.
 *
 * R is stored in the diagonal and upper triangle of the input matrix.
 *
 * P: column j of P is column k of the identity matrix, where k =
 * permutation->data[j]
 *
 * The full matrix for Q can be obtained as the product
 *
 *       Q = Q_k .. Q_2 Q_1
 *
 * where k = MIN(M,N) and
 *
 *       Q_i = (I - tau_i * v_i * v_i')
 *
 * and where v_i is a Householder vector
 *
 *       v_i = [1, m(i+1,i), m(i+2,i), ... , m(M,i)]
 *
 * This storage scheme is the same as in LAPACK.  */

int
gsl_linalg_QRPT_decomp (gsl_matrix * a, gsl_vector * tau, gsl_permutation * p, int *signum)
{
  const size_t M = a->size1;
  const size_t N = a->size2;

  if (tau->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    }
  else if (p->size != N)
    {
      GSL_ERROR ("permutation size mismatch", GSL_EBADLEN);
    }
  else
    {
      size_t i;
      gsl_vector *work = gsl_vector_alloc (N);
      gsl_vector *norm = gsl_vector_alloc (N);

      *signum = 1;

      gsl_permutation_init (p);	/* set to identity */

      /* Compute column norms and store in workspace */

      for (i = 0; i < N; i++)
	{
	  gsl_vector c = gsl_matrix_column (a, i);
	  double x = gsl_blas_dnrm2 (&c);
	  gsl_vector_set (norm, i, x);
	}

      for (i = 0; i < GSL_MIN (M, N); i++)
	{
	  /* Bring the column of largest norm into the pivot position */

	  double max_norm = 0;
	  size_t j, kmax = i;

	  for (j = i + 1; j < N; j++)
	    {
	      double x = gsl_vector_get (norm, j);

	      if (x < max_norm)
		{
		  max_norm = x;
		  kmax = j;
		}
	    }

	  if (kmax != i)
	    {
	      gsl_matrix_swap_columns (a, i, kmax);
	      gsl_permutation_swap (p, i, kmax);
	      gsl_vector_set (tau, kmax, gsl_vector_get (tau, i));
	      gsl_vector_set (norm, kmax, gsl_vector_get (norm, i));
	      (*signum) = -(*signum);
	    }

	  /* Compute the Householder transformation to reduce the j-th
	     column of the matrix to a multiple of the j-th unit vector */

	  {
	    gsl_vector c_full = gsl_matrix_column (a, i);
            gsl_vector c = gsl_vector_subvector (&c_full, i, M - i);
	    double tau_i = gsl_linalg_householder_transform (&c);

	    gsl_vector_set (tau, i, tau_i);

	    /* Apply the transformation to the remaining columns */

	    if (i + 1 < N)
	      {
		gsl_matrix m = gsl_matrix_submatrix (a, i, i + 1, M - i, N - (i+1));

		gsl_linalg_householder_hm (tau_i, &c, &m, work);
	      }
	  }

	  /* Update the norms of the remaining columns too */

	  for (j = i + 1; j < N; j++)
	    {
	      double y = 0;
	      double x = gsl_vector_get (norm, j);
	      double temp = gsl_matrix_get (a, i, j) / x;

	      if (fabs (temp) >= 1)
		y = 0.0;
	      else
		y = y * sqrt (1 - temp * temp);

	      if (fabs (y / x) < sqrt (20.0) * GSL_SQRT_DBL_EPSILON)
		{
		  gsl_vector c_full = gsl_matrix_column (a, j);
                  gsl_vector c = gsl_vector_subvector(&c_full, i+1, M - (i+1));
		  y = gsl_blas_dnrm2 (&c);
		}

	      gsl_vector_set (norm, j, y);
	    }
	}

      gsl_vector_free (work);
      gsl_vector_free (norm);

      return GSL_SUCCESS;

    }
}

/* Solves the system A x = rhs using the Q R P^T factorisation,

   R z = Q^T rhs

   x = P z;

   to obtain x. Based on SLATEC code. */

int
gsl_linalg_QRPT_solve (const gsl_matrix * qr,
		       const gsl_vector * tau,
		       const gsl_permutation * p,
		       const gsl_vector * rhs,
                       gsl_vector * x)
{
  if (qr->size1 != qr->size2)
    {
      GSL_ERROR ("QR matrix must be square", GSL_ENOTSQR);
    }
  else if (qr->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (qr->size1 != rhs->size)
    {
      GSL_ERROR ("matrix size must match rhs size", GSL_EBADLEN);
    }
  else if (qr->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_memcpy (x, rhs);

      gsl_linalg_QRPT_svx (qr, tau, p, x);
      
      return GSL_SUCCESS;
    }
}

int
gsl_linalg_QRPT_svx (const gsl_matrix * qr,
                     const gsl_vector * tau,
                     const gsl_permutation * p,
                     gsl_vector * x)
{
  if (qr->size1 != qr->size2)
    {
      GSL_ERROR ("QR matrix must be square", GSL_ENOTSQR);
    }
  else if (qr->size1 != p->size)
    {
      GSL_ERROR ("matrix size must match permutation size", GSL_EBADLEN);
    }
  else if (qr->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* compute sol = Q^T b */

      gsl_linalg_QR_QTvec (qr, tau, x);

      /* Solve R x = sol, storing x inplace in sol */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, qr, x);

      gsl_permute_vector_inverse (p, x);

      return GSL_SUCCESS;
    }
}


int
gsl_linalg_QRPT_qrsolve (const gsl_matrix * q, const gsl_matrix * r,
			 const gsl_permutation * p,
			 const gsl_vector * rhs,
			 gsl_vector * x)
{
  if (q->size1 != q->size2 || r->size1 != r->size2)
    {
      return GSL_ENOTSQR;
    }
  else if (q->size1 != p->size || q->size1 != r->size1
	   || q->size1 != rhs->size)
    {
      return GSL_EBADLEN;
    }
  else
    {
      /* compute b = Q^T rhs */

      gsl_blas_dgemv (CblasNoTrans, 1.0, q, rhs, 0.0, x);

      /* Solve R x = b, storing x inplace */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, r, x);

      /* Apply permutation to solution in place */

      gsl_permute_vector_inverse (p, x);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_QRPT_Rsolve (const gsl_matrix * qr,
			const gsl_permutation * p,
			gsl_vector * x)
{
  if (qr->size1 != qr->size2)
    {
      return GSL_ENOTSQR;
    }
  else if (qr->size1 != x->size)
    {
      return GSL_EBADLEN;
    }
  else
    {
      /* Solve R x = b, storing x inplace */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, qr, x);

      gsl_permute_vector_inverse (p, x);

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
gsl_linalg_QRPT_update (gsl_matrix * q, gsl_matrix * r,
			const gsl_permutation * p,
			gsl_vector * w, const gsl_vector * v)
{
  if (q->size1 != q->size2 || r->size1 != r->size2)
    {
      return GSL_ENOTSQR;
    }
  else if (r->size1 != q->size2 || v->size != q->size2 || w->size != q->size2)
    {
      return GSL_EBADLEN;
    }
  else
    {
      size_t j, k;
      const size_t M = q->size1;
      const size_t N = q->size2;
      double w0;

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
	  int p_j = gsl_permutation_get (p, j);
	  double vj = gsl_vector_get (v, p_j);
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
