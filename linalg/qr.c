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

#include "givens.c"

/* Factorise a general M x N matrix A into
 *  
 *   A = Q R
 *
 * where Q is orthogonal (M x M) and R is upper triangular (M x N).
 *
 * Q is stored as a packed set of Householder transformations in the
 * strict lower triangular part of the input matrix.
 *
 * R is stored in the diagonal and upper triangle of the input matrix.
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
gsl_linalg_QR_decomp (gsl_matrix * a, gsl_vector * tau)
{
  const size_t M = a->size1;
  const size_t N = a->size2;

  if (tau->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      gsl_vector *work = gsl_vector_alloc (N);

      for (i = 0; i < GSL_MIN (M, N); i++)
	{
	  /* Compute the Householder transformation to reduce the j-th
	     column of the matrix to a multiple of the j-th unit vector */

	  gsl_vector c_full = gsl_matrix_column (a, i);
          gsl_vector c = gsl_vector_subvector (&c_full, i, M-i);

	  double tau_i = gsl_linalg_householder_transform (&c);

	  gsl_vector_set (tau, i, tau_i);

	  /* Apply the transformation to the remaining columns and
	     update the norms */

	  if (i + 1 < N)
	    {
	      gsl_matrix m = gsl_matrix_submatrix (a, i, i + 1, M - i, N - (i + 1));

	      gsl_linalg_householder_hm (tau_i, &c, &m, work);
	    }
	}

      gsl_vector_free (work);

      return GSL_SUCCESS;
    }
}

/* Solves the system A x = b using the QR factorisation,

 *  R x = Q^T b
 *
 * to obtain x. Based on SLATEC code. 
 */

int
gsl_linalg_QR_solve (const gsl_matrix * qr, const gsl_vector * tau, const gsl_vector * b, gsl_vector * x)
{
  if (qr->size1 != qr->size2)
    {
      GSL_ERROR ("QR matrix must be square", GSL_ENOTSQR);
    }
  else if (qr->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (qr->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Copy x <- b */

      gsl_vector_memcpy (x, b);

      /* Solve for x */

      gsl_linalg_QR_svx (qr, tau, x);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_QR_svx (const gsl_matrix * qr, const gsl_vector * tau, gsl_vector * x)
{

  if (qr->size1 != qr->size2)
    {
      GSL_ERROR ("QR matrix must be square", GSL_ENOTSQR);
    }
  else if (qr->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match x/rhs size", GSL_EBADLEN);
    }
  else
    {
      /* compute rhs = Q^T b */

      gsl_linalg_QR_QTvec (qr, tau, x);

      /* Solve R x = rhs, storing x in-place */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, qr, x);

      return GSL_SUCCESS;
    }
}


int
gsl_linalgQR_Rsolve (const gsl_matrix * qr, gsl_vector * x)
{
  if (qr->size1 != qr->size2)
    {
      GSL_ERROR ("QR matrix must be square", GSL_ENOTSQR);
    }
  else if (qr->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match rhs size", GSL_EBADLEN);
    }
  else
    {
      /* Solve R x = b, storing x in-place */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, qr, x);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_R_solve (const gsl_matrix * r, const gsl_vector * b, gsl_vector *x)
{
  if (r->size1 != r->size2)
    {
      GSL_ERROR ("R matrix must be square", GSL_ENOTSQR);
    }
  else if (r->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (r->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Copy x <- b */

      gsl_vector_memcpy (x, b);

      /* Solve R x = b, storing x inplace in b */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, r, x);

      return GSL_SUCCESS;
    }
}


/* Form the product Q^T v  from a QR factorized matrix 
 */

int
gsl_linalg_QR_QTvec (const gsl_matrix * qr, const gsl_vector * tau, gsl_vector * v)
{
  const size_t M = qr->size1;
  const size_t N = qr->size2;

  if (tau->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    }
  else if (v->size != M)
    {
      GSL_ERROR ("vector size must be N", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* compute Q^T v */

      for (i = 0; i < GSL_MIN (M, N); i++)
	{
	  gsl_vector c = gsl_matrix_column (qr, i);
          gsl_vector h = gsl_vector_subvector (&c, i, M - i);
	  gsl_vector w = gsl_vector_subvector (v, i, M - i);
	  double ti = gsl_vector_get (tau, i);
	  gsl_linalg_householder_hv (ti, &h, &w);
	}
      return GSL_SUCCESS;
    }
}

/*  Form the orthogonal matrix Q from the packed QR matrix */

int
gsl_linalg_QR_unpack (const gsl_matrix * qr, const gsl_vector * tau, gsl_matrix * q, gsl_matrix * r)
{
  const size_t M = qr->size1;
  const size_t N = qr->size2;

  if (q->size1 != M || q->size2 != M)
    {
      GSL_ERROR ("Q matrix must be M x M", GSL_ENOTSQR);
    }
  else if (r->size1 != M || r->size2 != N)
    {
      GSL_ERROR ("R matrix must be M x N", GSL_ENOTSQR);
    }
  else if (tau->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    }
  else
    {
      size_t i, j, k;

      gsl_vector *work = gsl_vector_alloc (M);

      /* Initialize Q to the identity */

      gsl_matrix_set_identity (q);

      for (k = GSL_MIN (M, N); k > 0; k--)
	{
	  i = k - 1;
	  {
            gsl_vector c = gsl_matrix_column (qr, i);
            gsl_vector h = gsl_vector_subvector (&c, i, M - i);
	    gsl_matrix m = gsl_matrix_submatrix (q, i, i, M - i, M - i);
	    double ti = gsl_vector_get (tau, i);
	    gsl_linalg_householder_hm (ti, &h, &m, work);
	  }
	}

      /*  Form the right triangular matrix R from a packed QR matrix */

      for (i = 0; i < M; i++)
	{
	  for (j = 0; j < i && j < N; j++)
	    gsl_matrix_set (r, i, j, 0.0);

	  for (j = i; j < N; j++)
	    gsl_matrix_set (r, i, j, gsl_matrix_get (qr, i, j));
	}

      gsl_vector_free (work);

      return GSL_SUCCESS;
    }
}


/* Update a QR factorisation for A= Q R ,  A' = A + u v^T,

 * Q' R' = QR + u v^T
 *       = Q (R + Q^T u v^T)
 *       = Q (R + w v^T)
 *
 * where w = Q^T u.
 *
 * Algorithm from Golub and Van Loan, "Matrix Computations", Section
 * 12.5 (Updating Matrix Factorizations, Rank-One Changes)  
 */

int
gsl_linalg_QR_update (gsl_matrix * q, gsl_matrix * r,
		      gsl_vector * w, const gsl_vector * v)
{
  const size_t M = r->size1;
  const size_t N = r->size2;

  if (q->size1 != M || q->size2 != M)
    {
      GSL_ERROR ("Q matrix must be M x M if R is M x N", GSL_ENOTSQR);
    }
  else if (w->size != M)
    {
      GSL_ERROR ("w must be length M if R is M x N", GSL_EBADLEN);
    }
  else if (v->size != N)
    {
      GSL_ERROR ("v must be length N if R is M x N", GSL_EBADLEN);
    }
  else
    {
      size_t j, k;
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
	  double vj = gsl_vector_get (v, j);
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

	  gsl_matrix_set (r, k, k - 1, 0.0);	/* exact zero of G^T */
	}

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_QR_qrsolve (gsl_matrix * q, gsl_matrix * r, gsl_vector * b, gsl_vector * x)
{
  const size_t M = r->size1;
  const size_t N = r->size2;

  if (M != N)
    {
      return GSL_ENOTSQR;
    }
  else if (q->size1 != M || b->size != M || x->size != M)
    {
      return GSL_EBADLEN;
    }
  else
    {
      /* compute sol = Q^T b */

      gsl_blas_dgemv (CblasNoTrans, 1.0, q, b, 0.0, x);

      /* Solve R x = sol, storing x in-place */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, r, x);

      return GSL_SUCCESS;
    }
}
