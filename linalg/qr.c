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
#include "givens.c"

/* Factorise a general NxN matrix A into 

   A = Q R

   where Q is orthogonal and R is upper triangular.

   Q: diagonal and lower triangle of matrix contains a packed set of
   Householder transformations (see unpack_QR subroutine for unpacking 
   procedure)

   R: strict upper triangle of matrix, with diagonal elements rdiag

   From SLATEC, qrfac.f */

int
gsl_la_decomp_QR_impl (gsl_matrix * matrix, gsl_vector * rdiag)
{
  if (matrix == 0)
    {
      return GSL_EFAULT;
    }
  else if (matrix->size1 != matrix->size2)
    {
      return GSL_ENOTSQR;
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

      for (j = 0; j < N; j++)
	{
	  /* Compute the Householder transformation to reduce the j-th
	     column of the matrix to a multiple of the j-th unit vector */

	  REAL ajnorm = column_norm (matrix, j, M, j);

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

	  gsl_vector_set (rdiag, j, -ajnorm);

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
	    }
	}

      return GSL_SUCCESS;
    }
}

/* Solves the system A x = rhs using the QR factorisation,

   R x = Q^T rhs

   to obtain x. Based on SLATEC code. */

int
gsl_la_solve_QR_impl (const gsl_matrix * qr_matrix,
		      const gsl_vector * rdiag,
		      const gsl_vector * rhs,
		      gsl_vector * solution)
{
  if (qr_matrix == 0 || rhs == 0 || solution == 0)
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
  else if (qr_matrix->size1 != rhs->size
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

      /* Solve R x = sol, storing x inplace in sol */

      gsl_la_Rsolve_QR_impl (qr_matrix, rdiag, solution);

      return GSL_SUCCESS;
    }
}

int
gsl_la_qrsolve_QR_impl (const gsl_matrix * q, const gsl_matrix * r,
                        const gsl_vector * rhs,
                        gsl_vector * solution)
{
  if (q == 0 || r == 0 || rhs == 0 || solution == 0)
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
  else if (q->size1 != r->size1
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

      /* Solve R x = sol, storing x inplace in sol */

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

      return GSL_SUCCESS;
    }
}

int
gsl_la_Rsolve_QR_impl (const gsl_matrix * qr_matrix,
		       const gsl_vector * rdiag,
		       gsl_vector * solution)
{
  if (qr_matrix == 0 || solution == 0)
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
  else if (qr_matrix->size1 == 0)
    {
      return GSL_SUCCESS;	/* FIXME: dumb case */
    }
  else
    {
      const size_t N = qr_matrix->size1;
      size_t n;

      /* Solve R x = sol, storing x inplace in sol */

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

      return GSL_SUCCESS;
    }
}

int
gsl_la_Rsolve_impl (const gsl_matrix * r, gsl_vector * solution)
{
  if (r == 0 || solution == 0)
    {
      return GSL_EFAULT;
    }
  else if (solution->data == 0)
    {
      return GSL_EFAULT;
    }
  else if (r->size1 != r->size2)
    {
      return GSL_ENOTSQR;
    }
  else if (r->size1 == 0)
    {
      return GSL_SUCCESS;	/* FIXME: dumb case */
    }
  else
    {
      const size_t N = r->size1;
      size_t n;

      /* Solve R x = sol, storing x inplace in sol */

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

      return GSL_SUCCESS;
    }
}

/* Form the product Q^T v  from a QR factorized matrix */

int
gsl_la_QTvec_QR_impl (const gsl_matrix * qr_matrix,
                      const gsl_vector * v,
                      gsl_vector * result)
{
  if (qr_matrix == 0 || v == 0 || result == 0)
    {
      return GSL_EFAULT;
    }
  else if (v->data == 0 || result->data == 0)
    {
      return GSL_EFAULT;
    }
  else if (qr_matrix->size1 != qr_matrix->size2)
    {
      return GSL_ENOTSQR;
    }
  else if (qr_matrix->size1 != v->size
	   || qr_matrix->size1 != result->size)
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
      size_t i, j;
      int k;

      /* Copy rhs into solution */


      for (k = 0; k < N; k++)
	{
	  gsl_vector_set (result, k, gsl_vector_get (v, k));
	}

      /* compute Q^T v */

      for (j = 0; j < N; j++)
	{

	  REAL t, sum = 0.0;

	  for (i = j; i < N; i++)
	    {
	      REAL qij = gsl_matrix_get (qr_matrix, i, j);
	      REAL si = gsl_vector_get (result, i);

	      sum += qij * si;
	    }

	  t = -sum / gsl_matrix_get (qr_matrix, j, j);

	  for (i = j; i < N; i++)
	    {
	      REAL qij = gsl_matrix_get (qr_matrix, i, j);
	      REAL si = gsl_vector_get (result, i);

	      gsl_vector_set (result, i, si + t * qij);
	    }
	}
      return GSL_SUCCESS;
    }
}

/*  Form the orthogonal matrix Q from the packed QR matrix */

int
gsl_la_unpack_QR_impl (const gsl_matrix * qr, const gsl_vector * rdiag,
                       gsl_matrix * q, gsl_matrix * r)
{
  if (qr->size1 != qr->size2 || q->size1 != q->size2 || r->size1 != r->size2)
    {
      return GSL_ENOTSQR;
    }
  else if (q->size1 != qr->size2 || q->size1 != r->size1 
           || q->size1 != rdiag->size)
    {
      return GSL_EBADLEN;
    }
  else if (qr->size1 == 0)
    {
      return GSL_SUCCESS;	/* FIXME: what to do with this idiot case ? */
    }
  else
    {
      const int N = qr->size1;
      const int M = qr->size2;
      int i, j, k, L;

      for (j = 0; j < M; j++)
	{
	  for (i = 0; i < j; i++)
	    gsl_matrix_set (q, i, j, 0.0);

	  for (i = j; i < M; i++)
	    gsl_matrix_set (q, i, j, gsl_matrix_get (qr, i, j));
	}

      for (L = 0; L < M; L++)
	{
	  k = (M - 1) - L;

	  for (i = k; i < M; i++)
	    {
	      double ti = gsl_matrix_get (q, i, k); 
              gsl_matrix_set (r, 0, i, ti); /* use as temp storage */
	      gsl_matrix_set (q, i, k, 0.0);
	    }

	  gsl_matrix_set (q, k, k, 1.0);

	  if (gsl_matrix_get(r,0,k) == 0)
	    continue;

	  for (j = k; j < M; j++)
	    {
	      double temp, sum = 0.0;
	      for (i = k; i < M; i++)
		{
                  double ti = gsl_matrix_get(r,0,i);
		  sum += gsl_matrix_get (q, i, j) * ti;
		}
	      temp = sum / gsl_matrix_get(r,0,k);
	      for (i = k; i < M; i++)
		{
		  REAL qij = gsl_matrix_get (q, i, j);
                  double ti = gsl_matrix_get(r,0,i);
		  gsl_matrix_set (q, i, j, qij - temp * ti);
		}
	    }
	}

      /*  Form the right triangular matrix R from a packed QR matrix */

      for (i = 0; i < N; i++)
	{
	  for (j = 0; j < i; j++)
	    gsl_matrix_set (r, i, j, 0.0);

	  gsl_matrix_set (r, i, i, gsl_vector_get (rdiag, i));

	  for (j = i + 1; j < N; j++)
	    gsl_matrix_set (r, i, j, gsl_matrix_get (qr, i, j));
	}

      return GSL_SUCCESS;
    }
}


/* Update a QR factorisation for A= Q R ,  A' = A + u v^T,

   Q' R' = QR + u v^T
         = Q (R + Q^T u v^T)
         = Q (R + w v^T)

   where w = Q^T u.

   Algorithm from Golub and Van Loan, "Matrix Computations", Section
   12.5 (Updating Matrix Factorizations, Rank-One Changes)  */

int
gsl_la_update_QR_impl (gsl_matrix * q, gsl_matrix * r,
		       gsl_vector * w, const gsl_vector * v)
{
  if (q->size1 != q->size2 || r->size1 != r-> size2)
    {
      return GSL_ENOTSQR;
    }
  else if (r->size1 != q->size2  || v->size != q->size2 || w->size != q->size2 )
    {
      return GSL_EBADLEN;
    }
  else if (q->size1 == 0)
    {
      return GSL_SUCCESS;	/* FIXME: what to do with this idiot case ? */
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
        }
      
      return GSL_SUCCESS;
    }
}
