/* fit/multiple.c
 * 
 * Copyright (C) 2000 Brian Gough
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

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/* Fit

   y = X c

   where X is an M x N matrix of M observations for N variables.

   */

#ifdef JUNK  
int
gsl_fit_linear (gsl_matrix * X,
                gsl_vector * y,
                gsl_vector * c,
                gsl_matrix * cov)
{

  /* Decompose X into U S Q^T */

  gsl_linalg_SV_decomp (X, Q, S);
  
  /* Solve y = X c */

  gsl_linalg_SV_solve (X, Q, S, y, c);

  /* Form covariance matrix cov = Q Q^T */

}
#endif

int
gsl_fit_wmultilinear (gsl_matrix * X,
                 const gsl_vector * w,
                 const gsl_vector * y,
                 gsl_vector * c,
                 gsl_matrix * cov,
                 double * chisq)
{
  if (X->size1 != y->size) 
    {
      GSL_ERROR ("number of observations in y does not match rows of matrix X",
                 GSL_EBADLEN);
    }
  else if (X->size2 != c->size) 
    {
      GSL_ERROR ("number of parameters c does not match columns of matrix X",
                 GSL_EBADLEN);
    }
  else if (w->size != y->size)
    {
      GSL_ERROR ("number of weights does not match number of observations",
                 GSL_EBADLEN);
    }
  else if (cov->size1 != cov->size2)
    {
      GSL_ERROR ("covariance matrix is not square", GSL_ENOTSQR);
    }
  else if (c->size != cov->size1) 
    {
      GSL_ERROR ("number of parameters does not match size of covariance matrix",
                 GSL_EBADLEN);
    }
  else 
    {
      const size_t M = X->size1;
      const size_t N = X->size2;
      
      size_t i, j;

      gsl_matrix * Q = gsl_matrix_alloc (N, N);
      gsl_vector * S = gsl_vector_alloc (N);
      gsl_vector * t = gsl_vector_alloc (M);
      gsl_vector * xt = gsl_vector_alloc (N);
     
      /* Scale X,  A = sqrt(w) X */
      
      for (i = 0; i < M; i++)
        {
          double wi = gsl_vector_get (w, i);
          
          if (wi < 0) wi = 0;
          
          for (j = 0; j < N; j++)
            {
              gsl_vector row = gsl_matrix_row (X, i);
              gsl_vector_scale (&row, sqrt(wi));
            }
        }
      
      /* Decompose A into U S Q^T */
     
      gsl_linalg_SV_decomp (X, Q, S);
      
      /* Solve sqrt(w) y = A c for c, by first computing t = sqrt(w) y */

      for (i = 0; i < M; i++)
        {
          double wi = gsl_vector_get (w, i);
          double yi = gsl_vector_get (y, i);
          if (wi < 0) wi = 0;
          gsl_vector_set (t, i, sqrt(wi) * yi);
        }
      
      gsl_blas_dgemv (CblasTrans, 1.0, X, t, 0.0, xt) ;
      
      /* Scale the matrix Q,  Q' = Q S^-1 */

      for (j = 0 ; j < N ; j++)
        {
          gsl_vector column = gsl_matrix_column (Q, j);
          double alpha = gsl_vector_get (S, j);
          if (alpha != 0) alpha = 1.0 / alpha;
          gsl_vector_scale (&column, alpha);
        }
      
      gsl_blas_dgemv (CblasNoTrans, 1.0, Q, xt, 0.0, c);
      
      /* Form covariance matrix cov = (Q S^-1) (Q S^-1)^T */

      for (i = 0; i < N; i++) 
        {
          gsl_vector row_i = gsl_matrix_row (Q, i);

          for (j = i; j < N ; j++)
            {
              gsl_vector row_j = gsl_matrix_row (Q, j);
              
              double s;
              
              gsl_blas_ddot (&row_i,  &row_j, &s);

              gsl_matrix_set (cov, i, j, s);
              gsl_matrix_set (cov, j, i, s);
            }
        }

      return GSL_SUCCESS;
    }
}
  
  
