#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

/* Perform a SVD decomposition on the least squares matrix X = U S Q^T
 *
 * Inputs: X       - least squares matrix
 *                   X may point to work->A in the case of ridge
 *                   regression
 *         balance - 1 to perform column balancing
 *         work    - workspace
 *
 * Notes:
 * 1) On output,
 *    work->A contains the matrix U
 *    work->Q contains the matrix Q
 *    work->S contains the vector of singular values
 */

static int
multifit_linear_svd2 (const gsl_matrix * X,
                     const int balance,
                     gsl_multifit_linear_workspace * work)
{
  if (X->size1 != work->n || X->size2 != work->p)
    {
      GSL_ERROR
        ("size of workspace does not match size of observation matrix",
         GSL_EBADLEN);
    }
  else
    {
      gsl_matrix *A = work->A;
      gsl_matrix *Q = work->Q;
      gsl_matrix *QSI = work->QSI;
      gsl_vector *S = work->S;
      gsl_vector *xt = work->xt;
      gsl_vector *D = work->D;

      /* Copy X to workspace,  A <= X */

      if (X != A)
        gsl_matrix_memcpy (A, X);

      /* Balance the columns of the matrix A if requested */

      if (balance) 
        {
          gsl_linalg_balance_columns (A, D);
        }
      else
        {
          gsl_vector_set_all (D, 1.0);
        }

      /* Decompose A into U S Q^T */

      gsl_linalg_SV_decomp_mod (A, QSI, Q, S, xt);

      return GSL_SUCCESS;
    }
}

/* Fit
 *
 * y = X c
 *
 * where X is an M x N matrix of M observations for N variables.
 *
 * The solution includes a possible standard form Tikhonov regularization:
 *
 * c = (X^T X + lambda^2 I)^{-1} X^T y
 *
 * where lambda^2 is the Tikhonov regularization parameter.
 *
 * The function multifit_linear_svd() must first be called to
 * compute the SVD decomposition of X
 *
 * Inputs: y        - right hand side vector
 *         tol      - singular value tolerance
 *         lambda   - Tikhonov regularization parameter lambda
 *         rank     - (output) effective rank
 *         c        - (output) model coefficient vector
 *         cov      - (output) covariance matrix
 *         rnormsq  - (output) residual norm squared ||X c - y||^2
 *         snormsq  - (output) solution norm squared ||lambda c||^2
 *         chisq    - (output) residual chi^2
 *         work     - workspace
 */

static int
multifit_linear_solve (const gsl_vector * y,
                       const double tol,
                       const double lambda,
                       size_t * rank,
                       gsl_vector * c,
                       gsl_matrix * cov,
                       double *rnormsq,
                       double *snormsq,
                       double *chisq,
                       gsl_multifit_linear_workspace * work)
{
  if (work->n != y->size)
    {
      GSL_ERROR
        ("number of observations in y does not match workspace",
         GSL_EBADLEN);
    }
  else if (work->p != c->size)
    {
      GSL_ERROR ("number of parameters c does not match workspace",
                 GSL_EBADLEN);
    }
  else if (cov->size1 != cov->size2)
    {
      GSL_ERROR ("covariance matrix is not square", GSL_ENOTSQR);
    }
  else if (c->size != cov->size1)
    {
      GSL_ERROR
        ("number of parameters does not match size of covariance matrix",
         GSL_EBADLEN);
    }
  else if (tol <= 0)
    {
      GSL_ERROR ("tolerance must be positive", GSL_EINVAL);
    }
  else if (lambda < 0.0)
    {
      GSL_ERROR ("lambda must be positive", GSL_EINVAL);
    }
  else
    {
      const size_t n = work->n;
      const size_t p = work->p;

      size_t i, j, p_eff;

      gsl_matrix *A = work->A;
      gsl_matrix *Q = work->Q;
      gsl_matrix *QSI = work->QSI;
      gsl_vector *S = work->S;
      gsl_vector *xt = work->xt;
      gsl_vector *D = work->D;

      /*
       * Solve y = A c for c
       * c = Q diag(s_i / (s_i^2 + lambda_i^2)) U^T y
       */

      /* compute xt = U^T y */
      gsl_blas_dgemv (CblasTrans, 1.0, A, y, 0.0, xt);

      if (lambda > 0.0)
        {
          const double lambda_sq = lambda * lambda;

          /* xt <-- [ s(i) / (s(i)^2 + lambda^2) ] .* U^T y */
          for (j = 0; j < p; ++j)
            {
              double sj = gsl_vector_get(S, j);
              double *ptr = gsl_vector_ptr(xt, j);

              *ptr *= sj / (sj*sj + lambda_sq);
            }

          /* compute regularized solution vector */
          gsl_blas_dgemv (CblasNoTrans, 1.0, Q, xt, 0.0, c);
        }
      else
        {
          /* Scale the matrix Q,
           * QSI = Q (S^2 + lambda^2 I)^{-1} S
           *     = Q diag(s_i / (s_i^2 + lambda^2))
           * For standard least squares, lambda = 0 and QSI = Q S^{-1}
           */

          gsl_matrix_memcpy (QSI, Q);

          {
            double s0 = gsl_vector_get (S, 0);
            p_eff = 0;

            for (j = 0; j < p; j++)
              {
                gsl_vector_view column = gsl_matrix_column (QSI, j);
                double sj = gsl_vector_get (S, j);
                double alpha;

                if (sj <= tol * s0)
                  {
                    alpha = 0.0;
                  }
                else
                  {
                    alpha = 1.0 / sj;
                    p_eff++;
                  }

                gsl_vector_scale (&column.vector, alpha);
              }

            *rank = p_eff;
          }

          gsl_vector_set_zero (c);

          gsl_blas_dgemv (CblasNoTrans, 1.0, QSI, xt, 0.0, c);

          /* Unscale the balancing factors */
          gsl_vector_div (c, D);
        }

      /* Compute chisq, from residual r = y - X c */

      {
        double s2 = 0, r2 = 0.0, sn2 = 0.0;

#if 0
        for (i = 0; i < n; i++)
          {
            double yi = gsl_vector_get (y, i);
            gsl_vector_const_view row = gsl_matrix_const_row (X, i);
            double y_est, ri;
            gsl_blas_ddot (&row.vector, c, &y_est);
            ri = yi - y_est;
            r2 += ri * ri;
          }

        /* compute || L c ||^2 contribution to chi^2 */
        for (i = 0; i < p; ++i)
          {
            double ci = gsl_vector_get(c, i);
            sn2 += lambda_sq * ci * ci;
          }

        s2 = r2 / (n - p_eff);   /* p_eff == rank */

        *chisq = r2 + sn2;
        *snormsq = sn2;
        *rnormsq = r2;
#endif

        /* Form variance-covariance matrix cov = s2 * (Q S^-1) (Q S^-1)^T */

        for (i = 0; i < p; i++)
          {
            gsl_vector_view row_i = gsl_matrix_row (QSI, i);
            double d_i = gsl_vector_get (D, i);

            for (j = i; j < p; j++)
              {
                gsl_vector_view row_j = gsl_matrix_row (QSI, j);
                double d_j = gsl_vector_get (D, j);
                double s;

                gsl_blas_ddot (&row_i.vector, &row_j.vector, &s);

                gsl_matrix_set (cov, i, j, s * s2 / (d_i * d_j));
                gsl_matrix_set (cov, j, i, s * s2 / (d_i * d_j));
              }
          }
      }

      return GSL_SUCCESS;
    }
}
