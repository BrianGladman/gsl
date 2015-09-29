#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>

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
multifit_linear_svd (const gsl_matrix * X,
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
 *         rnorm    - (output) residual norm ||y - X c||
 *         snorm    - (output) solution norm ||lambda c||
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
                       double *rnorm,
                       double *snorm,
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

      double rho_ls = 0.0;     /* contribution to rnorm from OLS */

      size_t i, j, p_eff;

      gsl_matrix *A = work->A;
      gsl_matrix *Q = work->Q;
      gsl_matrix *QSI = work->QSI;
      gsl_vector *S = work->S;
      gsl_vector *xt = work->xt;
      gsl_vector *D = work->D;
      gsl_vector *t = work->t;

      /*
       * Solve y = A c for c
       * c = Q diag(s_i / (s_i^2 + lambda_i^2)) U^T y
       */

      /* compute xt = U^T y */
      gsl_blas_dgemv (CblasTrans, 1.0, A, y, 0.0, xt);

      if (n > p)
        {
          /*
           * compute OLS residual norm = || y - U U^T y ||;
           * for n = p, U U^T = I, so no need to calculate norm
           */
          gsl_vector_memcpy(t, y);
          gsl_blas_dgemv(CblasNoTrans, -1.0, A, xt, 1.0, t);
          rho_ls = gsl_blas_dnrm2(t);
        }

      if (lambda > 0.0)
        {
          const double lambda_sq = lambda * lambda;

          /* xt <-- [ s(i) / (s(i)^2 + lambda^2) ] .* U^T y */
          for (j = 0; j < p; ++j)
            {
              double sj = gsl_vector_get(S, j);
              double f = (sj * sj) / (sj * sj + lambda_sq);
              double *ptr = gsl_vector_ptr(xt, j);

              /* use D as workspace for residual norm */
              gsl_vector_set(D, j, (1.0 - f) * (*ptr));

              *ptr *= sj / (sj*sj + lambda_sq);
            }

          /* compute regularized solution vector */
          gsl_blas_dgemv (CblasNoTrans, 1.0, Q, xt, 0.0, c);

          /* compute solution norm */
          *snorm = lambda * gsl_blas_dnrm2(c);

          /* compute residual norm */
          *rnorm = gsl_blas_dnrm2(D);

          if (n > p)
            {
              /* add correction to residual norm (see eqs 6-7 of [1]) */
              *rnorm = sqrt((*rnorm) * (*rnorm) + rho_ls * rho_ls);
            }

          /* reset D vector */
          gsl_vector_set_all(D, 1.0);
        }
      else
        {
          /* Scale the matrix Q, QSI = Q S^{-1} */

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

          gsl_blas_dgemv (CblasNoTrans, 1.0, QSI, xt, 0.0, c);

          /* Unscale the balancing factors */
          gsl_vector_div (c, D);

          *snorm = 0.0;
          *rnorm = rho_ls;

          /* variance-covariance matrix cov = s2 * (Q S^-1) (Q S^-1)^T */
          {
            double r2 = (*rnorm) * (*rnorm);
            double s2 = r2 / (n - p_eff);   /* p_eff == rank */

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
        }

      /* compute chisq */
      *chisq = (*rnorm) * (*rnorm) + (*snorm) * (*snorm);

      return GSL_SUCCESS;
    }
}
