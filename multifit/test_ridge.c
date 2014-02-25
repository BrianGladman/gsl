/* test linear ridge regression */
static void
test_ridge(void)
{
  const size_t n = 100;
  const size_t p = 10;
  const double xmin = -1.0;
  const double xmax = 1.0;
  const double dx = (xmax - xmin) / (n - 1.0);
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  double *x = malloc(n * sizeof(double));
  double *y = malloc(n * sizeof(double));
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  size_t i, j;

  /* construct artificial data */
  for (i = 0; i < n; ++i)
    {
      double ei = 0.2 * gsl_rng_uniform(r);

      x[i] = xmin + dx * i;
      y[i] = 1.0 / (1.0 + 25.0*x[i]*x[i]) + ei;
    }

  /* construct least squares matrix with polynomial model */
  for (i = 0; i < n; ++i)
    {
      double Xij = 1.0;

      for (j = 0; j < p; ++j)
        {
          gsl_matrix_set(X, i, j, Xij);
          Xij *= x[i];
        }
    }

  /* least squares fits */
  {
    gsl_multifit_linear_workspace *w = gsl_multifit_linear_alloc(n, p);
    gsl_vector_view yv = gsl_vector_view_array(y, n);
    gsl_vector *c0 = gsl_vector_alloc(p);
    gsl_vector *c1 = gsl_vector_alloc(p);
    gsl_vector *c2 = gsl_vector_alloc(p);
    gsl_vector *c3 = gsl_vector_alloc(p);
    gsl_vector *lambda_vec = gsl_vector_calloc(p);
    gsl_matrix *cov = gsl_matrix_alloc(p, p);
    gsl_matrix *XTX = gsl_matrix_alloc(p, p); /* X^T X + lambda^2 I */
    gsl_vector *XTy = gsl_vector_alloc(p);    /* X^T y */
    gsl_vector_view xtx_diag = gsl_matrix_diagonal(XTX);
    gsl_permutation *perm = gsl_permutation_alloc(p);
    int signum;
    double chisq;

    /* construct XTy = X^T y */
    gsl_blas_dgemv(CblasTrans, 1.0, X, &yv.vector, 0.0, XTy);

    /* test that ridge equals OLS solution for lambda = 0 */
    gsl_multifit_linear(X, &yv.vector, c0, cov, &chisq, w);
    gsl_multifit_linear_ridge(0.0, X, &yv.vector, c1, cov, &chisq, w);

    /* test c0 = c1 */
    for (j = 0; j < p; ++j)
      {
        double c0j = gsl_vector_get(c0, j);
        double c1j = gsl_vector_get(c1, j);

        gsl_test_rel(c1j, c0j, 1.0e-10, "test_ridge: lambda = 0, c0/c1");
      }

    for (i = 0; i < 7; ++i)
      {
        double lambda = pow(10.0, -(double) i);

        gsl_multifit_linear_ridge(lambda, X, &yv.vector, c1, cov,
                                  &chisq, w);

        gsl_vector_set_all(lambda_vec, lambda);
        gsl_multifit_linear_ridge2(lambda_vec, X, &yv.vector, c2, cov,
                                   &chisq, w);

        /* construct XTX = X^T X + lamda^2 I */
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XTX);
        gsl_vector_add_constant(&xtx_diag.vector, lambda*lambda);

        /* solve XTX c = XTy with LU decomp */
        gsl_linalg_LU_decomp(XTX, perm, &signum);
        gsl_linalg_LU_solve(XTX, perm, XTy, c3);

        /* test c1 = c2 = c3 */
        for (j = 0; j < p; ++j)
          {
            double c1j = gsl_vector_get(c1, j);
            double c2j = gsl_vector_get(c2, j);
            double c3j = gsl_vector_get(c3, j);

            gsl_test_rel(c2j, c1j, 1.0e-10, "test_ridge: c2 lambda = %.1e", lambda);
            gsl_test_rel(c3j, c1j, 1.0e-9, "test_ridge: c3 lambda = %.1e", lambda);
          }

        /* now test a simple nontrivial L = diag(0.1,0.2,...) */

        /* XTX = X^T X */
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XTX);

        /* construct diag(L) and X^T X + L^T L */
        for (i = 0; i < p; ++i)
          {
            double val = (i + 1.0) / 10.0;

            gsl_vector_set(lambda_vec, i, val);
            *gsl_matrix_ptr(XTX, i, i) += val * val;
          }

        /* solve XTX c = XTy with LU decomp */
        gsl_linalg_LU_decomp(XTX, perm, &signum);
        gsl_linalg_LU_solve(XTX, perm, XTy, c1);

        /* solve with ridge routine */
        gsl_multifit_linear_ridge2(lambda_vec, X, &yv.vector, c2, cov,
                                   &chisq, w);

        /* test c1 = c2 */
        for (i = 0; i < p; ++i)
          {
            double c1i = gsl_vector_get(c1, i);
            double c2i = gsl_vector_get(c2, i);

            gsl_test_rel(c2i, c1i, 1.0e-12, "test_ridge: general L");
          }
      }

    gsl_multifit_linear_free(w);
    gsl_vector_free(c0);
    gsl_vector_free(c1);
    gsl_vector_free(c2);
    gsl_vector_free(c3);
    gsl_vector_free(lambda_vec);
    gsl_matrix_free(cov);
    gsl_matrix_free(XTX);
    gsl_vector_free(XTy);
    gsl_permutation_free(perm);
  }

  gsl_rng_free(r);
  free(x);
  free(y);
  gsl_matrix_free(X);
}
