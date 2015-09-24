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
    gsl_multifit_linear_workspace *w2 = gsl_multifit_linear_alloc(n, p);
    gsl_vector_view yv = gsl_vector_view_array(y, n);
    gsl_vector *c0 = gsl_vector_alloc(p);
    gsl_vector *c1 = gsl_vector_alloc(p);
    gsl_vector *c2 = gsl_vector_alloc(p);
    gsl_vector *c3 = gsl_vector_alloc(p);
    gsl_vector *g = gsl_vector_calloc(p);
    gsl_matrix *cov = gsl_matrix_alloc(p, p);
    gsl_matrix *XTX = gsl_matrix_alloc(p, p); /* X^T X + lambda^2 I */
    gsl_vector *XTy = gsl_vector_alloc(p);    /* X^T y */
    gsl_vector_view xtx_diag = gsl_matrix_diagonal(XTX);
    gsl_permutation *perm = gsl_permutation_alloc(p);
    int signum;
    double chisq, rnorm, snorm;

    /* construct XTy = X^T y */
    gsl_blas_dgemv(CblasTrans, 1.0, X, &yv.vector, 0.0, XTy);

    /* test that ridge equals OLS solution for lambda = 0 */
    gsl_multifit_linear(X, &yv.vector, c0, cov, &chisq, w);

    gsl_multifit_linear_svd(X, w);
    gsl_multifit_linear_ridge_solve(0.0, &yv.vector, c1, cov,
                                    &rnorm, &snorm, w);

    gsl_test_rel(rnorm*rnorm + snorm*snorm, chisq, 1.0e-10,
                 "test_ridge: lambda = 0, chisq");

    /* test c0 = c1 */
    for (j = 0; j < p; ++j)
      {
        double c0j = gsl_vector_get(c0, j);
        double c1j = gsl_vector_get(c1, j);

        gsl_test_rel(c1j, c0j, 1.0e-10, "test_ridge: lambda = 0, c0/c1");
      }

    /* compute SVD of X */
    gsl_multifit_linear_svd(X, w);

    /* solve regularized standard form systems with different lambdas */
    for (i = 0; i < 7; ++i)
      {
        double lambda = pow(10.0, -(double) i);

        gsl_multifit_linear_ridge_solve(lambda, &yv.vector, c1, cov,
                                        &rnorm, &snorm, w);

        /* compute SVD of transformed system */
        gsl_vector_set_all(g, 1.0);
        gsl_multifit_linear_ridge_svd(g, X, w2);
        gsl_multifit_linear_ridge_solve(lambda, &yv.vector, c2, cov,
                                        &rnorm, &snorm, w2);
        gsl_multifit_linear_ridge_trans(g, c2, w2);

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

            gsl_test_rel(c2j, c1j, 1.0e-10, "test_ridge: lambda=%g c2", lambda);
            gsl_test_rel(c3j, c1j, 1.0e-9, "test_ridge: lambda=%g c3", lambda);
          }

        /* now test a simple nontrivial system:
         * L = lambda * diag(0.1,0.2,...)
         */

        /* XTX = X^T X */
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XTX);

        /* construct diag(G) and X^T X + lambda * G^T G */
        for (j = 0; j < p; ++j)
          {
            double val = (j + 1.0) / 10.0;

            gsl_vector_set(g, j, val);
            *gsl_matrix_ptr(XTX, j, j) += pow(lambda * val, 2.0);
          }

        /* solve XTX c = XTy with LU decomp */
        gsl_linalg_LU_decomp(XTX, perm, &signum);
        gsl_linalg_LU_solve(XTX, perm, XTy, c1);

        /* solve with ridge routine */
        gsl_multifit_linear_ridge_svd(g, X, w2);
        gsl_multifit_linear_ridge_solve(lambda, &yv.vector, c2, cov,
                                        &rnorm, &snorm, w2);
        gsl_multifit_linear_ridge_trans(g, c2, w2);

        /* test c1 = c2 */
        for (j = 0; j < p; ++j)
          {
            double c1j = gsl_vector_get(c1, j);
            double c2j = gsl_vector_get(c2, j);

            gsl_test_rel(c2j, c1j, 1.0e-9, "test_ridge: lambda=%g general L", lambda);
          }
      }

    gsl_multifit_linear_free(w);
    gsl_multifit_linear_free(w2);
    gsl_vector_free(c0);
    gsl_vector_free(c1);
    gsl_vector_free(c2);
    gsl_vector_free(c3);
    gsl_vector_free(g);
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
