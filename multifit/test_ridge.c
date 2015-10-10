/* solve system with lambda = 0 and test against OLS solution */
static void
test_ridge1(const gsl_matrix * X, const gsl_vector * y,
            const double tol, gsl_multifit_linear_workspace * w)
{
  const size_t p = X->size2;
  double rnorm, snorm, chisq;
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  gsl_matrix *cov = gsl_matrix_alloc(p, p);
  size_t j;

  /* test that ridge equals OLS solution for lambda = 0 */
  gsl_multifit_linear(X, y, c0, cov, &chisq, w);

  gsl_multifit_linear_svd(X, w);
  gsl_multifit_linear_ridge_solve(0.0, y, c1, cov, &rnorm, &snorm, w);

  gsl_test_rel(rnorm*rnorm + snorm*snorm, chisq, tol,
               "test_ridge1: lambda = 0, chisq");

  /* test c0 = c1 */
  for (j = 0; j < p; ++j)
    {
      double c0j = gsl_vector_get(c0, j);
      double c1j = gsl_vector_get(c1, j);

      gsl_test_rel(c1j, c0j, tol, "test_ridge1: lambda = 0, c0/c1");
    }

  gsl_vector_free(c0);
  gsl_vector_free(c1);
  gsl_matrix_free(cov);
}

/* solve standard form system with given lambda and test against
 * normal equations solution */
static void
test_ridge2(const double lambda, const gsl_matrix * X, const gsl_vector * y,
            const double tol, gsl_multifit_linear_workspace * w)
{
  const size_t p = X->size2;
  double rnorm, snorm;
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  gsl_matrix *cov = gsl_matrix_alloc(p, p);
  gsl_matrix *XTX = gsl_matrix_alloc(p, p); /* X^T X + lambda^2 I */
  gsl_vector *XTy = gsl_vector_alloc(p);    /* X^T y */
  gsl_vector_view xtx_diag = gsl_matrix_diagonal(XTX);
  gsl_permutation *perm = gsl_permutation_alloc(p);
  int signum;
  size_t j;

  /* construct XTy = X^T y */
  gsl_blas_dgemv(CblasTrans, 1.0, X, y, 0.0, XTy);

  /* construct XTX = X^T X + lambda^2 I */
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XTX);
  gsl_vector_add_constant(&xtx_diag.vector, lambda*lambda);

  /* solve XTX c = XTy with LU decomp */
  gsl_linalg_LU_decomp(XTX, perm, &signum);
  gsl_linalg_LU_solve(XTX, perm, XTy, c0);

  /* compute SVD of X */
  gsl_multifit_linear_svd(X, w);

  /* solve regularized standard form system with lambda */
  gsl_multifit_linear_ridge_solve(lambda, y, c1, cov, &rnorm, &snorm, w);

  /* test c0 = c1 */
  for (j = 0; j < p; ++j)
    {
      double c0j = gsl_vector_get(c0, j);
      double c1j = gsl_vector_get(c1, j);

      gsl_test_rel(c1j, c0j, tol, "test_ridge2: lambda=%g", lambda);
    }

  gsl_matrix_free(XTX);
  gsl_vector_free(XTy);
  gsl_vector_free(c0);
  gsl_vector_free(c1);
  gsl_matrix_free(cov);
  gsl_permutation_free(perm);
}

/* solve system with given lambda and L = diag(L) and test against
 * normal equations solution */
static void
test_ridge3(const double lambda, const gsl_vector * L, const gsl_matrix * X,
            const gsl_vector * y, const double tol, gsl_multifit_linear_workspace * w)
{
  const size_t p = X->size2;
  double rnorm, snorm;
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  gsl_matrix *cov = gsl_matrix_alloc(p, p);
  gsl_matrix *XTX = gsl_matrix_alloc(p, p); /* X^T X + lambda^2 L^T L */
  gsl_vector *XTy = gsl_vector_alloc(p);    /* X^T y */
  gsl_permutation *perm = gsl_permutation_alloc(p);
  int signum;
  size_t j;

  /* construct XTy = X^T y */
  gsl_blas_dgemv(CblasTrans, 1.0, X, y, 0.0, XTy);

  /* construct XTX = X^T X + lambda^2 L^T L */
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XTX);

  for (j = 0; j < p; ++j)
    {
      double lj = gsl_vector_get(L, j);
      *gsl_matrix_ptr(XTX, j, j) += pow(lambda * lj, 2.0);
    }

  /* solve XTX c = XTy with LU decomp */
  gsl_linalg_LU_decomp(XTX, perm, &signum);
  gsl_linalg_LU_solve(XTX, perm, XTy, c0);

  /* solve with ridge routine */
  gsl_multifit_linear_ridge_svd(X, L, w);
  gsl_multifit_linear_ridge_solve(lambda, y, c1, cov, &rnorm, &snorm, w);
  gsl_multifit_linear_ridge_transform(L, c1, w);

  /* test c0 = c1 */
  for (j = 0; j < p; ++j)
    {
      double c0j = gsl_vector_get(c0, j);
      double c1j = gsl_vector_get(c1, j);

      gsl_test_rel(c1j, c0j, tol, "test_ridge3: lambda=%g diagonal L", lambda);
    }

  gsl_matrix_free(XTX);
  gsl_vector_free(XTy);
  gsl_vector_free(c0);
  gsl_vector_free(c1);
  gsl_matrix_free(cov);
  gsl_permutation_free(perm);
}

/* solve system with given lambda and L and test against
 * normal equations solution */
static void
test_ridge4(const double lambda, const gsl_matrix * L, const gsl_matrix * X,
            const gsl_vector * y, const double tol, gsl_multifit_linear_workspace * w)
{
  const size_t p = X->size2;
  double rnorm, snorm;
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  gsl_matrix *cov = gsl_matrix_alloc(p, p);
  gsl_matrix *QR = gsl_matrix_alloc(p, p);
  gsl_vector *tau = gsl_vector_alloc(p);
  gsl_matrix *LTL = gsl_matrix_alloc(p, p); /* L^T L */
  gsl_matrix *XTX = gsl_matrix_alloc(p, p); /* X^T X + lambda^2 L^T L */
  gsl_vector *XTy = gsl_vector_alloc(p);    /* X^T y */
  gsl_permutation *perm = gsl_permutation_alloc(p);
  int signum;
  size_t j;

  /* construct XTy = X^T y */
  gsl_blas_dgemv(CblasTrans, 1.0, X, y, 0.0, XTy);

  /* construct XTX = X^T X + lambda^2 L^T L */
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L, L, 0.0, LTL);
  gsl_matrix_scale(LTL, lambda * lambda);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XTX);
  gsl_matrix_add(XTX, LTL);

  /* solve XTX c = XTy with LU decomp */
  gsl_linalg_LU_decomp(XTX, perm, &signum);
  gsl_linalg_LU_solve(XTX, perm, XTy, c0);

  /* solve with ridge routine */
  gsl_matrix_memcpy(QR, L);
  gsl_multifit_linear_ridge_svd2(X, QR, tau, w);
  gsl_multifit_linear_ridge_solve(lambda, y, c1, cov, &rnorm, &snorm, w);
  gsl_multifit_linear_ridge_transform2(QR, tau, c1, w);

  /* test c0 = c1 */
  for (j = 0; j < p; ++j)
    {
      double c0j = gsl_vector_get(c0, j);
      double c1j = gsl_vector_get(c1, j);

      gsl_test_rel(c1j, c0j, tol, "test_ridge4: lambda=%g general L", lambda);
    }

  gsl_matrix_free(QR);
  gsl_vector_free(tau);
  gsl_matrix_free(LTL);
  gsl_matrix_free(XTX);
  gsl_vector_free(XTy);
  gsl_vector_free(c0);
  gsl_vector_free(c1);
  gsl_matrix_free(cov);
  gsl_permutation_free(perm);
}

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
  gsl_vector_view yv = gsl_vector_view_array(y, n);
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_multifit_linear_workspace *w = gsl_multifit_linear_alloc(n, p);
  gsl_vector *Ldiag = gsl_vector_alloc(p);
  gsl_matrix *L = gsl_matrix_alloc(p, p);
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

  test_ridge1(X, &yv.vector, 1.0e-10, w);

  /* build diag(L) vector */
  gsl_vector_set_all(Ldiag, 2.5);

  /* build general matrix L as second difference operator */
  gsl_matrix_set_zero(L);

  /* first row [2 -1 0 ... ] */
  gsl_matrix_set(L, 0, 0, 2.0);
  gsl_matrix_set(L, 0, 1, -1.0);

  for (i = 1; i < p - 1; ++i)
    {
      gsl_matrix_set(L, i, i - 1, -1.0);
      gsl_matrix_set(L, i, i, 2.0);
      gsl_matrix_set(L, i, i + 1, -1.0);
    }

  /* last row [0 ... 0 -1 2] */
  gsl_matrix_set(L, p - 1, p - 2, -1.0);
  gsl_matrix_set(L, p - 1, p - 1, 2.0);

  for (i = 0; i < 7; ++i)
    {
      double lambda = pow(10.0, -(double) i);

      test_ridge2(lambda, X, &yv.vector, 1.0e-9, w);
      test_ridge3(lambda, Ldiag, X, &yv.vector, 1.0e-9, w);
      test_ridge4(lambda, L, X, &yv.vector, 1.0e-8, w);
    }

  gsl_multifit_linear_free(w);
  gsl_rng_free(r);
  free(x);
  free(y);
  gsl_matrix_free(X);
  gsl_vector_free(Ldiag);
  gsl_matrix_free(L);
}
