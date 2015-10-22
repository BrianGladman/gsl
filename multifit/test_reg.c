#include "oct.c"

/* generate random square orthogonal matrix via QR decomposition */
static void
test_random_matrix_orth(gsl_matrix *m, const gsl_rng *r)
{
  const size_t M = m->size1;
  gsl_matrix *A = gsl_matrix_alloc(M, M);
  gsl_vector *tau = gsl_vector_alloc(M);
  gsl_matrix *R = gsl_matrix_alloc(M, M);

  test_random_matrix(A, r, -1.0, 1.0);
  gsl_linalg_QR_decomp(A, tau);
  gsl_linalg_QR_unpack(A, tau, m, R);

  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_vector_free(tau);
}

/* construct ill-conditioned matrix via SVD */
static void
test_random_matrix_ill(gsl_matrix *m, const gsl_rng *r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  gsl_matrix *U = gsl_matrix_alloc(M, M);
  gsl_matrix *V = gsl_matrix_alloc(N, N);
  gsl_vector *S = gsl_vector_alloc(N);
  gsl_matrix_view Uv = gsl_matrix_submatrix(U, 0, 0, M, N);
  const double smin = 16.0 * GSL_DBL_EPSILON;
  const double smax = 10.0;
  const double ratio = pow(smin / smax, 1.0 / (N - 1.0));
  double s;
  size_t j;

  test_random_matrix_orth(U, r);
  test_random_matrix_orth(V, r);

  /* compute U * S */

  s = smax;
  for (j = 0; j < N; ++j)
    {
      gsl_vector_view uj = gsl_matrix_column(U, j);

      gsl_vector_scale(&uj.vector, s);
      s *= ratio;
    }

  /* compute m = (U * S) * V' */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &Uv.matrix, V, 0.0, m);

  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(S);
}

/* solve system with lambda = 0 and test against OLS solution */
static void
test_reg1(const gsl_matrix * X, const gsl_vector * y,
          const double tol, gsl_multifit_linear_workspace * w)
{
  const size_t p = X->size2;
  double rnorm, snorm, chisq;
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  gsl_matrix *cov = gsl_matrix_alloc(p, p);
  size_t j;

  /* test that reg. solution equals OLS solution for lambda = 0 */
  gsl_multifit_linear(X, y, c0, cov, &chisq, w);

  gsl_multifit_linear_svd(X, w);
  gsl_multifit_linear_solve(0.0, X, y, c1, &rnorm, &snorm, w);

  gsl_test_rel(rnorm*rnorm + snorm*snorm, chisq, tol,
               "test_reg1: lambda = 0, chisq");

  /* test c0 = c1 */
  for (j = 0; j < p; ++j)
    {
      double c0j = gsl_vector_get(c0, j);
      double c1j = gsl_vector_get(c1, j);

      gsl_test_rel(c1j, c0j, tol, "test_reg1: lambda = 0, c0/c1");
    }

  gsl_vector_free(c0);
  gsl_vector_free(c1);
  gsl_matrix_free(cov);
}

/* solve standard form system with given lambda and test against
 * normal equations solution */
static void
test_reg2(const double lambda, const gsl_matrix * X, const gsl_vector * y,
          const double tol, gsl_multifit_linear_workspace * w)
{
  const size_t p = X->size2;
  double rnorm, snorm;
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
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
  gsl_multifit_linear_solve(lambda, X, y, c1, &rnorm, &snorm, w);

  /* test c0 = c1 */
  for (j = 0; j < p; ++j)
    {
      double c0j = gsl_vector_get(c0, j);
      double c1j = gsl_vector_get(c1, j);

      gsl_test_rel(c1j, c0j, tol, "test_reg2: lambda=%g", lambda);
    }

  gsl_matrix_free(XTX);
  gsl_vector_free(XTy);
  gsl_vector_free(c0);
  gsl_vector_free(c1);
  gsl_permutation_free(perm);
}

/* solve system with given lambda and L = diag(L) and test against
 * normal equations solution */
static void
test_reg3(const double lambda, const gsl_vector * L, const gsl_matrix * X,
          const gsl_vector * y, const double tol, gsl_multifit_linear_workspace * w)
{
  const size_t n = X->size1;
  const size_t p = X->size2;
  double rnorm, snorm;
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  gsl_matrix *XTX = gsl_matrix_alloc(p, p); /* X^T X + lambda^2 L^T L */
  gsl_vector *XTy = gsl_vector_alloc(p);    /* X^T y */
  gsl_matrix *X2 = gsl_matrix_alloc(n, p);
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

  /* solve with reg routine */
  gsl_matrix_memcpy(X2, X);
  gsl_multifit_linear_stdform1(L, X2, w);
  gsl_multifit_linear_svd(X2, w);
  gsl_multifit_linear_solve(lambda, X2, y, c1, &rnorm, &snorm, w);
  gsl_multifit_linear_genform1(L, c1, w);

  /* test c0 = c1 */
  for (j = 0; j < p; ++j)
    {
      double c0j = gsl_vector_get(c0, j);
      double c1j = gsl_vector_get(c1, j);

      gsl_test_rel(c1j, c0j, tol, "test_reg3: lambda=%g diagonal L", lambda);
    }

  gsl_matrix_free(X2);
  gsl_matrix_free(XTX);
  gsl_vector_free(XTy);
  gsl_vector_free(c0);
  gsl_vector_free(c1);
  gsl_permutation_free(perm);
}

/* solve system with given lambda and L and test against
 * normal equations solution */
static void
test_reg4(const double lambda, const gsl_matrix * L, const gsl_matrix * X,
          const gsl_vector * y, const double tol, gsl_multifit_linear_workspace * w)
{
  const size_t p = X->size2;
  double rnorm, snorm;
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  gsl_matrix *LTL = gsl_matrix_alloc(p, p); /* L^T L */
  gsl_matrix *XTX = gsl_matrix_alloc(p, p); /* X^T X + lambda^2 L^T L */
  gsl_vector *XTy = gsl_vector_alloc(p);    /* X^T y */
  gsl_permutation *perm = gsl_permutation_alloc(p);
  gsl_matrix *Xs, *M;
  gsl_vector *ys, *cs;
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

  /* solve with reg routine */
  gsl_multifit_linear_stdform2(L, X, y, &Xs, &ys, &M, w);
  gsl_multifit_linear_svd(Xs, w);
  cs = gsl_vector_alloc(Xs->size2);
  gsl_multifit_linear_solve(lambda, Xs, ys, cs, &rnorm, &snorm, w);
  gsl_multifit_linear_genform2(L, X, y, cs, M, c1, w);

  /* test c0 = c1 */
  for (j = 0; j < p; ++j)
    {
      double c0j = gsl_vector_get(c0, j);
      double c1j = gsl_vector_get(c1, j);

      gsl_test_rel(c1j, c0j, tol, "test_reg4: lambda=%g general L", lambda);
    }

  gsl_matrix_free(LTL);
  gsl_matrix_free(XTX);
  gsl_vector_free(XTy);
  gsl_vector_free(c0);
  gsl_vector_free(c1);
  gsl_permutation_free(perm);
  gsl_matrix_free(Xs);
  gsl_vector_free(ys);
  gsl_vector_free(cs);

  if (M)
    gsl_matrix_free(M);
}

static void
test_reg_system(const size_t n, const size_t p, const gsl_rng *r)
{
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *c = gsl_vector_alloc(p);
  gsl_multifit_linear_workspace *w = gsl_multifit_linear_alloc(n, p);
  gsl_vector *diagL = gsl_vector_alloc(p);
  gsl_matrix *L1 = gsl_matrix_alloc(p, p);
  gsl_matrix *L2 = gsl_multifit_linear_L(p, 2);
  size_t i;

  /* generate well-conditioned system and test against OLS solution */
  test_random_matrix(X, r, -1.0, 1.0);
  test_random_vector(y, r, -1.0, 1.0);
  test_reg1(X, y, 1.0e-12, w);

  /* generate ill-conditioned system */
  test_random_matrix_ill(X, r);
  test_random_vector(c, r, -1.0, 1.0);

  /* compute y = X c + noise */
  gsl_blas_dgemv(CblasNoTrans, 1.0, X, c, 0.0, y);
  test_random_vector_noise(r, y);

  /* random diag(L) vector */
  test_random_vector(diagL, r, -2.0, 2.0);

  /* random square L matrix */
  test_random_matrix(L1, r, -2.0, 2.0);

  for (i = 0; i < 3; ++i)
    {
      /*
       * can't make lambda too small or normal equations
       * approach won't work well
       */
      double lambda = pow(10.0, -(double) i);

      test_reg2(lambda, X, y, 1.0e-7, w);
      test_reg3(lambda, diagL, X, y, 1.0e-7, w);
      test_reg4(lambda, L1, X, y, 1.0e-8, w);
      test_reg4(lambda, L2, X, y, 1.0e-7, w);
    }

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_vector_free(diagL);
  gsl_matrix_free(L1);
  gsl_matrix_free(L2);
  gsl_multifit_linear_free(w);
}

/* test linear regularized regression */
static void
test_reg(void)
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  test_reg_system(100, 10, r);
  test_reg_system(100, 50, r);
  test_reg_system(100, 100, r);

  gsl_rng_free(r);
}
