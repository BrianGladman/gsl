static void
test_tsqr_system(const size_t n, const size_t p, const gsl_rng *r)
{
  const double tol = 1.0e-10;
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *c = gsl_vector_alloc(p);
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  char str[2048];
  size_t i;

  /* generate well-conditioned system and test against OLS solution */
  test_random_matrix(X, r, -1.0, 1.0);
  test_random_vector(c, r, -1.0, 1.0);

  /* compute y = X c + noise */
  gsl_blas_dgemv(CblasNoTrans, 1.0, X, c, 0.0, y);
  test_random_vector_noise(r, y);

  for (i = 0; i < 3; ++i)
    {
      double lambda = pow(10.0, -(double) i);

      /* solve system with multifit SVD approach */
      test_multifit_solve(lambda, X, y, c0);

      /* solve system with TSQR approach */
      test_multilarge_solve(gsl_multilarge_linear_tsqr, lambda, X, y, c1);

      sprintf(str, "tsqr n=%zu p=%zu lambda=%g", n, p, lambda);
      test_compare_vectors(tol, c0, c1, str);
    }

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_vector_free(c0);
  gsl_vector_free(c1);
}

/* test TSQR module */
static void
test_tsqr(const gsl_rng * r)
{
  test_tsqr_system(500, 350, r);
}
