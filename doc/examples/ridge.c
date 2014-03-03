#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>

#define N   50 /* number of data */

int
main()
{
  const size_t n = N;
  const size_t p = 2;
  size_t i;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *y = gsl_vector_alloc(n);

  for (i = 0; i < n; ++i)
    {
      /* generate first random variable u */
      double ui = gsl_ran_gaussian(r, 1.0);

      /* set v = u + noise */
      double vi = ui + gsl_ran_gaussian(r, 0.001);

      /* set y = u + v + noise */
      double yi = ui + vi + gsl_ran_gaussian(r, 1.0);

      /* since u =~ v, the matrix X is ill-conditioned */
      gsl_matrix_set(X, i, 0, ui);
      gsl_matrix_set(X, i, 1, vi);

      /* rhs vector */
      gsl_vector_set(y, i, yi);
    }

  {
    gsl_multifit_linear_workspace *w =
      gsl_multifit_linear_alloc(n, p);
    gsl_vector *c = gsl_vector_alloc(p);
    gsl_vector *c_ridge = gsl_vector_alloc(p);
    gsl_matrix *cov = gsl_matrix_alloc(p, p);
    double chisq;

    /* unregularized (standard) least squares fit, lambda = 0 */
    gsl_multifit_linear_ridge(0.0, X, y, c, cov, &chisq, w);

    fprintf(stderr, "=== Unregularized fit ===\n");
    fprintf(stderr, "best fit: y = %g u + %g v\n",
      gsl_vector_get(c, 0), gsl_vector_get(c, 1));
    fprintf(stderr, "chisq/dof = %g\n", chisq / (n - p));

    /* regularize with lambda = 1 */
    gsl_multifit_linear_ridge(1.0, X, y, c_ridge, cov, &chisq, w);

    fprintf(stderr, "=== Regularized fit ===\n");
    fprintf(stderr, "best fit: y = %g u + %g v\n",
      gsl_vector_get(c_ridge, 0), gsl_vector_get(c_ridge, 1));
    fprintf(stderr, "chisq/dof = %g\n", chisq / (n - p));

    gsl_multifit_linear_free(w);
    gsl_matrix_free(cov);
    gsl_vector_free(c);
    gsl_vector_free(c_ridge);
  }

  gsl_rng_free(r);
  gsl_matrix_free(X);
  gsl_vector_free(y);

  return 0;
}
