/*
 * test_shaw.c
 *
 * Test L-curve (Tikhonov) regression routines using Shaw
 * problem. See example 1.10 of
 *
 * [1] R.C. Aster, B. Borchers and C. H. Thurber,
 *     Parameter Estimation and Inverse Problems (2nd ed), 2012.
 */

#include <gsl/gsl_sf_trig.h>

#include "oct.c"

static size_t shaw_n = 20;
static size_t shaw_p = 20;

static double shaw_y[] = {
 8.7547455124379323e-04, 5.4996835885761936e-04, 1.7527999407005367e-06,
 1.9552372913117047e-03, 1.4411645433785081e-02, 5.2800013336393704e-02,
 1.3609152023257112e-01, 2.7203484587635818e-01, 4.3752225136193390e-01,
 5.7547667319875240e-01, 6.2445052213539942e-01, 5.6252658286441348e-01,
 4.2322239923561566e-01, 2.6768469219560631e-01, 1.4337901162734543e-01,
 6.5614569346074361e-02, 2.6013851831752945e-02, 9.2336933089481269e-03,
 3.2269066658993694e-03, 1.3999201459261811e-03
};

/* construct design matrix for Shaw problem */
static int
shaw_matrix(gsl_matrix * X)
{
  int s = GSL_SUCCESS;
  const size_t n = X->size1;
  const size_t p = X->size2;
  const double dtheta = M_PI / (double) p;
  size_t i, j;

  for (i = 0; i < n; ++i)
    {
      double si = (i + 0.5) * M_PI / n - M_PI / 2.0;
      double csi = cos(si);
      double sni = sin(si);

      for (j = 0; j < p; ++j)
        {
          double thetaj = (j + 0.5) * M_PI / p - M_PI / 2.0;
          double term1 = csi + cos(thetaj);
          double term2 = gsl_sf_sinc(sni + sin(thetaj));
          double Xij = term1 * term1 * term2 * term2 * dtheta;

          gsl_matrix_set(X, i, j, Xij);
        }
    }

  return s;
}

void
test_shaw(void)
{
  const size_t npoints = 1000; /* number of points on L-curve */
  gsl_vector * reg_param = gsl_vector_alloc(npoints);
  gsl_vector * rho = gsl_vector_alloc(npoints);
  gsl_vector * eta = gsl_vector_alloc(npoints);

  gsl_matrix * X = gsl_matrix_alloc(shaw_n, shaw_p);
  gsl_vector_view y = gsl_vector_view_array(shaw_y, shaw_n);
  gsl_multifit_linear_workspace * work = 
    gsl_multifit_linear_alloc (shaw_n, shaw_p);

  size_t reg_idx;
  double lambda;

  /* build design matrix */
  shaw_matrix(X);

  print_octave(X, "Shaw_X");

  /* SVD decomposition */
  gsl_multifit_linear_ridge_svd(X, work);

  /* calculate L-curve */
  gsl_multifit_linear_ridge_lcurve(&y.vector, reg_param, rho, eta, work);

  /* calculate corner of L-curve */
  gsl_multifit_linear_ridge_lcorner(rho, eta, &reg_idx);

  lambda = gsl_vector_get(reg_param, reg_idx);

  /* test against value from [1] */
  gsl_test_rel(lambda, 5.793190958069266e-06, 1.0e-12, "shaw: L-curve corner lambda");

  gsl_matrix_free(X);
  gsl_vector_free(reg_param);
  gsl_vector_free(rho);
  gsl_vector_free(eta);
  gsl_multifit_linear_free(work);
} /* test_shaw() */
