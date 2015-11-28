#define gaussian_N         15
#define gaussian_P         3

#define gaussian_NTRIES    2

static double gaussian_x0[gaussian_P] = { 0.4, 1.0, 0.0 };
static double gaussian_epsrel = 1.0e-10;

static double gaussian_f[gaussian_N];
static double gaussian_J[gaussian_N * gaussian_P];

static double gaussian_Y[gaussian_N] = {
0.0009, 0.0044, 0.0175, 0.0540, 0.1295, 0.2420, 0.3521, 0.3989,
0.3521, 0.2420, 0.1295, 0.0540, 0.0175, 0.0044, 0.0009
};

static void
gaussian_checksol(const double x[], const double sumsq,
                  const double epsrel, const char *sname,
                  const char *pname)
{
  size_t i;
  const double sumsq_exact = 1.12793276961871985e-08;
  const double gaussian_x[gaussian_P] = { 0.398956137838762825,
                                          1.00001908448786647,
                                          0.0 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < gaussian_P; ++i)
    {
      gsl_test_rel(x[i], gaussian_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
gaussian_fdf (const gsl_vector * x, gsl_matrix * JTJ,
              gsl_vector * JTf, double * normf, void *params)
{
  gsl_matrix_view J = gsl_matrix_view_array(gaussian_J, gaussian_N, gaussian_P);
  gsl_vector_view f = gsl_vector_view_array(gaussian_f, gaussian_N);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  size_t i;

  for (i = 0; i < gaussian_N; ++i)
    {
      double ti = (7.0 - i) / 2.0;
      double yi = gaussian_Y[i];
      double term1 = ti - x3;
      double term2 = exp(-x2*term1*term1/2.0);
      double fi = x1 * exp(-x2*term1*term1/2.0) - yi;

      gsl_vector_set(&f.vector, i, fi);

      if (JTJ)
        {
          gsl_matrix_set(&J.matrix, i, 0, term2);
          gsl_matrix_set(&J.matrix, i, 1, -0.5*x1*term2*term1*term1);
          gsl_matrix_set(&J.matrix, i, 2, x1*x2*term1*term2);
        }
    }

  *normf = gsl_blas_dnrm2(&f.vector);

  if (JTJ)
    {
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);
      gsl_blas_dgemv(CblasTrans, 1.0, &J.matrix, &f.vector, 0.0, JTf);
    }

  return GSL_SUCCESS;
}

static gsl_multilarge_function_fdf gaussian_func =
{
  &gaussian_fdf,
  gaussian_P,
  NULL,
  0,
  0
};

static test_fdf_problem gaussian_problem =
{
  "gaussian",
  gaussian_x0,
  NULL,
  &gaussian_epsrel,
  gaussian_NTRIES,
  &gaussian_checksol,
  &gaussian_func
};
