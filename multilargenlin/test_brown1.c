#define brown1_N         20
#define brown1_P         4

#define brown1_NTRIES    3

static double brown1_x0[brown1_P] = { 25, 5, -5, -1 };
static double brown1_epsrel = 1.0e-4;

static double brown1_f[brown1_N];
static double brown1_J[brown1_N * brown1_P];

static void
brown1_checksol(const double x[], const double sumsq,
                const double epsrel, const char *sname,
                const char *pname)
{
  size_t i;
  const double sumsq_exact = 8.582220162635628e+04;
  const double brown1_x[brown1_P] = {
    -1.159443990239263e+01, 1.320363005221244e+01,
    -4.034395456782477e-01, 2.367789088597534e-01 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < brown1_P; ++i)
    {
      gsl_test_rel(x[i], brown1_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
brown1_fdf (const gsl_vector * x, gsl_matrix * JTJ,
            gsl_vector * JTf, double * normf, void *params)
{
  gsl_matrix_view J = gsl_matrix_view_array(brown1_J, brown1_N, brown1_P);
  gsl_vector_view f = gsl_vector_view_array(brown1_f, brown1_N);
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);
  double x2 = gsl_vector_get (x, 2);
  double x3 = gsl_vector_get (x, 3);
  size_t i;

  for (i = 0; i < brown1_N; i++)
    {
      double ti = 0.2 * (i + 1);
      double ui = x0 + x1 * ti - exp (ti);
      double vi = x2 + x3 * sin (ti) - cos (ti);

      gsl_vector_set (&f.vector, i, ui * ui + vi * vi);

      if (JTJ)
        {
          gsl_matrix_set (&J.matrix, i, 0, 2 * ui);
          gsl_matrix_set (&J.matrix, i, 1, 2 * ui * ti);
          gsl_matrix_set (&J.matrix, i, 2, 2 * vi);
          gsl_matrix_set (&J.matrix, i, 3, 2 * vi * sin (ti));
        }
    }

  *normf = gsl_blas_dnrm2(&f.vector);

  if (JTJ)
    {
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);
      gsl_blas_dgemv(CblasTrans, 1.0, &J.matrix, &f.vector, 0.0, JTf);
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multilarge_function_fdf brown1_func =
{
  &brown1_fdf,
  brown1_P,
  NULL,
  0,
  0
};

static test_fdf_problem brown1_problem =
{
  "brown_dennis",
  brown1_x0,
  NULL,
  &brown1_epsrel,
  brown1_NTRIES,
  &brown1_checksol,
  &brown1_func
};
