#define brown3_N         3
#define brown3_P         2

#define brown3_NTRIES    3

static double brown3_x0[brown3_P] = { 1.0, 1.0 };
static double brown3_epsrel = 1.0e-12;

static double brown3_f[brown3_N];
static double brown3_J[brown3_N * brown3_P];

static void
brown3_checksol(const double x[], const double sumsq,
                const double epsrel, const char *sname,
                const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;
  const double brown3_x[brown3_P] = { 1.0e6, 2.0e-6 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < brown3_P; ++i)
    {
      gsl_test_rel(x[i], brown3_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
brown3_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(brown3_J, brown3_N, brown3_P);
  gsl_vector_view f = gsl_vector_view_array(brown3_f, brown3_N);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);

  gsl_vector_set(&f.vector, 0, x1 - 1.0e6);
  gsl_vector_set(&f.vector, 1, x2 - 2.0e-6);
  gsl_vector_set(&f.vector, 2, x1*x2 - 2.0);

  if (evaldf)
    {
      gsl_matrix_set_zero(&J.matrix);
      gsl_matrix_set(&J.matrix, 0, 0, 1.0);
      gsl_matrix_set(&J.matrix, 1, 1, 1.0);
      gsl_matrix_set(&J.matrix, 2, 0, x2);
      gsl_matrix_set(&J.matrix, 2, 1, x1);
    }

  status = gsl_multilarge_nlinear_accumulate(&J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf brown3_func =
{
  &brown3_fdf,
  brown3_P,
  NULL,
  0,
  0
};

static test_fdf_problem brown3_problem =
{
  "brown_badly_scaled",
  brown3_x0,
  NULL,
  &brown3_epsrel,
  brown3_NTRIES,
  &brown3_checksol,
  &brown3_func
};
