#define beale_N         3
#define beale_P         2

#define beale_NTRIES    1

static double beale_x0[beale_P] = { 1.0, 1.0 };
static double beale_epsrel = 1.0e-12;

static double beale_f[beale_N];
static double beale_J[beale_N * beale_P];

static double beale_Y[beale_N] = { 1.5, 2.25, 2.625 };

static void
beale_checksol(const double x[], const double sumsq,
               const double epsrel, const char *sname,
               const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;
  const double beale_x[beale_P] = { 3.0, 0.5 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < beale_P; ++i)
    {
      gsl_test_rel(x[i], beale_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
beale_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(beale_J, beale_N, beale_P);
  gsl_vector_view f = gsl_vector_view_array(beale_f, beale_N);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  size_t i;

  for (i = 0; i < beale_N; ++i)
    {
      double yi = beale_Y[i];
      double term = pow(x2, (double) i);
      double fi = yi - x1*(1.0 - term*x2);

      gsl_vector_set(&f.vector, i, fi);

      if (evaldf)
        {
          gsl_matrix_set(&J.matrix, i, 0, term*x2 - 1.0);
          gsl_matrix_set(&J.matrix, i, 1, (i + 1.0) * x1 * term);
        }
    }

  status = test_accumulate(2, &J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf beale_func =
{
  &beale_fdf,
  beale_P,
  NULL,
  0,
  0
};

static test_fdf_problem beale_problem =
{
  "beale",
  beale_x0,
  NULL,
  &beale_epsrel,
  beale_NTRIES,
  &beale_checksol,
  &beale_func
};
