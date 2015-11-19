#define jennrich_N         10
#define jennrich_P         2

#define jennrich_NTRIES    1

static double jennrich_x0[jennrich_P] = { 0.3, 0.4 };
static double jennrich_epsrel = 1.0e-8;

static double jennrich_f[jennrich_N];
static double jennrich_J[jennrich_N * jennrich_P];

static void
jennrich_checksol(const double x[], const double sumsq,
                  const double epsrel, const char *sname,
                  const char *pname)
{
  size_t i;
  const double sumsq_exact = 1.243621823556148e+02;
  const double jennrich_x[jennrich_P] = { 2.578252139935855e-01,
                                          2.578252133471426e-01 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < jennrich_P; ++i)
    {
      gsl_test_rel(x[i], jennrich_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
jennrich_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(jennrich_J, jennrich_N, jennrich_P);
  gsl_vector_view f = gsl_vector_view_array(jennrich_f, jennrich_N);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  size_t i;

  for (i = 0; i < jennrich_N; ++i)
    {
      double ip1 = i + 1.0;
      double fi = 2.0*(i + 2.0) - (exp(x1*ip1) + exp(x2*ip1));

      gsl_vector_set(&f.vector, i, fi);

      if (evaldf)
        {
          gsl_matrix_set(&J.matrix, i, 0, -ip1*exp(ip1*x1));
          gsl_matrix_set(&J.matrix, i, 1, -ip1*exp(ip1*x2));
        }
    }

  status = test_accumulate(2, &J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf jennrich_func =
{
  &jennrich_fdf,
  jennrich_P,
  NULL,
  0,
  0
};

static test_fdf_problem jennrich_problem =
{
  "jennrich",
  jennrich_x0,
  NULL,
  &jennrich_epsrel,
  jennrich_NTRIES,
  &jennrich_checksol,
  &jennrich_func
};
