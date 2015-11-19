#define rosenbrock_N         2
#define rosenbrock_P         2

#define rosenbrock_NTRIES    4

static double rosenbrock_x0[rosenbrock_P] = { -1.2, 1.0 };
static double rosenbrock_epsrel = 1.0e-12;

static double rosenbrock_f[rosenbrock_N];
static double rosenbrock_J[rosenbrock_N * rosenbrock_P];

static void
rosenbrock_checksol(const double x[], const double sumsq,
                    const double epsrel, const char *sname,
                    const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < rosenbrock_P; ++i)
    {
      gsl_test_rel(x[i], 1.0, epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
rosenbrock_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(rosenbrock_J, rosenbrock_N, rosenbrock_P);
  gsl_vector_view f = gsl_vector_view_array(rosenbrock_f, rosenbrock_N);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);

  gsl_vector_set(&f.vector, 0, 10.0 * (x2 - x1*x1));
  gsl_vector_set(&f.vector, 1, 1.0 - x1);

  if (evaldf)
    {
      gsl_matrix_set(&J.matrix, 0, 0, -20.0*x1);
      gsl_matrix_set(&J.matrix, 0, 1, 10.0);
      gsl_matrix_set(&J.matrix, 1, 0, -1.0);
      gsl_matrix_set(&J.matrix, 1, 1, 0.0);
    }

  status = test_accumulate(2, &J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf rosenbrock_func =
{
  &rosenbrock_fdf,
  rosenbrock_P,
  NULL,
  0,
  0
};

static test_fdf_problem rosenbrock_problem =
{
  "rosenbrock",
  rosenbrock_x0,
  NULL,
  &rosenbrock_epsrel,
  rosenbrock_NTRIES,
  &rosenbrock_checksol,
  &rosenbrock_func
};
