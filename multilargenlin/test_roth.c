#define roth_N         2
#define roth_P         2

#define roth_NTRIES    3

static double roth_x0[roth_P] = { 0.5, -2.0 };
static double roth_epsrel = 1.0e-8;

static double roth_f[roth_N];
static double roth_J[roth_N * roth_P];

static void
roth_checksol(const double x[], const double sumsq,
              const double epsrel, const char *sname,
              const char *pname)
{
  size_t i;
  const double sumsq_exact1 = 0.0;
  const double roth_x1[roth_P] = { 5.0, 4.0 };
  const double sumsq_exact2 = 48.9842536792400;
  const double roth_x2[roth_P] = { 11.4127789869021, -0.896805253274477 };
  const double *roth_x;
  double sumsq_exact;

  if (fabs(sumsq) < 0.1)
    {
      sumsq_exact = sumsq_exact1;
      roth_x = roth_x1;
    }
  else
    {
      sumsq_exact = sumsq_exact2;
      roth_x = roth_x2;
    }

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < roth_P; ++i)
    {
      gsl_test_rel(x[i], roth_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
roth_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(roth_J, roth_N, roth_P);
  gsl_vector_view f = gsl_vector_view_array(roth_f, roth_N);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);

  gsl_vector_set(&f.vector, 0, x1 - x2*(2.0 - x2*(5.0 - x2)) - 13.0);
  gsl_vector_set(&f.vector, 1, x1 - x2*(14.0 - x2*(1.0 + x2)) - 29.0);

  if (evaldf)
    {
      gsl_matrix_set(&J.matrix, 0, 0, 1.0);
      gsl_matrix_set(&J.matrix, 0, 1, -2.0 + x2*(10.0 - 3.0*x2));
      gsl_matrix_set(&J.matrix, 1, 0, 1.0);
      gsl_matrix_set(&J.matrix, 1, 1, -14.0 + x2*(2.0 + 3.0*x2));
    }

  status = test_accumulate(2, &J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf roth_func =
{
  &roth_fdf,
  roth_P,
  NULL,
  0,
  0
};

static test_fdf_problem roth_problem =
{
  "roth_freudenstein",
  roth_x0,
  NULL,
  &roth_epsrel,
  roth_NTRIES,
  &roth_checksol,
  &roth_func
};
