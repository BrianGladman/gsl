#define powell2_N        2
#define powell2_P        2

#define powell2_NTRIES   3

static double powell2_x0[powell2_P] = { 3.0, 1.0 };
static double powell2_epsrel = 1.0e-6;

static double powell2_f[powell2_N];
static double powell2_J[powell2_N * powell2_P];

static void
powell2_checksol(const double x[], const double sumsq,
                 const double epsrel, const char *sname,
                 const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < powell2_P; ++i)
    {
      gsl_test_rel(x[i], 0.0, epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
powell2_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(powell2_J, powell2_N, powell2_P);
  gsl_vector_view f = gsl_vector_view_array(powell2_f, powell2_N);
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);
  double term = x0 + 0.1;

  gsl_vector_set(&f.vector, 0, x0);
  gsl_vector_set(&f.vector, 1, 10.0*x0/(x0 + 0.1) + 2.0*x1*x1);

  if (evaldf)
    {
      gsl_matrix_set(&J.matrix, 0, 0, 1.0);
      gsl_matrix_set(&J.matrix, 0, 1, 0.0);
      gsl_matrix_set(&J.matrix, 1, 0, 1.0 / (term * term));
      gsl_matrix_set(&J.matrix, 1, 1, 4.0 * x1);
    }

  status = gsl_multilarge_nlinear_accumulate(&J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf powell2_func =
{
  &powell2_fdf,
  powell2_P,
  NULL,
  0,
  0
};

static test_fdf_problem powell2_problem =
{
  "powell2",
  powell2_x0,
  NULL,
  &powell2_epsrel,
  powell2_NTRIES,
  &powell2_checksol,
  &powell2_func
};
