#define helical_N         3
#define helical_P         3

#define helical_NTRIES    4

static double helical_x0[helical_P] = { -1.0, 0.0, 0.0 };
static double helical_x[helical_P] = { 1.0, 0.0, 0.0 };

static double helical_epsrel = 1.0e-12;

static double helical_f[helical_N];
static double helical_J[helical_N * helical_P];

static void
helical_checksol(const double x[], const double sumsq,
                 const double epsrel, const char *sname,
                 const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < helical_P; ++i)
    {
      gsl_test_rel(x[i], helical_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
helical_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(helical_J, helical_N, helical_P);
  gsl_vector_view f = gsl_vector_view_array(helical_f, helical_N);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double theta = (x1 >= 0.0) ? 0.0 : 5.0;
  double nx = gsl_hypot(x1, x2);

  gsl_vector_set(&f.vector, 0, 10.0 * (x3 - 5.0/M_PI*atan(x2 / x1) - theta));
  gsl_vector_set(&f.vector, 1, 10.0*(nx - 1.0));
  gsl_vector_set(&f.vector, 2, x3);

  if (evaldf)
    {
      double nx_sq = nx * nx;
      double term1 = 50.0 / (M_PI * nx_sq);
      double term2 = 10.0 / nx;

      gsl_matrix_set(&J.matrix, 0, 0, term1*x2);
      gsl_matrix_set(&J.matrix, 0, 1, -term1*x1);
      gsl_matrix_set(&J.matrix, 0, 2, 10.0);

      gsl_matrix_set(&J.matrix, 1, 0, term2*x1);
      gsl_matrix_set(&J.matrix, 1, 1, term2*x2);
      gsl_matrix_set(&J.matrix, 1, 2, 0.0);

      gsl_matrix_set(&J.matrix, 2, 0, 0.0);
      gsl_matrix_set(&J.matrix, 2, 1, 0.0);
      gsl_matrix_set(&J.matrix, 2, 2, 1.0);
    }

  status = test_accumulate(1, &J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf helical_func =
{
  &helical_fdf,
  helical_P,
  NULL,
  0,
  0
};

static test_fdf_problem helical_problem =
{
  "helical",
  helical_x0,
  NULL,
  &helical_epsrel,
  helical_NTRIES,
  &helical_checksol,
  &helical_func
};
