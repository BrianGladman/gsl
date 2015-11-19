#define wood_N         6
#define wood_P         4

#define wood_NTRIES    3

static double wood_x0[wood_P] = { -3.0, -1.0, -3.0, -1.0 };
static double wood_epsrel = 1.0e-12;

static double wood_f[wood_N];
static double wood_J[wood_N * wood_P];

static void
wood_checksol(const double x[], const double sumsq,
              const double epsrel, const char *sname,
              const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;
  const double wood_x[wood_P] = { 1.0, 1.0, 1.0, 1.0 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < wood_P; ++i)
    {
      gsl_test_rel(x[i], wood_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
wood_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(wood_J, wood_N, wood_P);
  gsl_vector_view f = gsl_vector_view_array(wood_f, wood_N);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  double s90 = sqrt(90.0);
  double s10 = sqrt(10.0);

  gsl_vector_set(&f.vector, 0, 10.0*(x2 - x1*x1));
  gsl_vector_set(&f.vector, 1, 1.0 - x1);
  gsl_vector_set(&f.vector, 2, s90*(x4 - x3*x3));
  gsl_vector_set(&f.vector, 3, 1.0 - x3);
  gsl_vector_set(&f.vector, 4, s10*(x2 + x4 - 2.0));
  gsl_vector_set(&f.vector, 5, (x2 - x4) / s10);

  if (evaldf)
    {
      gsl_matrix_set_zero(&J.matrix);

      gsl_matrix_set(&J.matrix, 0, 0, -20.0*x1);
      gsl_matrix_set(&J.matrix, 0, 1, 10.0);
      gsl_matrix_set(&J.matrix, 1, 0, -1.0);
      gsl_matrix_set(&J.matrix, 2, 2, -2.0*s90*x3);
      gsl_matrix_set(&J.matrix, 2, 3, s90);
      gsl_matrix_set(&J.matrix, 3, 2, -1.0);
      gsl_matrix_set(&J.matrix, 4, 1, s10);
      gsl_matrix_set(&J.matrix, 4, 3, s10);
      gsl_matrix_set(&J.matrix, 5, 1, 1.0/s10);
      gsl_matrix_set(&J.matrix, 5, 3, -1.0/s10);
    }

  status = test_accumulate(2, &J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf wood_func =
{
  &wood_fdf,
  wood_P,
  NULL,
  0,
  0
};

static test_fdf_problem wood_problem =
{
  "wood",
  wood_x0,
  NULL,
  &wood_epsrel,
  wood_NTRIES,
  &wood_checksol,
  &wood_func
};
