#define biggs_N         6  /* >= p */
#define biggs_P         6

static double biggs_x0[biggs_P] = { 1.0, 2.0, 1.0, 1.0, 1.0, 1.0 };
static double biggs_x[biggs_P] = { 1.0, 10.0, 1.0, 5.0, 4.0, 3.0 };

static double biggs_sumsq = 0.0;
static double biggs_epsrel = 1.0e-9;

static int
biggs_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  double x5 = gsl_vector_get(x, 4);
  double x6 = gsl_vector_get(x, 5);
  size_t i;

  for (i = 0; i < biggs_N; ++i)
    {
      double ti = 0.1 * (i + 1.0);
      double yi = exp(-ti) - 5*exp(-10*ti) + 3*exp(-4*ti);
      double fi = x3*exp(-ti*x1) - x4*exp(-ti*x2) + x6*exp(-ti*x5) - yi;

      gsl_vector_set(f, i, fi);
    }

  return GSL_SUCCESS;
}

static int
biggs_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  double x5 = gsl_vector_get(x, 4);
  double x6 = gsl_vector_get(x, 5);
  size_t i;

  for (i = 0; i < biggs_N; ++i)
    {
      double ti = 0.1 * (i + 1.0);

      gsl_matrix_set(J, i, 0, -ti*x3*exp(-ti*x1));
      gsl_matrix_set(J, i, 1, ti*x4*exp(-ti*x2));
      gsl_matrix_set(J, i, 2, exp(-ti*x1));
      gsl_matrix_set(J, i, 3, -exp(-ti*x2));
      gsl_matrix_set(J, i, 4, -ti*x6*exp(-ti*x5));
      gsl_matrix_set(J, i, 5, exp(-ti*x5));
    }

  return GSL_SUCCESS;
}

static gsl_multifit_function_fdf biggs_func =
{
  &biggs_f,
  &biggs_df,
  NULL,
  biggs_N,
  biggs_P,
  NULL,
  0,
  0
};

static test_fdf_problem biggs_problem =
{
  "biggs",
  biggs_x0,
  biggs_x,
  &biggs_sumsq,
  NULL,
  &biggs_epsrel,
  &biggs_func
};
