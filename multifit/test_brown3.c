#define brown3_N         3
#define brown3_P         2

static double brown3_x0[brown3_P] = { 1.0, 1.0 };
static double brown3_x[brown3_P] = { 1.0e6, 2.0e-6 };

static double brown3_sumsq = 0.0;
static double brown3_epsrel = 1.0e-12;

static int
brown3_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);

  gsl_vector_set(f, 0, x1 - 1.0e6);
  gsl_vector_set(f, 1, x2 - 2.0e-6);
  gsl_vector_set(f, 2, x1*x2 - 2.0);

  return GSL_SUCCESS;
}

static int
brown3_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);

  gsl_matrix_set_zero(J);

  gsl_matrix_set(J, 0, 0, 1.0);
  gsl_matrix_set(J, 1, 1, 1.0);
  gsl_matrix_set(J, 2, 0, x2);
  gsl_matrix_set(J, 2, 1, x1);

  return GSL_SUCCESS;
}

static gsl_multifit_function_fdf brown3_func =
{
  &brown3_f,
  &brown3_df,
  NULL,
  brown3_N,
  brown3_P,
  NULL,
  0,
  0
};

static test_fdf_problem brown3_problem =
{
  "brown_badly_scaled",
  brown3_x0,
  brown3_x,
  &brown3_sumsq,
  NULL,
  &brown3_epsrel,
  &brown3_func
};
