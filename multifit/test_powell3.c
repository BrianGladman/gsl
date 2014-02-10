#define powell3_N         2
#define powell3_P         2

static double powell3_x0[powell3_P] = { 0.0, 1.0 };
static double powell3_x[powell3_P] = { 1.09815932969975976e-05, 9.10614673986700218 };

static double powell3_sumsq = 0.0;
static double powell3_epsrel = 1.0e-12;

static int
powell3_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);

  gsl_vector_set(f, 0, 1.0e4*x1*x2 - 1.0);
  gsl_vector_set(f, 1, exp(-x1) + exp(-x2) - 1.0001);

  return GSL_SUCCESS;
}

static int
powell3_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);

  gsl_matrix_set(J, 0, 0, 1.0e4*x2);
  gsl_matrix_set(J, 0, 1, 1.0e4*x1);

  gsl_matrix_set(J, 1, 0, -exp(-x1));
  gsl_matrix_set(J, 1, 1, -exp(-x2));

  return GSL_SUCCESS;
}

static gsl_multifit_function_fdf powell3_func =
{
  &powell3_f,
  &powell3_df,
  NULL,
  powell3_N,
  powell3_P,
  NULL,
  0,
  0
};

static test_fdf_problem powell3_problem =
{
  "powell_badly_scaled",
  powell3_x0,
  powell3_x,
  &powell3_sumsq,
  NULL,
  &powell3_epsrel,
  &powell3_func
};
