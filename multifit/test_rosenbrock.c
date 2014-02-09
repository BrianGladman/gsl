#define rosenbrock_N         2
#define rosenbrock_P         2

static double rosenbrock_x0[rosenbrock_P] = { -1.2, 1.0 };
static double rosenbrock_x[rosenbrock_P] = { 1.0, 1.0 };

static double rosenbrock_sumsq = 0.0;

static double rosenbrock_epsrel = 1.0e-12;

static int
rosenbrock_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);

  gsl_vector_set(f, 0, 10.0 * (x2 - x1*x1));
  gsl_vector_set(f, 1, 1.0 - x1);

  return GSL_SUCCESS;
}

static int
rosenbrock_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  double x1 = gsl_vector_get(x, 0);

  gsl_matrix_set(J, 0, 0, -20.0*x1);
  gsl_matrix_set(J, 0, 1, 10.0);
  gsl_matrix_set(J, 1, 0, -1.0);
  gsl_matrix_set(J, 1, 1, 0.0);

  return GSL_SUCCESS;
}

static int
rosenbrock_fdf (const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * J)
{
  rosenbrock_f (x, params, f);
  rosenbrock_df (x, params, J);

  return GSL_SUCCESS;
}

static gsl_multifit_function_fdf rosenbrock_func =
{
  &rosenbrock_f,
  &rosenbrock_df,
  &rosenbrock_fdf,
  rosenbrock_N,
  rosenbrock_P,
  NULL,
  0,
  0
};

static test_fdf_problem rosenbrock_problem =
{
  "rosenbrock",
  rosenbrock_x0,
  rosenbrock_x,
  &rosenbrock_sumsq,
  &rosenbrock_epsrel,
  &rosenbrock_func
};
