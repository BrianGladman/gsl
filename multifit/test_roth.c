#define roth_N         2
#define roth_P         2

static double roth_x0[roth_P] = { 0.5, -2.0 };
static double roth_x[roth_P] = { 11.4127789869021, -0.896805253274477 };

static double roth_sumsq = 48.9842536792400;

static double roth_epsrel = 1.0e-9;

static int
roth_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);

  gsl_vector_set(f, 0, x1 - x2*(2.0 - x2*(5.0 - x2)) - 13.0);
  gsl_vector_set(f, 1, x1 - x2*(14.0 - x2*(1.0 + x2)) - 29.0);

  return GSL_SUCCESS;
}

static int
roth_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  double x2 = gsl_vector_get(x, 1);

  gsl_matrix_set(J, 0, 0, 1.0);
  gsl_matrix_set(J, 0, 1, -2.0 + x2*(10.0 - 3.0*x2));
  gsl_matrix_set(J, 1, 0, 1.0);
  gsl_matrix_set(J, 1, 1, -14.0 + x2*(2.0 + 3.0*x2));

  return GSL_SUCCESS;
}

static int
roth_fdf (const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * J)
{
  roth_f (x, params, f);
  roth_df (x, params, J);

  return GSL_SUCCESS;
}

static gsl_multifit_function_fdf roth_func =
{
  &roth_f,
  &roth_df,
  &roth_fdf,
  roth_N,
  roth_P,
  NULL,
  0,
  0
};

static test_fdf_problem roth_problem =
{
  "roth_freudenstein",
  roth_x0,
  roth_x,
  &roth_sumsq,
  &roth_epsrel,
  &roth_func
};
