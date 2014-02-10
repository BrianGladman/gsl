#define powell1_N     4
#define powell1_P     4

static double powell1_x0[powell1_P] = { 3.0, -1.0, 0.0, 1.0 };
static double powell1_x[powell1_P] = { 0.0, 0.0, 0.0, 0.0 };

static double powell1_sumsq = 0.0;

static double powell1_epsrel = 1.0e-5;

static int
powell1_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get (x, 0);
  double x2 = gsl_vector_get (x, 1);
  double x3 = gsl_vector_get (x, 2);
  double x4 = gsl_vector_get (x, 3);

  gsl_vector_set(f, 0, x1 + 10.0*x2);
  gsl_vector_set(f, 1, sqrt(5.0) * (x3 - x4));
  gsl_vector_set(f, 2, pow(x2 - 2.0*x3, 2.0));
  gsl_vector_set(f, 3, sqrt(10.0) * pow((x1 - x4), 2.0));

  return GSL_SUCCESS;
}

static int
powell1_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  double x1 = gsl_vector_get (x, 0);
  double x2 = gsl_vector_get (x, 1);
  double x3 = gsl_vector_get (x, 2);
  double x4 = gsl_vector_get (x, 3);
  double term1 = x2 - 2.0*x3;
  double term2 = x1 - x4;

  gsl_matrix_set(J, 0, 0, 1.0);
  gsl_matrix_set(J, 0, 1, 10.0);
  gsl_matrix_set(J, 0, 2, 0.0);
  gsl_matrix_set(J, 0, 3, 0.0);

  gsl_matrix_set(J, 1, 0, 0.0);
  gsl_matrix_set(J, 1, 1, 0.0);
  gsl_matrix_set(J, 1, 2, sqrt(5.0));
  gsl_matrix_set(J, 1, 3, -sqrt(5.0));

  gsl_matrix_set(J, 2, 0, 0.0);
  gsl_matrix_set(J, 2, 1, 2.0*term1);
  gsl_matrix_set(J, 2, 2, -4.0*term1);
  gsl_matrix_set(J, 2, 3, 0.0);

  gsl_matrix_set(J, 3, 0, 2.0*sqrt(10.0)*term2);
  gsl_matrix_set(J, 3, 1, 0.0);
  gsl_matrix_set(J, 3, 2, 0.0);
  gsl_matrix_set(J, 3, 3, -2.0*sqrt(10.0)*term2);

  return GSL_SUCCESS;
}

static int
powell1_fdf (const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * J)
{
  powell1_f (x, params, f);
  powell1_df (x, params, J);

  return GSL_SUCCESS;
}

static gsl_multifit_function_fdf powell1_func =
{
  &powell1_f,
  &powell1_df,
  &powell1_fdf,
  powell1_N,
  powell1_P,
  NULL,
  0,
  0
};

static test_fdf_problem powell1_problem =
{
  "powell1",
  powell1_x0,
  powell1_x,
  &powell1_sumsq,
  NULL,
  &powell1_epsrel,
  &powell1_func
};
