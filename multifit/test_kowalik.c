#define kowalik_N         11
#define kowalik_P         4

static double kowalik_x0[kowalik_P] = { 0.25, 0.39, 0.415, 0.39 };
static double kowalik_x[kowalik_P] = { 1.928069345723978e-01, 1.912823290344599e-01,
                                       1.230565070690708e-01, 1.360623308065148e-01 };

static double kowalik_sumsq = 3.075056038492370e-04;
static double kowalik_epsrel = 1.0e-8;

static double kowalik_Y[kowalik_N] = {
0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627,
0.0456, 0.0342, 0.0323, 0.0235, 0.0246
};

static double kowalik_U[kowalik_N] = {
4.0000, 2.0000, 1.0000, 0.5000, 0.2500, 0.1670,
0.1250, 0.1000, 0.0833, 0.0714, 0.0625
};

static int
kowalik_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  size_t i;

  for (i = 0; i < kowalik_N; ++i)
    {
      double yi = kowalik_Y[i];
      double ui = kowalik_U[i];
      double fi = yi - (x1*ui*(ui+x2)) / (x4 + ui*(ui + x3));
      gsl_vector_set(f, i, fi);
    }

  return GSL_SUCCESS;
}

static int
kowalik_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  size_t i;

  for (i = 0; i < kowalik_N; ++i)
    {
      double ui = kowalik_U[i];
      double term1 = ui*(ui + x2);
      double term2 = ui*(ui + x3) + x4;

      gsl_matrix_set(J, i, 0, -term1 / term2);
      gsl_matrix_set(J, i, 1, -ui*x1/term2);
      gsl_matrix_set(J, i, 2, ui*term1*x1 / (term2*term2));
      gsl_matrix_set(J, i, 3, term1*x1 / (term2*term2));
    }

  return GSL_SUCCESS;
}

static gsl_multifit_function_fdf kowalik_func =
{
  &kowalik_f,
  &kowalik_df,
  NULL,
  kowalik_N,
  kowalik_P,
  NULL,
  0,
  0
};

static test_fdf_problem kowalik_problem =
{
  "kowalik",
  kowalik_x0,
  kowalik_x,
  &kowalik_sumsq,
  NULL,
  &kowalik_epsrel,
  &kowalik_func
};
