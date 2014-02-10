#define lin1_N         11  /* can be anything >= p */
#define lin1_P         5

static double lin1_x0[lin1_P] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
static double lin1_x[lin1_P] = { -1.0, -1.0, -1.0, -1.0, -1.0 };

static double lin1_sumsq = (double) (lin1_N - lin1_P);

static double lin1_epsrel = 1.0e-12;

static int
lin1_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  size_t i, j;

  for (i = 0; i < lin1_N; ++i)
    {
      double fi = 0.0;

      for (j = 0; j < lin1_P; ++j)
        {
          double xj = gsl_vector_get(x, j);
          double Aij = (i == j) ? 1.0 : 0.0;
          Aij -= 2.0 / lin1_N;
          fi += Aij * xj;
        }

      fi -= 1.0;
      gsl_vector_set(f, i, fi);
    }

  return GSL_SUCCESS;
}

static int
lin1_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  size_t i, j;

  for (i = 0; i < lin1_N; ++i)
    {
      for (j = 0; j < lin1_P; ++j)
        {
          double Jij = (i == j) ? 1.0 : 0.0;
          Jij -= 2.0 / lin1_N;
          gsl_matrix_set(J, i, j, Jij);
        }
    }

  return GSL_SUCCESS;
}

static gsl_multifit_function_fdf lin1_func =
{
  &lin1_f,
  &lin1_df,
  NULL,
  lin1_N,
  lin1_P,
  NULL,
  0,
  0
};

static test_fdf_problem lin1_problem =
{
  "lin1",
  lin1_x0,
  lin1_x,
  &lin1_sumsq,
  NULL,
  &lin1_epsrel,
  &lin1_func
};
