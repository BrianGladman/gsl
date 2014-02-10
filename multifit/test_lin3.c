#define lin3_N         50  /* can be anything >= p */
#define lin3_P         10  /* >= 3 */

static double lin3_x0[lin3_P] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
/* no unique solution vector */

/* F(x*) = [m^2 + 3m - 6] / [2(2m - 3)] */
static double lin3_sumsq = 1.362886597938144e+01;
static double lin3_epsrel = 1.0e-12;

static int
lin3_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  size_t i, j;

  gsl_vector_set(f, 0, -1.0);
  gsl_vector_set(f, lin3_N - 1, -1.0);

  for (i = 1; i < lin3_N - 1; ++i)
    {
      double fi = 0.0;

      for (j = 1; j < lin3_P - 1; ++j)
        {
          double xj = gsl_vector_get(x, j);
          fi += (j + 1) * xj;
        }

      fi = i * fi - 1.0;
      gsl_vector_set(f, i, fi);
    }

  return GSL_SUCCESS;
}

static int
lin3_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  size_t i, j;

  gsl_matrix_set_zero(J);

  for (i = 1; i < lin3_N - 1; ++i)
    {
      for (j = 1; j < lin3_P - 1; ++j)
        {
          gsl_matrix_set(J, i, j, i * (j + 1.0));
        }
    }

  return GSL_SUCCESS;
}

static gsl_multifit_function_fdf lin3_func =
{
  &lin3_f,
  &lin3_df,
  NULL,
  lin3_N,
  lin3_P,
  NULL,
  0,
  0
};

static test_fdf_problem lin3_problem =
{
  "linear_rank1zeros",
  lin3_x0,
  NULL, /* no unique solution vector */
  &lin3_sumsq,
  NULL,
  &lin3_epsrel,
  &lin3_func
};
