#define lin2_N         20  /* can be anything >= p */
#define lin2_P         5

static double lin2_x0[lin2_P] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
/* no unique solution vector */

/* F(x*) = [m(m-1}] / [2(2m + 1)] */
static double lin2_sumsq = 4.634146341463414e+00;
static double lin2_epsrel = 1.0e-12;

static int
lin2_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  size_t i, j;

  for (i = 0; i < lin2_N; ++i)
    {
      double fi = 0.0;

      for (j = 0; j < lin2_P; ++j)
        {
          double xj = gsl_vector_get(x, j);
          fi += (j + 1) * xj;
        }

      fi = (i + 1) * fi - 1.0;
      gsl_vector_set(f, i, fi);
    }

  return GSL_SUCCESS;
}

static int
lin2_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  size_t i, j;

  for (i = 0; i < lin2_N; ++i)
    {
      for (j = 0; j < lin2_P; ++j)
        {
          gsl_matrix_set(J, i, j, (i + 1.0) * (j + 1.0));
        }
    }

  return GSL_SUCCESS;
}

static gsl_multifit_function_fdf lin2_func =
{
  &lin2_f,
  &lin2_df,
  NULL,
  lin2_N,
  lin2_P,
  NULL,
  0,
  0
};

static test_fdf_problem lin2_problem =
{
  "linear_rank1",
  lin2_x0,
  NULL, /* no unique solution vector */
  &lin2_sumsq,
  NULL,
  &lin2_epsrel,
  &lin2_func
};
