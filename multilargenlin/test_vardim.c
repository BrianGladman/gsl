#define vardim_N         7 /* p + 2 */
#define vardim_P         5

#define vardim_NTRIES    4

static double vardim_x0[vardim_P] = { 0.8, 0.6, 0.4, 0.2, 0.0 };
static double vardim_epsrel = 1.0e-12;

static double vardim_f[vardim_N];
static double vardim_J[vardim_N * vardim_P];

static void
vardim_checksol(const double x[], const double sumsq,
                const double epsrel, const char *sname,
                const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < vardim_P; ++i)
    {
      gsl_test_rel(x[i], 1.0, epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
vardim_fdf (const gsl_vector * x, gsl_matrix * JTJ,
            gsl_vector * JTf, double * normf, void *params)
{
  gsl_matrix_view J = gsl_matrix_view_array(vardim_J, vardim_N, vardim_P);
  gsl_vector_view f = gsl_vector_view_array(vardim_f, vardim_N);
  size_t i;
  double sum = 0.0;
  gsl_matrix_view m = gsl_matrix_submatrix(&J.matrix, 0, 0, vardim_P, vardim_P);

  gsl_matrix_set_identity(&m.matrix);

  for (i = 0; i < vardim_P; ++i)
    {
      double xi = gsl_vector_get(x, i);
      sum += (i + 1.0) * (xi - 1.0);
      gsl_vector_set(&f.vector, i, xi - 1.0);
    }

  gsl_vector_set(&f.vector, vardim_P, sum);
  gsl_vector_set(&f.vector, vardim_P + 1, sum*sum);

  if (JTJ)
    {
      for (i = 0; i < vardim_P; ++i)
        {
          gsl_matrix_set(&J.matrix, vardim_P, i, i + 1.0);
          gsl_matrix_set(&J.matrix, vardim_P + 1, i, 2*(i + 1.0)*sum);
        }
    }

  *normf = gsl_blas_dnrm2(&f.vector);

  if (JTJ)
    {
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);
      gsl_blas_dgemv(CblasTrans, 1.0, &J.matrix, &f.vector, 0.0, JTf);
    }

  return GSL_SUCCESS;
}

static gsl_multilarge_function_fdf vardim_func =
{
  &vardim_fdf,
  vardim_P,
  NULL,
  0,
  0
};

static test_fdf_problem vardim_problem =
{
  "vardim",
  vardim_x0,
  NULL,
  &vardim_epsrel,
  vardim_NTRIES,
  &vardim_checksol,
  &vardim_func
};
