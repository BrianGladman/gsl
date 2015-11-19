#define biggs_N         6  /* >= p */
#define biggs_P         6

#define biggs_NTRIES    2

static double biggs_x0[biggs_P] = { 1.0, 2.0, 1.0, 1.0, 1.0, 1.0 };
static double biggs_epsrel = 1.0e-9;

static double biggs_f[biggs_N];
static double biggs_J[biggs_N * biggs_P];

static void
biggs_checksol(const double x[], const double sumsq,
               const double epsrel, const char *sname,
               const char *pname)
{
  const double sumsq_exact = 0.0;
  const double biggs_x[biggs_P] = { 1.0, 10.0, 1.0, 5.0, 4.0, 3.0 };
  const double norm_exact = 12.3288280059380;
  gsl_vector_const_view v = gsl_vector_const_view_array(biggs_x, biggs_P);
  double norm = gsl_blas_dnrm2(&v.vector);

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  /*
   * the solution vector is not unique due to permutations, so test
   * the norm instead of individual elements
   */
  gsl_test_rel(norm, norm_exact, epsrel, "%s/%s norm",
               sname, pname);
}

static int
biggs_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(biggs_J, biggs_N, biggs_P);
  gsl_vector_view f = gsl_vector_view_array(biggs_f, biggs_N);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  double x5 = gsl_vector_get(x, 4);
  double x6 = gsl_vector_get(x, 5);
  size_t i;

  for (i = 0; i < biggs_N; ++i)
    {
      double ti = 0.1 * (i + 1.0);
      double yi = exp(-ti) - 5*exp(-10*ti) + 3*exp(-4*ti);
      double fi = x3*exp(-ti*x1) - x4*exp(-ti*x2) + x6*exp(-ti*x5) - yi;

      gsl_vector_set(&f.vector, i, fi);

      if (evaldf)
        {
          gsl_matrix_set(&J.matrix, i, 0, -ti*x3*exp(-ti*x1));
          gsl_matrix_set(&J.matrix, i, 1, ti*x4*exp(-ti*x2));
          gsl_matrix_set(&J.matrix, i, 2, exp(-ti*x1));
          gsl_matrix_set(&J.matrix, i, 3, -exp(-ti*x2));
          gsl_matrix_set(&J.matrix, i, 4, -ti*x6*exp(-ti*x5));
          gsl_matrix_set(&J.matrix, i, 5, exp(-ti*x5));
        }
    }

  status = test_accumulate(2, &J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf biggs_func =
{
  &biggs_fdf,
  biggs_P,
  NULL,
  0,
  0
};

static test_fdf_problem biggs_problem =
{
  "biggs",
  biggs_x0,
  NULL,
  &biggs_epsrel,
  biggs_NTRIES,
  &biggs_checksol,
  &biggs_func
};
