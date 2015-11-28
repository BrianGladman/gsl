#define boxbod_N       6
#define boxbod_P       2

#define boxbod_NTRIES  2

static double boxbod_x0[boxbod_P] = { 100.0, 0.75 };
static double boxbod_epsrel = 1.0e-7;

static double boxbod_sigma[boxbod_P] = {
  1.2354515176E+01, 1.0455993237E-01
};

static double boxbod_f[boxbod_N];
static double boxbod_J[boxbod_N * boxbod_P];

static double boxbod_X[boxbod_N] = { 1.0, 2.0, 3.0, 5.0, 7.0, 10.0 };

static double boxbod_F[boxbod_N] = { 109.0, 149.0, 149.0, 191.0,
                                     213.0, 224.0 };

static void
boxbod_checksol(const double x[], const double sumsq,
                 const double epsrel, const char *sname,
                 const char *pname)
{
  size_t i;
  const double sumsq_exact = 1.1680088766E+03;
  const double boxbod_x[boxbod_P] = { 2.1380940889E+02,
                                      5.4723748542E-01 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < boxbod_P; ++i)
    {
      gsl_test_rel(x[i], boxbod_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
boxbod_fdf (const gsl_vector * x, gsl_matrix * JTJ,
            gsl_vector * JTf, double * normf, void *params)
{
  double b[boxbod_P];
  gsl_matrix_view J = gsl_matrix_view_array(boxbod_J, boxbod_N, boxbod_P);
  gsl_vector_view f = gsl_vector_view_array(boxbod_f, boxbod_N);
  size_t i;

  for (i = 0; i < boxbod_P; i++)
    {
      b[i] = gsl_vector_get(x, i);
    }

  for (i = 0; i < boxbod_N; i++)
    {
      double xi = boxbod_X[i];
      double term = exp(-b[1] * xi);
      double yi = b[0] * (1.0 - term);

      gsl_vector_set (&f.vector, i, yi - boxbod_F[i]);

      if (JTJ)
        {
          gsl_matrix_set (&J.matrix, i, 0, 1.0 - term);
          gsl_matrix_set (&J.matrix, i, 1, b[0] * term * xi);
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

static gsl_multilarge_function_fdf boxbod_func =
{
  &boxbod_fdf,
  boxbod_P,
  NULL
};

static test_fdf_problem boxbod_problem =
{
  "nist-boxbod",
  boxbod_x0,
  boxbod_sigma,
  &boxbod_epsrel,
  boxbod_NTRIES,
  &boxbod_checksol,
  &boxbod_func
};
