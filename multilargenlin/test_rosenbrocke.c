#define rosenbrocke_N         8 /* = p */
#define rosenbrocke_P         8 /* must be even */

#define rosenbrocke_NTRIES    4

static double rosenbrocke_x0[rosenbrocke_P] = { -1.2, 1.0, -1.2, 1.0,
                                                -1.2, 1.0, -1.2, 1.0 };
static double rosenbrocke_epsrel = 1.0e-12;

static double rosenbrocke_f[rosenbrocke_N];
static double rosenbrocke_J[rosenbrocke_N * rosenbrocke_P];

static void
rosenbrocke_checksol(const double x[], const double sumsq,
                     const double epsrel, const char *sname,
                     const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;
  const double rosenbrocke_x[rosenbrocke_P] = { 1.0, 1.0, 1.0, 1.0,
                                                1.0, 1.0, 1.0, 1.0 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < rosenbrocke_P; ++i)
    {
      gsl_test_rel(x[i], rosenbrocke_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
rosenbrocke_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(rosenbrocke_J, rosenbrocke_N, rosenbrocke_P);
  gsl_vector_view f = gsl_vector_view_array(rosenbrocke_f, rosenbrocke_N);
  size_t i;

  gsl_matrix_set_zero(&J.matrix);

  for (i = 0; i < rosenbrocke_N / 2; ++i)
    {
      double x2i = gsl_vector_get(x, 2*i + 1);
      double x2im1 = gsl_vector_get(x, 2*i);

      gsl_vector_set(&f.vector, 2*i, 10.0 * (x2i - x2im1*x2im1));
      gsl_vector_set(&f.vector, 2*i + 1, 1.0 - x2im1);

      if (evaldf)
        {
          gsl_matrix_set(&J.matrix, 2*i, 2*i, -20.0*x2im1);
          gsl_matrix_set(&J.matrix, 2*i, 2*i + 1, 10.0);
          gsl_matrix_set(&J.matrix, 2*i + 1, 2*i, -1.0);
        }
    }

  status = gsl_multilarge_nlinear_accumulate(&J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf rosenbrocke_func =
{
  &rosenbrocke_fdf,
  rosenbrocke_P,
  NULL,
  0,
  0
};

static test_fdf_problem rosenbrocke_problem =
{
  "rosenbrock_extended",
  rosenbrocke_x0,
  NULL,
  &rosenbrocke_epsrel,
  rosenbrocke_NTRIES,
  &rosenbrocke_checksol,
  &rosenbrocke_func
};
