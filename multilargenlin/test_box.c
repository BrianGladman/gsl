#define box_N         10 /* can be >= p */
#define box_P         3

#define box_NTRIES    1

static double box_x0[box_P] = { 0.0, 10.0, 20.0 };
static double box_epsrel = 1.0e-12;

static double box_f[box_N];
static double box_J[box_N * box_P];

static void
box_checksol(const double x[], const double sumsq,
             const double epsrel, const char *sname,
             const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;
  const double box_x[box_P] = { 1.0, 10.0, 1.0 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < box_P; ++i)
    {
      gsl_test_rel(x[i], box_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
box_fdf (const int evaldf, const gsl_vector * x, void *params, void * work)
{
  int status;
  gsl_matrix_view J = gsl_matrix_view_array(box_J, box_N, box_P);
  gsl_vector_view f = gsl_vector_view_array(box_f, box_N);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  size_t i;

  for (i = 0; i < box_N; ++i)
    {
      double ti = (i + 1.0) / 10.0;
      double term1 = exp(-x1*ti);
      double term2 = exp(-x2*ti);
      double term3 = exp(-10.0*ti) - exp(-ti);
      double fi = term1 - term2 + x3*term3;

      gsl_vector_set(&f.vector, i, fi);

      if (evaldf)
        {
          gsl_matrix_set(&J.matrix, i, 0, -ti*term1);
          gsl_matrix_set(&J.matrix, i, 1, ti*term2);
          gsl_matrix_set(&J.matrix, i, 2, term3);
        }
    }

  status = gsl_multilarge_nlinear_accumulate(&J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf box_func =
{
  &box_fdf,
  box_P,
  NULL,
  0,
  0
};

static test_fdf_problem box_problem =
{
  "box3d",
  box_x0,
  NULL,
  &box_epsrel,
  box_NTRIES,
  &box_checksol,
  &box_func
};
