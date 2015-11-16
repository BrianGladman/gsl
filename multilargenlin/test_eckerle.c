#define eckerle_N       35
#define eckerle_P       3

#define eckerle_NTRIES  1

static double eckerle_x0[eckerle_P] = { 1.0, 10.0, 500.0 };
static double eckerle_epsrel = 1.0e-7;

static double eckerle_f[eckerle_N];
static double eckerle_J[eckerle_N * eckerle_P];

static double eckerle_sigma[eckerle_P] = {
  1.5408051163E-02, 4.6803020753E-02, 4.6800518816E-02
};

static double eckerle_X[eckerle_N] = {
  400.000000, 405.000000, 410.000000, 415.000000, 420.000000,
  425.000000, 430.000000, 435.000000, 436.500000, 438.000000,
  439.500000, 441.000000, 442.500000, 444.000000, 445.500000,
  447.000000, 448.500000, 450.000000, 451.500000, 453.000000,
  454.500000, 456.000000, 457.500000, 459.000000, 460.500000,
  462.000000, 463.500000, 465.000000, 470.000000, 475.000000,
  480.000000, 485.000000, 490.000000, 495.000000, 500.000000 };

static double eckerle_F[eckerle_N] = {
  0.0001575, 0.0001699, 0.0002350, 0.0003102, 0.0004917,
  0.0008710, 0.0017418, 0.0046400, 0.0065895, 0.0097302,
  0.0149002, 0.0237310, 0.0401683, 0.0712559, 0.1264458,
  0.2073413, 0.2902366, 0.3445623, 0.3698049, 0.3668534,
  0.3106727, 0.2078154, 0.1164354, 0.0616764, 0.0337200,
  0.0194023, 0.0117831, 0.0074357, 0.0022732, 0.0008800,
  0.0004579, 0.0002345, 0.0001586, 0.0001143, 0.0000710 };

static void
eckerle_checksol(const double x[], const double sumsq,
                 const double epsrel, const char *sname,
                 const char *pname)
{
  size_t i;
  const double sumsq_exact = 1.4635887487E-03;
  const double eckerle_x[eckerle_P] = { 1.5543827178E+00,
                                        4.0888321754E+00,
                                        4.5154121844E+02 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < eckerle_P; ++i)
    {
      gsl_test_rel(x[i], eckerle_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
eckerle_fdf (const int eval_J, const gsl_vector * x, void *params, void * work)
{
  int status;
  double b[eckerle_P];
  gsl_matrix_view J = gsl_matrix_view_array(eckerle_J, eckerle_N, eckerle_P);
  gsl_vector_view f = gsl_vector_view_array(eckerle_f, eckerle_N);
  size_t i;

  for (i = 0; i < eckerle_P; i++)
    {
      b[i] = gsl_vector_get(x, i);
    }

  for (i = 0; i < eckerle_N; i++)
    {
      double xi = eckerle_X[i];
      double term1 = xi - b[2];
      double term2 = exp(-0.5 * term1 * term1 / (b[1] * b[1]));
      double yi = b[0] / b[1] * term2;

      gsl_vector_set (&f.vector, i, yi - eckerle_F[i]);

      if (eval_J)
        {
          gsl_matrix_set (&J.matrix, i, 0, term2 / b[1]);
          gsl_matrix_set (&J.matrix, i, 1,
            -b[0] * term2 / (b[1] * b[1]) + b[0] / pow(b[1], 4.0) * term2 * term1 * term1);
          gsl_matrix_set (&J.matrix, i, 2, b[0] / pow(b[1], 3.0) * term1 * term2);
        }
    }

  status = gsl_multilarge_nlinear_accumulate(&J.matrix, &f.vector, work);

  return status;
}

static gsl_multilarge_function_fdf eckerle_func =
{
  &eckerle_fdf,
  eckerle_P,
  NULL
};

static test_fdf_problem eckerle_problem =
{
  "nist-eckerle",
  eckerle_x0,
  eckerle_sigma,
  &eckerle_epsrel,
  eckerle_NTRIES,
  &eckerle_checksol,
  &eckerle_func
};
