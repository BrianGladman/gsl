/* compile with gcc -O2 -Wall -g testunsymm.c -I .. .libs/libgsleigen.a ../linalg/.libs/libgsllinalg.a -lgsl -lgslcblas -llapack -lg2c -lm */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define FLIP_FOR_LAPACK 1
/* #define BALANCE 1  */

void alken_free (void);
void alken_setup (const gsl_matrix * m);
int alken_unsymm (gsl_matrix * m, gsl_vector_complex * eval);
void dump_eval (gsl_vector_complex * eval);
int lapack_dgees (gsl_matrix * m, gsl_vector_complex * eval);
void lapack_free (void);
void lapack_setup (const gsl_matrix * m);
void make_random_matrix (gsl_matrix * m);
void test_eigen_unsymm (const char *desc, const gsl_matrix * m);
void sort_complex_vector (gsl_vector_complex * v);
int test_evals (gsl_vector_complex * e1, gsl_vector_complex * e);
void my_error_handler (const char *reason, const char *file, int line,
                       int err);
void make_start_matrix (gsl_matrix * m);
int inc_matrix (gsl_matrix * m);

void dgees_ (char *jobvs, char *sort, int *select, int *n,
             double *a, int *lda, int *sdim, double *wr, double *wi,
             double *vs, int *ldvs, double *work, int *lwork,
             int *bwork, int *info);

void dgebal_ (char *job, int * n, double *a, int *lda, int *ilo, int *ihi,
              double *scale, int *info);

gsl_vector_complex *eval1, *eval2;
gsl_matrix *A, *A1, *A2;
gsl_vector *wr, *wi, *work;
gsl_eigen_unsymm_workspace *w;
gsl_rng *r;

int lower = -10;
int upper = +10;
int incremental = 0;
unsigned long count = 0;

int
main (int argc, char *argv[])
{
  int N = 10, status, c;

  gsl_ieee_env_setup ();
  gsl_rng_env_setup() ;

  r = gsl_rng_alloc (gsl_rng_default);

  gsl_set_error_handler (&my_error_handler);

  while ((c = getopt (argc, argv, "in:l:u:")) != -1)
    {
      switch (c)
        {
        case 'i':
          incremental = 1;
          break;

        case 'n':
          N = strtol (optarg, NULL, 0);
          break;

        case 'l':
          lower = strtol (optarg, NULL, 0);
          break;

        case 'u':
          upper = strtol (optarg, NULL, 0);
          break;

        case '?':
          fprintf (stderr, "unrecognized options\n");
        default:
          printf
            ("usage: test [-i] [-n SIZE] [-l LOWER-BOUND] [-u UPPER-BOUND]\n");
          printf
            ("defaults to random testing, to get exhaustive testing over range use -i\n");
          abort ();
        }
    }

  A = gsl_matrix_alloc (N, N);

  alken_setup (A);
  lapack_setup (A);

  printf ("testing N=%d", N);
  if (incremental)
    printf (" incrementally");
  else
    printf (" randomly");

  printf (" on element range [%d,%d]\n", lower, upper);

  if (incremental)
    make_start_matrix (A);

  while (1)
    {
      count++;

      if (!incremental)
        make_random_matrix (A);

#ifdef FLIP_FOR_LAPACK 
      gsl_matrix_transpose_memcpy (A2, A);
#else
      gsl_matrix_memcpy (A2, A);
#endif
      status = lapack_dgees (A2, eval2);
      if (status)
        {
          continue;
        };
      sort_complex_vector (eval2);

      gsl_matrix_memcpy (A1, A);
      status = alken_unsymm (A1, eval1);
      sort_complex_vector (eval1);

      test_evals (eval1, eval2);

      if (incremental)
        {
          status = inc_matrix (A);
          if (status)
            break;              /* all done */
        }
    }

  alken_free ();
  lapack_free ();
  gsl_matrix_free (A);

  gsl_rng_free (r);
  exit (0);
}

void
make_random_matrix (gsl_matrix * m)
{
  size_t i, j;
  size_t N = A->size1;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      gsl_matrix_set (m, i, j, gsl_rng_uniform (r) * (upper - lower) + lower);
}

void
make_start_matrix (gsl_matrix * m)
{
  size_t i, j;
  size_t N = A->size1;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      gsl_matrix_set (m, i, j, lower);
}

int
inc_matrix (gsl_matrix * m)
{
  size_t i = 0;
  size_t N = A->size1 * A->size2;
  int carry = 1;

  for (i = 0; carry > 0 && i < N; i++)
    {
      double v = m->data[i] + carry;
      carry = (v > upper) ? 1 : 0;
      m->data[i] = (v > upper) ? lower : v;
    }

  return carry;
}


void
alken_setup (const gsl_matrix * m)
{
  size_t N = m->size1;
  A1 = gsl_matrix_alloc (N, N);
  eval1 = gsl_vector_complex_alloc (N);
  w = gsl_eigen_unsymm_alloc (N);
}

void
alken_free (void)
{
  gsl_eigen_unsymm_free (w);
  gsl_vector_complex_free (eval1);
  gsl_matrix_free (A1);
}

int
alken_unsymm (gsl_matrix * m, gsl_vector_complex * eval)
{
  int s;
  s = gsl_eigen_unsymm (m, eval, w);
  return s;
}

void
lapack_setup (const gsl_matrix * m)
{
  int info, lwork, size, sdim, ldvs, ilo, ihi;
  size_t N = m->size1;
  char jobvs = 'N', sort = 'N', jobg = 'B';

  A2 = gsl_matrix_alloc (N, N);
  eval2 = gsl_vector_complex_alloc (N);

  wr = gsl_vector_alloc (N);
  wi = gsl_vector_alloc (N);
  work = gsl_vector_calloc (3 * N);

  info = 0;
  lwork = -1;
  sdim = 0;
  ldvs = 1;

#ifdef BALANCE
  dgebal_ ((char *)&jobg, (int *)&N, A->data, (int *)&A->tda, 
           &ilo, &ihi, wr->data, (int *)&info);
#endif

  dgees_ ((char *) &jobvs, (char *) &sort, (int *) NULL, (int *) &N,
          A->data, (int *) &(A->tda), (int *) &sdim, wr->data, wi->data,
          (double *) NULL, (int *) &ldvs, work->data, (int *) &lwork,
          (int *) NULL, (int *) &info);

  size = work->data[0];

  //printf ("lwork = %d  info = %d\n", lwork, info);
  //printf ("recommended size is %d\n", size);

  if (size > work->size)
    {
      gsl_vector_free (work);
      work = gsl_vector_calloc (size);
    }
}

void
lapack_free (void)
{
  gsl_vector_free (wr);
  gsl_vector_free (wi);
  gsl_vector_free (work);
  gsl_matrix_free (A2);
  gsl_vector_complex_free (eval2);
}

int
lapack_dgees (gsl_matrix * m, gsl_vector_complex * eval)
{
  int info, lwork, sdim, ldvs;
  size_t N = m->size1, i;
  char jobvs = 'N', sort = 'N';

  info = 0;
  lwork = work->size;
  sdim = 0;
  ldvs = 1;

  dgees_ ((char *) &jobvs, (char *) &sort, (int *) NULL, (int *) &N,
          m->data, (int *) &(m->tda), (int *) &sdim, wr->data, wi->data,
          (double *) NULL, (int *) &ldvs, work->data, (int *) &lwork,
          (int *) NULL, (int *) &info);

  if (info != 0)
    printf ("info = %d\n", info);

  for (i = 0; i < N; i++)
    {
      gsl_complex z;
      GSL_SET_COMPLEX (&z, wr->data[i], wi->data[i]);
      gsl_vector_complex_set (eval, i, z);
    }

  return info;
}

void
dump_eval (gsl_vector_complex * eval)
{
  int i;
  size_t N = eval->size;

  for (i = 0; i < N; i++)
    {
      gsl_complex z = gsl_vector_complex_get (eval, i);
      printf ("% .18e % .18e\n", GSL_REAL (z), GSL_IMAG (z));
    }
  printf ("\n");
}

inline int
cmp (double a, double b)
{
  return ((a > b) ? 1 : ((a < b) ? -1 : 0));
}

int
compare (const void *a, const void *b)
{
  const double *x = a;
  const double *y = b;
  int r1 = cmp (y[0], x[0]);
  int r2 = cmp (y[1], x[1]);
  return r1 ? r1 : r2;
}

void
sort_complex_vector (gsl_vector_complex * v)
{
  qsort (v->data, v->size, 2 * sizeof (double), &compare);
}

int
test_evals (gsl_vector_complex * obs, gsl_vector_complex * expected)
{
  int i, j, k = 0, idx = -1;
  size_t N = expected->size;
  double max = 0.0, max_abserr = 0.0, max_relerr = 0.0;

  for (i = 0; i < N; i++)
    {
      gsl_complex z = gsl_vector_complex_get (expected, i);
      max = GSL_MAX_DBL (max, gsl_complex_abs (z));
    }

  for (i = 0; i < N; i++)
    {
      gsl_complex z_obs = gsl_vector_complex_get (obs, i);
      gsl_complex z_exp = gsl_vector_complex_get (expected, i);

      double x_obs = GSL_REAL (z_obs), y_obs = GSL_IMAG (z_obs);
      double x = GSL_REAL (z_exp), y = GSL_IMAG (z_exp);

      double abserr_x = fabs (x_obs - x);
      double abserr_y = fabs (y_obs - y);
      double noise = max * GSL_DBL_EPSILON * N * N;

      max_abserr = GSL_MAX_DBL(max_abserr, abserr_x + abserr_y);

      if (abserr_x < noise && abserr_y < noise)
        continue;

      if (abserr_x > 1e-6 || abserr_y > 1e-6) {
        k++;;
      }

      {
        double relerr_x = abserr_x / fabs (x);
        double relerr_y = abserr_y / fabs (y);

        max_relerr = GSL_MAX_DBL(max_relerr, relerr_x + relerr_y);

        if (relerr_x > 1e-6 || relerr_y > 1e-6)
          {
            k++;
          }
      }
    }

  if (k)
    {
      printf ("==== CASE %lu ==========================\n\n", count);

      printf (" A = [ \n");
      for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
          {
            printf ("% 10.18g%s", gsl_matrix_get (A, i, j),
                    (j < N - 1) ? "," : ";\n");
          }
      printf ("]\n\n");

      printf ("=== eval - lapack dgees ===\n");
      dump_eval (expected);

      printf ("=== eval - alken ===\n");
      dump_eval (obs);

      printf("max abserr = %g  max relerr = %g\n", max_abserr, max_relerr);

      printf ("==================================\n\n");
    }

  return k;
}

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  printf ("[caught: %s:%d: errno=%d %s]\n", file, line, err, reason);
}
