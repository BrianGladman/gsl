/*
 * testnonsymm.c
 * Patrick Alken
 *
 * Compile: gcc -g -O2 -Wall -o testnonsymm testnonsymm.c -lm -lgsl -lcblas -latlas
 *
 * Usage: testnonsymm [options]
 *
 * -i             : incremental matrices
 * -b             : balance the matrices
 * -z             : compute Schur vectors and test them
 * -n size        : size of matrices
 * -l lower-bound : lower bound for matrix elements
 * -u upper-bound : upper bound for matrix elements
 * -c num         : number of matrices to solve
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <getopt.h>
#include <sys/times.h>

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_test.h>

typedef struct
{
  gsl_eigen_nonsymm_workspace *nonsymm_p;
  gsl_eigen_nonsymmv_workspace *nonsymmv_p;
  gsl_matrix *A;
  gsl_matrix *Av;
  gsl_vector_complex *eval;
  gsl_vector_complex *evalv;
  gsl_matrix_complex *evec;

  int compute_z;
  gsl_matrix *Z;
  gsl_matrix *Zv;

  size_t n_evals;
} nonsymm_workspace;

nonsymm_workspace *nonsymm_alloc(size_t n, int compute_z, int do_balance);
void nonsymm_free(nonsymm_workspace *w);
int nonsymm_proc(nonsymm_workspace *w);

/*
 * Global variables
 */
unsigned long count = 0;

/*
 * Prototypes
 */

void make_random_matrix(gsl_matrix *m, gsl_rng *r, int lower,
                        int upper);
void make_start_matrix(gsl_matrix *m, int lower);
int inc_matrix (gsl_matrix *m, int lower, int upper);
void output_matrix(gsl_matrix *m);
void print_octave(gsl_matrix *m, const char *str);
void print_matrix(gsl_matrix *m, const char *str);
void print_hess(gsl_matrix *H, const char *str);
void print_vector(gsl_vector_complex *eval, const char *str);
int cmp(double a, double b);
int compare(const void *a, const void *b);
void sort_complex_vector(gsl_vector_complex *v);
int test_evals(gsl_vector_complex *obs, gsl_vector_complex *expected,
               gsl_matrix *A, const char *obsname, const char *expname);
int test_Z(gsl_matrix *A, gsl_matrix *Z, gsl_matrix *T);
int test_eigenvectors(gsl_matrix *A, gsl_vector_complex *eval,
                      gsl_matrix_complex *evec);
void my_error_handler(const char *reason, const char *file, int line,
                      int err);

nonsymm_workspace *
nonsymm_alloc(size_t n, int compute_z, int do_balance)

{
  nonsymm_workspace *w;

  w = (nonsymm_workspace *) malloc(sizeof(nonsymm_workspace));

  memset(w, '\0', sizeof(nonsymm_workspace));

  w->nonsymm_p = gsl_eigen_nonsymm_alloc(n);
  w->nonsymmv_p = gsl_eigen_nonsymmv_alloc(n);

  w->A = gsl_matrix_alloc(n, n);
  w->Av = gsl_matrix_alloc(n, n);

  if (compute_z)
    {
      w->Z = gsl_matrix_alloc(n, n);
      w->Zv = gsl_matrix_alloc(n, n);
      w->compute_z = 1;
      gsl_eigen_nonsymm_params(1, do_balance, w->nonsymm_p);
    }
  else
    gsl_eigen_nonsymm_params(0, do_balance, w->nonsymm_p);

  w->eval = gsl_vector_complex_alloc(n);
  w->evalv = gsl_vector_complex_alloc(n);

  w->evec = gsl_matrix_complex_alloc(n, n);

  return (w);
} /* nonsymm_alloc() */

void
nonsymm_free(nonsymm_workspace *w)

{
  if (w->nonsymm_p)
    gsl_eigen_nonsymm_free(w->nonsymm_p);

  if (w->nonsymmv_p)
    gsl_eigen_nonsymmv_free(w->nonsymmv_p);

  if (w->A)
    gsl_matrix_free(w->A);

  if (w->Av)
    gsl_matrix_free(w->Av);

  if (w->Z)
    gsl_matrix_free(w->Z);

  if (w->Zv)
    gsl_matrix_free(w->Zv);

  if (w->eval)
    gsl_vector_complex_free(w->eval);

  if (w->evalv)
    gsl_vector_complex_free(w->evalv);

  if (w->evec)
    gsl_matrix_complex_free(w->evec);

  free(w);
} /* nonsymm_free() */

int
nonsymm_proc(nonsymm_workspace *w)

{
  int s1, s2, s;

  if (w->compute_z)
    {
      s1 = gsl_eigen_nonsymm_Z(w->A, w->eval, w->Z, w->nonsymm_p);
      s2 = gsl_eigen_nonsymmv_Z(w->Av, w->evalv, w->evec, w->Zv, w->nonsymmv_p);
    }
  else
    {
      s1 = gsl_eigen_nonsymm(w->A, w->eval, w->nonsymm_p);
      s2 = gsl_eigen_nonsymmv(w->Av, w->evalv, w->evec, w->nonsymmv_p);
    }

  w->n_evals = w->nonsymm_p->n_evals;

  s = 0;
  if (s1)
    s = s1;
  else if (s2)
    s = s2;

  return s;
} /* nonsymm_proc() */

/**********************************************
 * General routines
 **********************************************/

void
make_random_matrix(gsl_matrix *m, gsl_rng *r, int lower, int upper)

{
  size_t i, j;
  size_t N = m->size1;

  for (i = 0; i < N; ++i)
  {
    for (j = 0; j < N; ++j)
    {
      gsl_matrix_set(m,
                     i,
                     j,
                     gsl_rng_uniform(r) * (upper - lower) + lower);
    }
  }
} /* make_random_matrix() */

void
make_start_matrix(gsl_matrix *m, int lower)

{
  size_t i, j;
  size_t N = m->size1;

  for (i = 0; i < N; ++i)
    for (j = 0; j < N; ++j)
      gsl_matrix_set(m, i, j, (double)lower);
} /* make_start_matrix() */

int
inc_matrix (gsl_matrix *m, int lower, int upper)
{
  size_t i = 0;
  size_t N = m->size1 * m->size2;
  int carry = 1;

  for (i = 0; carry > 0 && i < N; i++)
    {
      double v = m->data[i] + carry;
      carry = (v > upper) ? 1 : 0;
      m->data[i] = (v > upper) ? lower : v;
    }

  return carry;
} /* inc_matrix() */

void
output_matrix(gsl_matrix *m)
{
  size_t i, j;
  size_t N = m->size1;
  size_t M = m->size2;

  for (i = 0; i < N; ++i)
  {
    for (j = 0; j < M; ++j)
      {
        /*printf("%10.9f%s",*/
        printf("%10.18e%s",
               gsl_matrix_get(m, i, j),
               (j < M - 1) ? "," : ";\n");
      }
  }
}

void
print_octave(gsl_matrix *m, const char *str)
{
  FILE *fp;
  size_t i, j;
  const size_t N = m->size1;
  const size_t M = m->size2;

  fp = fopen(str, "w");

  if (!fp)
    return;

  fprintf(fp, "# Created by Octave 2.1.73, Tue Aug 01 15:00:27 2006 MDT <blah@blah>\n");
  fprintf(fp, "# name: %s\n", str);
  fprintf(fp, "# type: matrix\n");
  fprintf(fp, "# rows: %u\n", N);
  fprintf(fp, "# columns: %u\n", N);

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < M; ++j)
        {
          fprintf(fp,
                  "%10.9f%s",
                  gsl_matrix_get(m, i, j),
                  (j < M - 1) ? " " : "\n");
        }
    }

  fclose(fp);
}

void
print_matrix(gsl_matrix *m, const char *str)

{
  size_t i, j;
  size_t N = m->size1;
  size_t M = m->size2;
  gsl_matrix_view v;
  size_t rows, cols;
  size_t r, c;
  char buf[100];

  /*print_octave(m, str);
  return;*/

  /*rows = GSL_MIN(15, N);*/
  rows = N;
  cols = GSL_MIN(15, N);
  /*cols = N;*/

  for (i = 0; i < N; i += rows)
  {
    for (j = 0; j < M; j += cols)
    {
      r = GSL_MIN(rows, N - i);
      c = GSL_MIN(cols, N - j);

      v = gsl_matrix_submatrix(m,
                               i,
                               j,
                               r,
                               c);

      sprintf(buf, "%s(%u:%u,%u:%u)",
              str,
              i + 1,
              i + r,
              j + 1,
              j + c);

      printf("%s = [\n", buf);

      output_matrix(&v.matrix);

      printf("]\n");
    }
  }
} /* print_matrix() */

void
print_vector(gsl_vector_complex *eval, const char *str)

{
  size_t N = eval->size;
  size_t i;
  gsl_complex z;

  printf("%s = [\n", str);

  for (i = 0; i < N; ++i)
    {
      z = gsl_vector_complex_get(eval, i);
      printf("%.18e %.18e;\n", GSL_REAL(z), GSL_IMAG(z));
    }

  printf("\n]\n");
} /* print_vector() */

int
cmp(double a, double b)

{
  return ((a > b) ? 1 : ((a < b) ? -1 : 0));
} /* cmp() */

int
compare(const void *a, const void *b)

{
  const double *x = a;
  const double *y = b;
  int r1 = cmp(y[0], x[0]);
  int r2 = cmp(y[1], x[1]);

  if (fabs(x[0] - y[0]) < 1.0e-8)
    {
      /* real parts are very close to each other */
      return r2;
    }
  else
    {
      return r1 ? r1 : r2;
    }
} /* compare() */

void
sort_complex_vector(gsl_vector_complex *v)

{
  qsort(v->data, v->size, 2 * sizeof(double), &compare);
} /* sort_complex_vector() */

int
test_evals(gsl_vector_complex *obs, gsl_vector_complex *expected,
           gsl_matrix *A, const char *obsname, const char *expname)

{
  size_t N = expected->size;
  size_t i, k;
  double max, max_abserr, max_relerr;

  max = 0.0;
  max_abserr = 0.0;
  max_relerr = 0.0;
  k = 0;

  for (i = 0; i < N; ++i)
    {
      gsl_complex z = gsl_vector_complex_get(expected, i);
      max = GSL_MAX_DBL(max, gsl_complex_abs(z));
    }

  for (i = 0; i < N; ++i)
    {
      gsl_complex z_obs = gsl_vector_complex_get(obs, i);
      gsl_complex z_exp = gsl_vector_complex_get(expected, i);

      double x_obs = GSL_REAL(z_obs);
      double y_obs = GSL_IMAG(z_obs);
      double x_exp = GSL_REAL(z_exp);
      double y_exp = GSL_IMAG(z_exp);

      double abserr_x = fabs(x_obs - x_exp);
      double abserr_y = fabs(y_obs - y_exp);
      double noise = max * GSL_DBL_EPSILON * N * N;

      max_abserr = GSL_MAX_DBL(max_abserr, abserr_x + abserr_y);

      if (abserr_x < noise && abserr_y < noise)
        continue;

      if (abserr_x > 1.0e-6 || abserr_y > 1.0e-6)
        ++k;
    }

    if (k)
      {
        printf("==== CASE %lu ===========================\n\n", count);

        print_matrix(A, "A");

        printf("=== eval - %s ===\n", expname);
        print_vector(expected, expname);

        printf("=== eval - %s ===\n", obsname);
        print_vector(obs, obsname);

        printf("max abserr = %g  max relerr = %g\n", max_abserr, max_relerr);

        printf("=========================================\n\n");
      }

    return k;
} /* test_evals() */

/* test if A = ZTZ^t (or AZ = ZT) */

int
test_Z(gsl_matrix *A, gsl_matrix *Z, gsl_matrix *T)

{
  size_t N = A->size1;
  size_t i, j, k;
  double rhs, lhs;
  double abserr;
  gsl_matrix *T1, *T2;

  T1 = gsl_matrix_alloc(N, N);
  T2 = gsl_matrix_alloc(N, N);

  /* zero the lower triangle of T */
  if (N > 3)
    {
      for (i = 2; i < N; ++i)
        {
          for (j = 0; j < (i - 1); ++j)
            {
              gsl_matrix_set(T, i, j, 0.0);
            }
        }
    }
  else if (N == 3)
    {
      gsl_matrix_set(T, 2, 0, 0.0);
    }

  /* compute T1 = A Z */
  gsl_blas_dgemm(CblasNoTrans,
                 CblasNoTrans,
                 1.0,
                 A,
                 Z,
                 0.0,
                 T1);

  /* compute T2 = Z T */
  gsl_blas_dgemm(CblasNoTrans,
                 CblasNoTrans,
                 1.0,
                 Z,
                 T,
                 0.0,
                 T2);

  k = 0;
  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          lhs = gsl_matrix_get(T1, i, j);
          rhs = gsl_matrix_get(T2, i, j);

          abserr = fabs(lhs - rhs);

          if (abserr > 1.0e-6)
            ++k;
        }
    }

  if (k)
    {
      printf("==== CASE %lu ===========================\n\n", count);

      print_matrix(A, "A");

      printf("=== Schur Form matrix ===\n");
      print_matrix(T, "T");

      printf("=== Similarity matrix ===\n");
      print_matrix(Z, "Z");

      printf("=== A Z ===\n");
      print_matrix(T1, "A Z");

      printf("=== Z T ===\n");
      print_matrix(T2, "Z T");

      printf("=== A Z - Z T ===\n");
      gsl_matrix_sub(T1, T2);
      print_matrix(T1, "A Z - Z T");

      printf("=========================================\n\n");
    }

  gsl_matrix_free(T1);
  gsl_matrix_free(T2);

  return k;
} /* test_Z() */

int
test_eigenvectors(gsl_matrix *A, gsl_vector_complex *eval,
                  gsl_matrix_complex *evec)
{
  const size_t N = A->size1;
  size_t i, j;
  int k, s;
  gsl_matrix_complex *m;
  gsl_vector_complex *x, *y;
  gsl_complex z_one, z_zero;

  m = gsl_matrix_complex_alloc(N, N);
  y = gsl_vector_complex_alloc(N);
  x = gsl_vector_complex_alloc(N);

  /* m <- A */
  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex z;

          GSL_SET_COMPLEX(&z, gsl_matrix_get(A, i, j), 0.0);
          gsl_matrix_complex_set(m, i, j, z);
        }
    }

  GSL_SET_COMPLEX(&z_one, 1.0, 0.0);
  GSL_SET_COMPLEX(&z_zero, 0.0, 0.0);

  s = 0;

  /* check eigenvalues */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_complex_view vi = gsl_matrix_complex_column(evec, i);
      gsl_complex ei = gsl_vector_complex_get(eval, i);
      double norm = gsl_blas_dznrm2(&vi.vector);

      /* check that eigenvector is normalized */
      gsl_test_rel(norm, 1.0, N * GSL_DBL_EPSILON, "case %u, normalized",
                   count);

      /* compute x = lambda * v */
      gsl_vector_complex_memcpy(x, &vi.vector);
      gsl_blas_zscal(ei, x);

      /* compute y = A v */
      gsl_blas_zgemv(CblasNoTrans, z_one, m, &vi.vector, z_zero, y);

      k = 0;

      /* now test if y = x */
      for (j = 0; j < N; ++j)
        {
          gsl_complex z;
          double lhs_r, lhs_i;
          double rhs_r, rhs_i;

          z = gsl_vector_complex_get(y, j);
          lhs_r = GSL_REAL(z);
          lhs_i = GSL_IMAG(z);

          z = gsl_vector_complex_get(x, j);
          rhs_r = GSL_REAL(z);
          rhs_i = GSL_IMAG(z);

          if (fabs(lhs_r - rhs_r) > 1e8 * GSL_DBL_EPSILON)
            ++k;
          if (fabs(lhs_i - rhs_i) > 1e8 * GSL_DBL_EPSILON)
            ++k;
        }

      if (k)
        {
          s++;

          printf("==== CASE %lu ===========================\n\n", count);

          print_matrix(A, "A");

          printf("eval = %f + %fi\n", GSL_REAL(ei), GSL_IMAG(ei));

          print_vector(&vi.vector, "v");

          print_vector(y, "Av");

          print_vector(x, "ev_v");

          printf("=========================================\n\n");
        }
    }

  gsl_matrix_complex_free(m);
  gsl_vector_complex_free(y);
  gsl_vector_complex_free(x);

  return s;
} /* test_eigenvectors() */

void
my_error_handler(const char *reason, const char *file, int line,
                 int err)

{
  printf("[caught: %s:%d: errno=%d %s]\n", file, line, err, reason);
} /* my_error_handler() */

int
main(int argc, char *argv[])

{
  nonsymm_workspace *nonsymm_workspace_p;
  size_t N;
  int c;
  gsl_matrix *A;
  gsl_rng *r;
  int incremental; /* incremental/random matrices */
  int lower;  /* lower bound */
  int upper;  /* upper bound */
  unsigned int nmat;   /* number of matrices to solve */
  int status;
  int compute_z;
  int do_balance;

  gsl_ieee_env_setup();
  gsl_rng_env_setup();

  /*gsl_set_error_handler(&my_error_handler);*/

  N = 30;
  incremental = 0;
  compute_z = 0;
  do_balance = 0;
  lower = -10;
  upper = 10;
  nmat = 0;

  while ((c = getopt(argc, argv, "izbc:n:l:u:t:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            incremental = 1;
            break;

          case 'n':
            N = strtol(optarg, NULL, 0);
            break;

          case 'l':
            lower = strtol(optarg, NULL, 0);
            break;

          case 'u':
            upper = strtol(optarg, NULL, 0);
            break;

          case 'c':
            nmat = strtoul(optarg, NULL, 0);
            break;

          case 'b':
            do_balance = 1;
            break;

          case 'z':
            compute_z = 1;
            break;

          case '?':
          default:
            printf("usage: %s [-i] [-b] [-z] [-n size] [-l lower-bound] [-u upper-bound] [-c num]\n", argv[0]);
            exit(1);
            break;
        } /* switch (c) */
    }

  A = gsl_matrix_alloc(N, N);
  nonsymm_workspace_p = nonsymm_alloc(N, compute_z, do_balance);

  if (!incremental)
    r = gsl_rng_alloc(gsl_rng_default);
  else
  {
    r = 0;
    make_start_matrix(A, lower);
  }

  fprintf(stderr, "testing N = %d", N);
  if (incremental)
    fprintf(stderr, " incrementally");
  else
    fprintf(stderr, " randomly");

  fprintf(stderr, " on element range [%d, %d]", lower, upper);

  if (compute_z)
    fprintf(stderr, ", with Schur vectors");

  if (do_balance)
    fprintf(stderr, ", with balancing");

  fprintf(stderr, "\n");

  while (1)
    {
      if (nmat && (count >= nmat))
        break;

      ++count;

      if (!incremental)
        make_random_matrix(A, r, lower, upper);
      else
        {
          status = inc_matrix(A, lower, upper);
          if (status)
            break; /* all done */
        }

      /* make copies of the matrix */
      gsl_matrix_memcpy(nonsymm_workspace_p->A, A);
      gsl_matrix_memcpy(nonsymm_workspace_p->Av, A);

      status = nonsymm_proc(nonsymm_workspace_p);

      if (status)
        {
          printf("=========== CASE %lu ============\n", count);
          printf("Failed to converge: found %u eigenvalues\n",
                 nonsymm_workspace_p->n_evals);
          print_matrix(A, "A");
        }

      gsl_eigen_nonsymmv_sort(nonsymm_workspace_p->evalv,
                              nonsymm_workspace_p->evec,
                              GSL_EIGEN_SORT_ABS_ASC);

      status = test_eigenvectors(A,
                                 nonsymm_workspace_p->evalv,
                                 nonsymm_workspace_p->evec);

      sort_complex_vector(nonsymm_workspace_p->eval);
      sort_complex_vector(nonsymm_workspace_p->evalv);

      status = test_evals(nonsymm_workspace_p->eval,
                          nonsymm_workspace_p->evalv,
                          A,
                          "nonsymm",
                          "nonsymmv");

      if (compute_z)
        {
          test_Z(A, nonsymm_workspace_p->Z, nonsymm_workspace_p->A);
          test_Z(A, nonsymm_workspace_p->Zv, nonsymm_workspace_p->Av);
        }
    }

  gsl_matrix_free(A);
  nonsymm_free(nonsymm_workspace_p);

  if (r)
    gsl_rng_free(r);

  return 0;
} /* main() */
