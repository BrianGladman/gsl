/*
 * testgen.c
 * Patrick Alken
 *
 * Compile: gcc -g -O2 -Wall -o testgen testgen.c -lm -lgsl -llapack -lf77blas -lcblas -latlas -lg2c
 *
 * Usage: testgen [options]
 *
 * -i             : incremental matrices
 * -z             : compute Schur vectors and test them
 * -n size        : size of matrices
 * -l lower-bound : lower bound for matrix elements
 * -u upper-bound : upper bound for matrix elements
 * -c num         : number of matrices to solve
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

typedef struct
{
  gsl_eigen_gen_workspace *gen_p;
  gsl_matrix *A;
  gsl_matrix *B;
  gsl_vector_complex *alpha;
  gsl_vector *beta;
  gsl_vector_complex *evals;

  gsl_matrix *Q;
  gsl_matrix *Z;
  int compute_schur;

  size_t n_evals;
} gen_workspace;

gen_workspace *gen_alloc(size_t n, int compute_schur);
void gen_free(gen_workspace *w);
int gen_proc(gen_workspace *w);

typedef struct
{
  gsl_matrix *A;
  gsl_matrix *B;
  gsl_matrix *Q;
  gsl_matrix *Z;
  int N;

  char jobvsl;
  char jobvsr;
  char sort;
  int selctg;
  int lda;
  int ldb;
  int sdim;
  double *alphar;
  double *alphai;
  gsl_vector *beta;
  int ldvsr;
  int lwork;
  int info;
  double *work;
  
  gsl_vector_complex *evals;
  gsl_vector_complex *alpha;
  size_t n_evals;
} lapack_workspace;

lapack_workspace *lapack_alloc(const size_t n);
void lapack_free(lapack_workspace *w);
int lapack_proc(lapack_workspace *w);

void dgges_(char *jobvsl, char *jobvsr, char *sort, int *selctg, int *n,
            double *a, int *lda, double *b, int *ldb, int *sdim,
            double *alphar, double *alphai, double *beta, double *vsl,
            int *ldvsl, double *vsr, int *ldvsr, double *work, int *lwork,
            int *bwork, int *info);

/*
 * Global variables
 */
unsigned long count = 0;

/*
 * Prototypes
 */

void make_random_matrix(gsl_matrix *m, gsl_rng *r, int lower, int upper);
void make_random_integer_matrix(gsl_matrix *m, gsl_rng *r, int lower,
                                int upper);
void make_start_matrix(gsl_matrix *m, int lower);
int inc_matrix (gsl_matrix *m, int lower, int upper);
void output_matrix(gsl_matrix *m);
void print_matrix(gsl_matrix *m, const char *str);
int test_evals(gsl_vector_complex *obs, gsl_vector_complex *expected,
               gsl_matrix *A, gsl_matrix *B,
               const char *obsname, const char *expname);
int test_alpha(gsl_vector_complex *obs, gsl_vector_complex *expected,
               gsl_matrix *A, gsl_matrix *B, const char *obsname,
               const char *expname);
int test_beta(gsl_vector *obs, gsl_vector *expected,
              gsl_matrix *A, gsl_matrix *B, const char *obsname,
              const char *expname);
void test_schur(gsl_matrix *A, gsl_matrix *S, gsl_matrix *Q, gsl_matrix *Z);
void print_vector(gsl_vector_complex *eval, const char *str);
int cmp(double a, double b);
int compare(const void *a, const void *b);
void sort_complex_vector(gsl_vector_complex *v);

gen_workspace *
gen_alloc(size_t n, int compute_schur)
{
  gen_workspace *w;

  w = (gen_workspace *) calloc(1, sizeof(gen_workspace));

  w->gen_p = gsl_eigen_gen_alloc(n);

  w->A = gsl_matrix_alloc(n, n);
  w->B = gsl_matrix_alloc(n, n);
  w->alpha = gsl_vector_complex_alloc(n);
  w->beta = gsl_vector_alloc(n);
  w->evals = gsl_vector_complex_alloc(n);
  w->compute_schur = compute_schur;

  if (compute_schur)
    {
      w->Q = gsl_matrix_alloc(n, n);
      w->Z = gsl_matrix_alloc(n, n);
      gsl_eigen_gen_params(1, 1, 0, w->gen_p);
    }

  return (w);
} /* gen_alloc() */

void
gen_free(gen_workspace *w)
{
  if (w->gen_p)
    gsl_eigen_gen_free(w->gen_p);

  if (w->A)
    gsl_matrix_free(w->A);

  if (w->B)
    gsl_matrix_free(w->B);

  if (w->alpha)
    gsl_vector_complex_free(w->alpha);

  if (w->beta)
    gsl_vector_free(w->beta);

  if (w->evals)
    gsl_vector_complex_free(w->evals);

  if (w->Q)
    gsl_matrix_free(w->Q);

  if (w->Z)
    gsl_matrix_free(w->Z);

  free(w);
}

int
gen_proc(gen_workspace *w)
{
  int s;

  s = gsl_eigen_gen_QZ(w->A, w->B, w->alpha, w->beta, w->Q, w->Z, w->gen_p);

  w->n_evals = w->gen_p->n_evals;

  return s;
} /* gen_proc() */

lapack_workspace *
lapack_alloc(const size_t n)
{
  lapack_workspace *w;
  double work[1];

  w = (lapack_workspace *) calloc(1, sizeof(lapack_workspace));

  w->A = gsl_matrix_alloc(n, n);
  w->B = gsl_matrix_alloc(n, n);
  w->Q = gsl_matrix_alloc(n, n);
  w->Z = gsl_matrix_alloc(n, n);
  w->alphar = malloc(n * sizeof(double));
  w->alphai = malloc(n * sizeof(double));
  w->beta = gsl_vector_alloc(n);
  w->alpha = gsl_vector_complex_alloc(n);
  w->evals = gsl_vector_complex_alloc(n);

  w->N = (int) n;
  w->n_evals = 0;

  w->jobvsl = 'N';
  w->jobvsr = 'N';
  w->sort = 'N';
  w->info = 0;

  w->lwork = -1;
  dgges_(&w->jobvsl,
         &w->jobvsr,
         &w->sort,
         (int *) 0,
         &w->N,
         w->A->data,
         (int *) &w->A->tda,
         w->B->data,
         (int *) &w->B->tda,
         &w->sdim,
         w->alphar,
         w->alphai,
         w->beta->data,
         w->Q->data,
         (int *) &w->Q->tda,
         w->Z->data,
         (int *) &w->Z->tda,
         work,
         &w->lwork,
         (int *) 0,
         &w->info);

  w->lwork = (int) work[0];
  w->work = malloc(w->lwork * sizeof(double));

  return (w);
} /* lapack_alloc() */

void
lapack_free(lapack_workspace *w)
{
  if (w->A)
    gsl_matrix_free(w->A);

  if (w->B)
    gsl_matrix_free(w->B);

  if (w->Q)
    gsl_matrix_free(w->Q);

  if (w->Z)
    gsl_matrix_free(w->Z);

  if (w->work)
    free(w->work);

  if (w->alphar)
    free(w->alphar);

  if (w->alphai)
    free(w->alphai);

  if (w->beta)
    gsl_vector_free(w->beta);

  if (w->alpha)
    gsl_vector_complex_free(w->alpha);

  if (w->evals)
    gsl_vector_complex_free(w->evals);

  free(w);
} /* lapack_free() */

int
lapack_proc(lapack_workspace *w)
{
  dgges_(&w->jobvsl,
         &w->jobvsr,
         &w->sort,
         (int *) 0,
         &w->N,
         w->A->data,
         (int *) &w->A->tda,
         w->B->data,
         (int *) &w->B->tda,
         &w->sdim,
         w->alphar,
         w->alphai,
         w->beta->data,
         w->Q->data,
         (int *) &w->Q->tda,
         w->Z->data,
         (int *) &w->Z->tda,
         w->work,
         &w->lwork,
         (int *) 0,
         &w->info);

  return (w->info);
} /* lapack_proc() */

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
make_random_integer_matrix(gsl_matrix *m, gsl_rng *r, int lower, int upper)
{
  size_t i, j;
  size_t N = m->size1;

  for (i = 0; i < N; ++i)
  {
    for (j = 0; j < N; ++j)
    {
      double a = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_matrix_set(m, i, j, floor(a));
    }
  }
} /* make_random_integer_matrix() */

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
        printf("%10.18e%s",
        /*printf("%10.18e%s",*/
               gsl_matrix_get(m, i, j),
               (j < M - 1) ? "," : ";\n");
      }
  }
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

      v = gsl_matrix_submatrix(m, i, j, r, c);

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

int
test_evals(gsl_vector_complex *obs, gsl_vector_complex *expected,
           gsl_matrix *A, gsl_matrix *B, const char *obsname,
           const char *expname)
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
        print_matrix(B, "B");

        printf("=== eval - %s ===\n", expname);
        print_vector(expected, expname);

        printf("=== eval - %s ===\n", obsname);
        print_vector(obs, obsname);

        printf("max abserr = %g  max relerr = %g\n", max_abserr, max_relerr);

        printf("=========================================\n\n");
      }

    return k;
} /* test_evals() */

int
test_alpha(gsl_vector_complex *obs, gsl_vector_complex *expected,
           gsl_matrix *A, gsl_matrix *B, const char *obsname,
           const char *expname)
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
      double y_obs = fabs(GSL_IMAG(z_obs));
      double x_exp = GSL_REAL(z_exp);
      double y_exp = fabs(GSL_IMAG(z_exp));

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
        print_matrix(B, "B");

        printf("=== alpha - %s ===\n", expname);
        print_vector(expected, expname);

        printf("=== alpha - %s ===\n", obsname);
        print_vector(obs, obsname);

        printf("max abserr = %g  max relerr = %g\n", max_abserr, max_relerr);

        printf("=========================================\n\n");
      }

    return k;
} /* test_alpha() */

int
test_beta(gsl_vector *obs, gsl_vector *expected,
          gsl_matrix *A, gsl_matrix *B, const char *obsname,
          const char *expname)
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
      double z = gsl_vector_get(expected, i);
      max = GSL_MAX_DBL(max, fabs(z));
    }

  for (i = 0; i < N; ++i)
    {
      double v_obs = gsl_vector_get(obs, i);
      double v_exp = gsl_vector_get(expected, i);
      double abserr = fabs(v_obs - v_exp);
      double noise = max * GSL_DBL_EPSILON * N * N;

      max_abserr = GSL_MAX_DBL(max_abserr, abserr);

      if (abserr < noise)
        continue;

      if (abserr > 1.0e-6)
        ++k;
    }

    if (k)
      {
        printf("==== CASE %lu ===========================\n\n", count);

        print_matrix(A, "A");
        print_matrix(B, "B");

        printf("=== beta - %s ===\n", expname);
        printf("%s = [\n", expname);
        gsl_vector_fprintf(stdout, expected, "%.12e");
        printf("]\n");

        printf("=== beta - %s ===\n", obsname);
        printf("%s = [\n", obsname);
        gsl_vector_fprintf(stdout, obs, "%.12e");
        printf("]\n");

        printf("max abserr = %g  max relerr = %g\n", max_abserr, max_relerr);

        printf("=========================================\n\n");
      }

    return k;
} /* test_beta() */

/* test if A = Q S Z^t */
void
test_schur(gsl_matrix *A, gsl_matrix *S, gsl_matrix *Q, gsl_matrix *Z)
{
  const size_t N = A->size1;
  gsl_matrix *T1, *T2;
  size_t i, j, k;
  double lhs, rhs;
  double abserr;

  T1 = gsl_matrix_alloc(N, N);
  T2 = gsl_matrix_alloc(N, N);

  /* compute T1 = S Z^t */
  gsl_blas_dgemm(CblasNoTrans,
                 CblasTrans,
                 1.0,
                 S,
                 Z,
                 0.0,
                 T1);

  /* compute T2 = Q T1 = Q S Z^t */
  gsl_blas_dgemm(CblasNoTrans,
                 CblasNoTrans,
                 1.0,
                 Q,
                 T1,
                 0.0,
                 T2);

  k = 0;
  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          lhs = gsl_matrix_get(A, i, j);
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
      print_matrix(S, "S");

      printf("=== Left Schur matrix ===\n");
      print_matrix(Q, "Q");

      printf("=== Right Schur matrix ===\n");
      print_matrix(Z, "Z");

      printf("=== Q S Z^t ===\n");
      print_matrix(T1, "Q S Z^t");

      printf("=== A - Q S Z^t ===\n");
      gsl_matrix_sub(T2, A);
      print_matrix(T1, "A - Q S Z^t");

      printf("=========================================\n\n");
    }

  gsl_matrix_free(T1);
  gsl_matrix_free(T2);
} /* test_schur() */

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

  printf("]\n");
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

  if (!gsl_finite(x[0]))
    return 1;
  if (!gsl_finite(y[0]))
    return -1;

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
main(int argc, char *argv[])
{
  gen_workspace *gen_workspace_p;
  lapack_workspace *lapack_workspace_p;
  size_t N;
  int c;
  int lower;
  int upper;
  int incremental;
  size_t nmat;
  gsl_matrix *A, *B;
  gsl_rng *r;
  int s;
  int compute_schur;
  size_t i;

  gsl_ieee_env_setup();
  gsl_rng_env_setup();

  N = 30;
  lower = -10;
  upper = 10;
  incremental = 0;
  nmat = 0;
  compute_schur = 0;

  while ((c = getopt(argc, argv, "ic:n:l:u:z")) != (-1))
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

          case 'z':
            compute_schur = 1;
            break;

          case '?':
          default:
            printf("usage: %s [-i] [-z] [-n size] [-l lower-bound] [-u upper-bound] [-c num]\n", argv[0]);
            exit(1);
            break;
        } /* switch (c) */
    }

  A = gsl_matrix_alloc(N, N);
  B = gsl_matrix_alloc(N, N);
  gen_workspace_p = gen_alloc(N, compute_schur);
  lapack_workspace_p = lapack_alloc(N);

  r = gsl_rng_alloc(gsl_rng_default);

  if (incremental)
    {
      make_start_matrix(A, lower);

      /* we need B to be non-singular */
      make_random_integer_matrix(B, r, lower, upper);
    }

  fprintf(stderr, "testing N = %d", N);
  if (incremental)
    fprintf(stderr, " incrementally");
  else
    fprintf(stderr, " randomly");

  fprintf(stderr, " on element range [%d, %d]", lower, upper);

  if (compute_schur)
    fprintf(stderr, ", with Schur vectors");

  fprintf(stderr, "\n");

  while (1)
    {
      if (nmat && (count >= nmat))
        break;

      ++count;

      if (!incremental)
        {
          make_random_matrix(A, r, lower, upper);
          make_random_matrix(B, r, lower, upper);
        }
      else
        {
          s = inc_matrix(A, lower, upper);
          if (s)
            break; /* all done */

          make_random_integer_matrix(B, r, lower, upper);
        }

      /*if (count != 89120)
        continue;*/

      /* make copies of matrices */
      gsl_matrix_memcpy(gen_workspace_p->A, A);
      gsl_matrix_memcpy(gen_workspace_p->B, B);
      gsl_matrix_transpose_memcpy(lapack_workspace_p->A, A);
      gsl_matrix_transpose_memcpy(lapack_workspace_p->B, B);

      /* compute eigenvalues with LAPACK */
      s = lapack_proc(lapack_workspace_p);

      if (s != GSL_SUCCESS)
        {
          printf("LAPACK failed, case %lu\n", count);
          exit(1);
        }

#if 0
      print_matrix(A, "A");
      print_matrix(B, "B");
      gsl_matrix_transpose(lapack_workspace_p->A);
      gsl_matrix_transpose(lapack_workspace_p->B);
      print_matrix(lapack_workspace_p->A, "S_lapack");
      print_matrix(lapack_workspace_p->B, "T_lapack");
#endif

      /* compute eigenvalues with GSL */
      s = gen_proc(gen_workspace_p);

      if (s != GSL_SUCCESS)
        {
          printf("=========== CASE %lu ============\n", count);
          printf("Failed to converge: found %u eigenvalues\n",
                 gen_workspace_p->n_evals);
          print_matrix(A, "A");
          print_matrix(B, "B");
          print_matrix(gen_workspace_p->A, "Af");
          print_matrix(gen_workspace_p->B, "Bf");
          print_matrix(lapack_workspace_p->A, "Ae");
          print_matrix(lapack_workspace_p->B, "Be");
          exit(1);
        }

#if 0
      print_matrix(gen_workspace_p->A, "S_gsl");
      print_matrix(gen_workspace_p->B, "T_gsl");
#endif

      /* compute alpha / beta vectors */
      for (i = 0; i < N; ++i)
        {
          double beta;
          gsl_complex alpha, z;

          beta = gsl_vector_get(gen_workspace_p->beta, i);
          if (beta == 0.0)
            GSL_SET_COMPLEX(&z, GSL_POSINF, GSL_POSINF);
          else
            {
              alpha = gsl_vector_complex_get(gen_workspace_p->alpha, i);
              z = gsl_complex_div_real(alpha, beta);
            }

          gsl_vector_complex_set(gen_workspace_p->evals, i, z);

          beta = gsl_vector_get(lapack_workspace_p->beta, i);
          GSL_SET_COMPLEX(&alpha,
                          lapack_workspace_p->alphar[i],
                          lapack_workspace_p->alphai[i]);

          if (beta == 0.0)
            GSL_SET_COMPLEX(&z, GSL_POSINF, GSL_POSINF);
          else
            z = gsl_complex_div_real(alpha, beta);

          gsl_vector_complex_set(lapack_workspace_p->evals, i, z);
          gsl_vector_complex_set(lapack_workspace_p->alpha, i, alpha);
        }

#if 0
      gsl_sort_vector(gen_workspace_p->beta);
      gsl_sort_vector(lapack_workspace_p->beta);
      sort_complex_vector(gen_workspace_p->alpha);
      sort_complex_vector(lapack_workspace_p->alpha);

      s = test_alpha(gen_workspace_p->alpha,
                     lapack_workspace_p->alpha,
                     A,
                     B,
                     "gen",
                     "lapack");
      s = test_beta(gen_workspace_p->beta,
                    lapack_workspace_p->beta,
                    A,
                    B,
                    "gen",
                    "lapack");
#endif
#if 1
      sort_complex_vector(gen_workspace_p->evals);
      sort_complex_vector(lapack_workspace_p->evals);

      s = test_evals(gen_workspace_p->evals,
                     lapack_workspace_p->evals,
                     A,
                     B,
                     "gen",
                     "lapack");
#endif

      if (compute_schur)
        {
          test_schur(A,
                     gen_workspace_p->A,
                     gen_workspace_p->Q,
                     gen_workspace_p->Z);
          test_schur(B,
                     gen_workspace_p->B,
                     gen_workspace_p->Q,
                     gen_workspace_p->Z);
        }
    }

  gsl_matrix_free(A);
  gsl_matrix_free(B);
  gen_free(gen_workspace_p);
  lapack_free(lapack_workspace_p);

  if (r)
    gsl_rng_free(r);

  return 0;
} /* main() */
