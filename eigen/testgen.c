/*
 * testgen.c
 * Patrick Alken
 *
 * Compile: gcc -g -O2 -Wall -o testgen testgen.c -lm -lgsl -lcblas -latlas
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
#include <gsl/gsl_test.h>
#include <gsl/gsl_errno.h>

typedef struct
{
  gsl_eigen_gen_workspace *gen_p;
  gsl_matrix *A;
  gsl_matrix *B;
  gsl_vector_complex *alpha;
  gsl_vector *beta;
  gsl_vector_complex *evals;

  gsl_eigen_genv_workspace *genv_p;
  gsl_matrix *Av;
  gsl_matrix *Bv;
  gsl_vector_complex *alphav;
  gsl_vector *betav;
  gsl_matrix_complex *evec;

  gsl_matrix *Q;
  gsl_matrix *Z;
  int compute_schur;

  size_t n_evals;
} gen_workspace;

gen_workspace *gen_alloc(size_t n, int compute_schur);
void gen_free(gen_workspace *w);
int gen_proc(gen_workspace *w);

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
void output_matrix(const gsl_matrix *m);
void print_matrix(const gsl_matrix *m, const char *str);
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
void print_vector(const gsl_vector_complex *eval, const char *str);
int cmp(double a, double b);
int compare(const void *a, const void *b);
void sort_complex_vector(gsl_vector_complex *v);
int test_eigenvectors(const gsl_matrix *A, const gsl_matrix *B,
                      const gsl_vector_complex *alpha,
                      const gsl_vector *beta,
                      const gsl_matrix_complex *evec);

gen_workspace *
gen_alloc(size_t n, int compute_schur)
{
  gen_workspace *w;

  w = (gen_workspace *) calloc(1, sizeof(gen_workspace));

  w->gen_p = gsl_eigen_gen_alloc(n);
  w->genv_p = gsl_eigen_genv_alloc(n);

  w->A = gsl_matrix_alloc(n, n);
  w->B = gsl_matrix_alloc(n, n);
  w->alpha = gsl_vector_complex_alloc(n);
  w->beta = gsl_vector_alloc(n);
  w->evals = gsl_vector_complex_alloc(n);
  w->Av = gsl_matrix_alloc(n, n);
  w->Bv = gsl_matrix_alloc(n, n);
  w->alphav = gsl_vector_complex_alloc(n);
  w->betav = gsl_vector_alloc(n);
  w->evec = gsl_matrix_complex_alloc(n, n);
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

  if (w->genv_p)
    gsl_eigen_genv_free(w->genv_p);

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

  if (w->Av)
    gsl_matrix_free(w->Av);

  if (w->Bv)
    gsl_matrix_free(w->Bv);

  if (w->alphav)
    gsl_vector_complex_free(w->alphav);

  if (w->betav)
    gsl_vector_free(w->betav);

  if (w->evec)
    gsl_matrix_complex_free(w->evec);

  if (w->Q)
    gsl_matrix_free(w->Q);

  if (w->Z)
    gsl_matrix_free(w->Z);

  free(w);
}

int
gen_proc(gen_workspace *w)
{
  int s1, s2, s;

  s1 = gsl_eigen_gen_QZ(w->A, w->B, w->alpha, w->beta, w->Q, w->Z, w->gen_p);
  s2 = gsl_eigen_genv(w->Av, w->Bv, w->alphav, w->betav, w->evec, w->genv_p);

  w->n_evals = w->gen_p->n_evals;

  s = 0;
  if (s1)
    s = s1;
  else if (s2)
    s = s2;

  return s;
} /* gen_proc() */

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
output_matrix(const gsl_matrix *m)
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
print_matrix(const gsl_matrix *m, const char *str)

{
  size_t i, j;
  size_t N = m->size1;
  size_t M = m->size2;
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

      sprintf(buf, "%s(%u:%u,%u:%u)",
              str,
              i + 1,
              i + r,
              j + 1,
              j + c);

      printf("%s = [\n", buf);

      {
        gsl_matrix_const_view v = gsl_matrix_const_submatrix(m, i, j, r, c);
        output_matrix(&v.matrix);
      }

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
print_vector(const gsl_vector_complex *eval, const char *str)
{
  size_t N = eval->size;
  size_t i;
  gsl_complex z;

  printf("%s = [\n", str);

  for (i = 0; i < N; ++i)
    {
      z = gsl_vector_complex_get(eval, i);
      printf("%.18e + %.18ei;\n", GSL_REAL(z), GSL_IMAG(z));
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
test_eigenvectors(const gsl_matrix *A, const gsl_matrix *B,
                  const gsl_vector_complex *alpha, const gsl_vector *beta,
                  const gsl_matrix_complex *evec)
{
  const size_t N = A->size1;
  size_t i, j;
  int k, s;
  gsl_matrix_complex *ma, *mb;
  gsl_vector_complex *x, *y;
  gsl_complex z_one, z_zero;

  ma = gsl_matrix_complex_alloc(N, N);
  mb = gsl_matrix_complex_alloc(N, N);
  y = gsl_vector_complex_alloc(N);
  x = gsl_vector_complex_alloc(N);

  /* ma <- A, mb <- B */
  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex z;

          GSL_SET_COMPLEX(&z, gsl_matrix_get(A, i, j), 0.0);
          gsl_matrix_complex_set(ma, i, j, z);

          GSL_SET_COMPLEX(&z, gsl_matrix_get(B, i, j), 0.0);
          gsl_matrix_complex_set(mb, i, j, z);
        }
    }

  GSL_SET_COMPLEX(&z_one, 1.0, 0.0);
  GSL_SET_COMPLEX(&z_zero, 0.0, 0.0);

  s = 0;

  /* check eigenvalues */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_complex_const_view vi = gsl_matrix_complex_const_column(evec, i);
      gsl_complex ai = gsl_vector_complex_get(alpha, i);
      double bi = gsl_vector_get(beta, i);
      double norm = gsl_blas_dznrm2(&vi.vector);

      /* check that eigenvector is normalized */
      gsl_test_rel(norm, 1.0, N * GSL_DBL_EPSILON, "case %u, normalized",
                   count);

      /* compute x = alpha * B * v */
      gsl_blas_zgemv(CblasNoTrans, z_one, mb, &vi.vector, z_zero, x);
      gsl_blas_zscal(ai, x);

      /* compute y = beta * A v */
      gsl_blas_zgemv(CblasNoTrans, z_one, ma, &vi.vector, z_zero, y);
      gsl_blas_zdscal(bi, y);

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

          if (fabs(lhs_r - rhs_r) > 1e10 * GSL_DBL_EPSILON)
            ++k;
          if (fabs(lhs_i - rhs_i) > 1e10 * GSL_DBL_EPSILON)
            ++k;
        }

      if (k)
        {
          s++;

          printf("==== CASE %lu ===========================\n\n", count);

          print_matrix(A, "A");
          print_matrix(B, "B");

          printf("alpha = %.10e + %.10ei\n", GSL_REAL(ai), GSL_IMAG(ai));
          printf("beta = %.10e\n", bi);
          printf("alpha/beta = %.10e + %.10ei\n", GSL_REAL(ai)/bi, GSL_IMAG(ai)/bi);

          print_vector(&vi.vector, "v");

          print_vector(y, "beta*A*v");

          print_vector(x, "alpha*B*v");

          printf("=========================================\n\n");
        }
    }

  gsl_matrix_complex_free(ma);
  gsl_matrix_complex_free(mb);
  gsl_vector_complex_free(y);
  gsl_vector_complex_free(x);

  return s;
} /* test_eigenvectors() */

int
main(int argc, char *argv[])
{
  gen_workspace *gen_workspace_p;
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
  /*gsl_set_error_handler_off();*/

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

      /*if (count != 23441872)
        continue;*/

      /* make copies of matrices */
      gsl_matrix_memcpy(gen_workspace_p->A, A);
      gsl_matrix_memcpy(gen_workspace_p->B, B);
      gsl_matrix_memcpy(gen_workspace_p->Av, A);
      gsl_matrix_memcpy(gen_workspace_p->Bv, B);

      /* compute eigenvalues with GSL */
      s = gen_proc(gen_workspace_p);

      if (s != GSL_SUCCESS)
        {
          printf("=========== CASE %lu ============\n", count);
          printf("Failed to converge: found %u eigenvalues\n",
                 gen_workspace_p->n_evals);
          print_matrix(A, "A");
          print_matrix(B, "B");
          print_matrix(gen_workspace_p->Av, "S");
          print_matrix(gen_workspace_p->Bv, "T");
          continue;
        }

      /* compute alpha / beta vectors */
      for (i = 0; i < N; ++i)
        {
          double beta = gsl_vector_get(gen_workspace_p->beta, i);
          gsl_complex alpha =
            gsl_vector_complex_get(gen_workspace_p->alpha, i);
          gsl_complex z;

          if (!gsl_finite(beta) || !gsl_finite(GSL_REAL(alpha)) ||
              !gsl_finite(GSL_IMAG(alpha)))
            {
              printf("nan/inf in element %u of (alpha,beta): alphar = %g, alphai = %g, beta = %g\n",
                     i, GSL_REAL(alpha), GSL_IMAG(alpha), beta);
            }

          if (beta == 0.0)
            GSL_SET_COMPLEX(&z, GSL_POSINF, GSL_POSINF);
          else
            z = gsl_complex_div_real(alpha, beta);

          gsl_vector_complex_set(gen_workspace_p->evals, i, z);
        }

      test_eigenvectors(A,
                        B,
                        gen_workspace_p->alphav,
                        gen_workspace_p->betav,
                        gen_workspace_p->evec);

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
  gsl_rng_free(r);

  return 0;
} /* main() */
