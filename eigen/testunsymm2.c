/*
 * testunsymm.c
 * Patrick Alken
 *
 * Compile: gcc -g -O2 -Wall -lm -lgsl -llapack -lf77blas -lcblas -latlas -lg2c
 *
 * Usage: testunsymm [options]
 *
 * -i             : incremental matrices
 * -z             : compute Schur vectors and test them
 * -n size        : size of matrices
 * -l lower-bound : lower bound for matrix elements
 * -u upper-bound : upper bound for matrix elements
 * -t count       : print time comparisons every count matrices
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
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#undef TTQRE
#undef LAPACK_HQR

#define FOR_TRUE       1
#define FOR_FALSE      0

typedef struct
{
  gsl_eigen_unsymm_workspace *unsymm_p;
  gsl_matrix *A;
  gsl_vector_complex *eval;

  int compute_z;
  gsl_matrix *Z;
  gsl_matrix *T1;
  gsl_matrix *T2;
} unsymm_workspace;

unsymm_workspace *unsymm_alloc(size_t n, int compute_z);
void unsymm_free(unsymm_workspace *w);
int unsymm_proc(unsymm_workspace *w);

typedef struct
{
  gsl_matrix *A;
  gsl_vector *wr;
  gsl_vector *wi;
  gsl_vector_complex *eval;
  gsl_vector *work;

  gsl_matrix *Z;

  int info;
  int lwork;
  int sdim;
  int ldvs;
  int ilo;
  int ihi;

  char jobvs;
  char sort;
  char jobg;

  /* dlahqr */
  char jobhqr;
  int wantt;
  int wantz;
  int iloz;
  int ihiz;
  int ldz;
  gsl_vector *work2;
} lapack_workspace;

lapack_workspace *lapack_alloc(size_t n, int compute_z);
void lapack_free(lapack_workspace *w);
int lapack_proc(lapack_workspace *w);

void dgees_(char *jobvs, char *sort, int *sel, int *n,
            double *a, int *lda, int *sdim, double *wr, double *wi,
            double *vs, int *ldvs, double *work, int *lwork,
            int *bwork, int *info);
void dgebal_(char *job, int *n, double *a, int *lda, int *ilo,
             int *ihi, double *scale, int *info);
void dlahqr_(int *wantt, int *wantz, int *n, int *ilo,
             int *ihi, double *h, int *ldh, double *wr,
             double *wi, int *iloz, int *ihiz, double *z,
             int *ldz, int *info);

#ifdef TTQRE

typedef struct {
  gsl_matrix *A;
  gsl_vector *wr;
  gsl_vector *wi;
  gsl_vector_complex *eval;
  gsl_vector *work;
  gsl_vector *work2;

  int info;
  int lwork;
  int sdim;
  int ldvs;
  int ilo;
  int ihi;
  int ldz;

  char job;
  char compz;
  char sort;
  char jobg;
} ttqre_workspace;

ttqre_workspace *ttqre_alloc(size_t n);
void ttqre_free(ttqre_workspace *w);
int ttqre_proc(ttqre_workspace *w);

void ttqre_(char *job, char *compz, int *n, int *ilo, int *ihi,
            double *h, int *ldh, double *wr, double *wi,
            double *z, int *ldz, double *work, int *lwork,
            int *info);

#endif /* TTQRE */

#if defined(TTQRE) || defined(LAPACK_HQR)

void dgehrd_(int *n, int *ilo, int *ihi, double *a, int *lda,
             double *tau, double *work, int *lwork, int *info);

#endif

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
void print_evals(gsl_vector_complex *eval);
int cmp(double a, double b);
int compare(const void *a, const void *b);
void sort_complex_vector(gsl_vector_complex *v);
int test_evals(gsl_vector_complex *obs, gsl_vector_complex *expected,
               gsl_matrix *A, const char *obsname);
int test_Z(gsl_matrix *A, gsl_matrix *Z, gsl_matrix *T, gsl_matrix *T1,
           gsl_matrix *T2);
void my_error_handler(const char *reason, const char *file, int line,
                      int err);

/**********************************************
 * GSL routines
 **********************************************/

unsymm_workspace *
unsymm_alloc(size_t n, int compute_z)

{
  unsymm_workspace *w;

  w = (unsymm_workspace *) malloc(sizeof(unsymm_workspace));

  memset(w, '\0', sizeof(unsymm_workspace));

  w->unsymm_p = gsl_eigen_unsymm_alloc(n);

  w->A = gsl_matrix_alloc(n, n);

  if (compute_z)
    {
      w->Z = gsl_matrix_alloc(n, n);
      w->T1 = gsl_matrix_alloc(n, n);
      w->T2 = gsl_matrix_alloc(n, n);
      w->compute_z = 1;
    }

  w->eval = gsl_vector_complex_alloc(n);

  return (w);
} /* unsymm_alloc() */

void
unsymm_free(unsymm_workspace *w)

{
  if (!w)
    return;

  if (w->unsymm_p)
    gsl_eigen_unsymm_free(w->unsymm_p);

  if (w->A)
    gsl_matrix_free(w->A);

  if (w->Z)
    gsl_matrix_free(w->Z);

  if (w->T1)
    gsl_matrix_free(w->T1);

  if (w->T2)
    gsl_matrix_free(w->T2);

  if (w->eval)
    gsl_vector_complex_free(w->eval);

  free(w);
} /* unsymm_free() */

int
unsymm_proc(unsymm_workspace *w)

{
  int s;

  if (w->compute_z)
  {
    gsl_eigen_unsymm_params(1, 0, w->unsymm_p);
    s = gsl_eigen_unsymm_Z(w->A, w->eval, w->Z, w->unsymm_p);
  }
  else
    s = gsl_eigen_unsymm(w->A, w->eval, w->unsymm_p);

  return s;
} /* unsymm_proc() */

/**********************************************
 * LAPACK routines
 **********************************************/

lapack_workspace *
lapack_alloc(size_t n, int compute_z)

{
  lapack_workspace *w;
  size_t size;

  w = (lapack_workspace *) malloc(sizeof(lapack_workspace));

  memset(w, '\0', sizeof(lapack_workspace));

  w->A = gsl_matrix_alloc(n, n);
  w->wr = gsl_vector_alloc(n);
  w->wi = gsl_vector_alloc(n);
  w->eval = gsl_vector_complex_alloc(n);
  w->work = gsl_vector_calloc(3 * n);
  w->work2 = gsl_vector_calloc(3 * n);
  w->Z = gsl_matrix_alloc(n, n);

  w->info = 0;
  w->lwork = -1;
  w->sdim = 0;
  w->ldvs = 1;
  w->sort = 'N';
  w->jobg = 'B';

  if (compute_z)
    {
      w->jobvs = 'V';
      w->ldvs = w->Z->tda;
    }
  else
    {
      w->jobvs = 'N';
    }

  w->wantt = FOR_FALSE;
  w->wantz = compute_z ? FOR_TRUE : FOR_FALSE;

  dgebal_((char *) &w->jobg,
          (int *) &n,
          &(w->A->data[0]),
          (int *) &w->A->tda,
          &w->ilo,
          &w->ihi,
          &(w->wr->data[0]),
          (int *) &w->info);

  dgees_((char *) &w->jobvs,
         (char *) &w->sort,
         (int *) NULL,
         (int *) &n,
         w->A->data,
         (int *) &w->A->tda,
         (int *) &w->sdim,
         w->wr->data,
         w->wi->data,
         w->Z->data,
         (int *) &w->Z->tda,
         w->work->data,
         (int *) &w->lwork,
         (int *) NULL,
         (int *) &w->info);

  size = w->work->data[0];

  if (size > w->work->size)
    {
      gsl_vector_free(w->work);
      gsl_vector_free(w->work2);
      w->work = gsl_vector_calloc(size);
      w->work2 = gsl_vector_calloc(size);
    }

  return (w);
}

void
lapack_free(lapack_workspace *w)

{
  gsl_matrix_free(w->A);
  gsl_vector_free(w->wr);
  gsl_vector_free(w->wi);
  gsl_vector_complex_free(w->eval);
  gsl_vector_free(w->work);
  gsl_vector_free(w->work2);
  gsl_matrix_free(w->Z);

  free(w);
}

int
lapack_proc(lapack_workspace *w)

{
  size_t N = w->A->size1;
  size_t i;
  gsl_complex z;

  w->info = 0;
  w->lwork = w->work->size;
  w->sdim = 0;
  w->ldvs = 1;

  w->wantt = FOR_FALSE;

#ifdef LAPACK_HQR

  /* balance */
  dgebal_((char *) &w->jobg,
          (int *) &N,
          w->A->data,
          (int *) &w->A->tda,
          &w->ilo,
          &w->ihi,
          w->wr->data,
          (int *) &w->info);

  /* convert to hessenberg */
  dgehrd_((int *) &N,
          &w->ilo,
          &w->ihi,
          w->A->data,
          (int *) &w->A->tda,
          w->work2->data,
          w->work->data,
          (int *) &w->lwork,
          (int *) &w->info);

  w->iloz = w->ilo;
  w->ihiz = w->ihi;

  dlahqr_(&w->wantt,
          &w->wantz,
          (int *) &N,
          &w->ilo,
          &w->ihi,
          w->A->data,
          (int *) &w->A->tda,
          w->wr->data,
          w->wi->data,
          &w->iloz,
          &w->ihiz,
          w->Z->data,
          (int *) &w->ldvs,
          (int *) &w->info);

#else

  dgees_((char *) &w->jobvs,
         (char *) &w->sort,
         (int *) NULL,
         (int *) &N,
         w->A->data,
         (int *) &w->A->tda,
         (int *) &w->sdim,
         w->wr->data,
         w->wi->data,
         w->Z->data,
         (int *) &(w->Z->tda),
         w->work->data,
         (int *) &w->lwork,
         (int *) NULL,
         (int *) &w->info);

#endif

  if (w->info != 0)
    printf("lapack_proc: info = %d\n", w->info);

  for (i = 0; i < N; ++i)
    {
      GSL_SET_COMPLEX(&z, w->wr->data[i], w->wi->data[i]);
      gsl_vector_complex_set(w->eval, i, z);
    }

  return (w->info);
} /* lapack_proc() */

#ifdef TTQRE

/**********************************************
 * TTQRE routines
 **********************************************/

ttqre_workspace *
ttqre_alloc(size_t n)

{
  ttqre_workspace *w;
  size_t size;

  w = (ttqre_workspace *) malloc(sizeof(ttqre_workspace));

  memset(w, '\0', sizeof(ttqre_workspace));

  w->A = gsl_matrix_alloc(n, n);
  w->wr = gsl_vector_alloc(n);
  w->wi = gsl_vector_alloc(n);
  w->eval = gsl_vector_complex_alloc(n);
  w->work = gsl_vector_calloc(3 * n);
  w->work2 = gsl_vector_calloc(3 * n);

  w->info = 0;
  w->lwork = -1;
  w->job = 'E';
  w->compz = 'N';
  w->jobg = 'B';
  w->ldz = 1;

  dgebal_((char *) &w->jobg,
          (int *) &n,
          w->A->data,
          (int *) &w->A->tda,
          &w->ilo,
          &w->ihi,
          w->wr->data,
          (int *) &w->info);

  ttqre_((char *) &w->job,
         (char *) &w->compz,
         (int *) &n,
         &w->ilo,
         &w->ihi,
         w->A->data,
         (int *) &w->A->tda,
         w->wr->data,
         w->wi->data,
         (double *) NULL,
         (int *) &w->ldz,
         w->work->data,
         (int *) &w->lwork,
         (int *) &w->info);

  size = w->work->data[0];

  if (size > w->work->size)
    {
      gsl_vector_free(w->work);
      gsl_vector_free(w->work2);
      w->work = gsl_vector_calloc(size);
      w->work2 = gsl_vector_calloc(size);
    }

  return (w);
} /* ttqre_alloc() */

void
ttqre_free(ttqre_workspace *w)

{
  gsl_matrix_free(w->A);
  gsl_vector_free(w->wr);
  gsl_vector_free(w->wi);
  gsl_vector_complex_free(w->eval);
  gsl_vector_free(w->work);
  gsl_vector_free(w->work2);

  free(w);
} /* ttqre_free() */

int
ttqre_proc(ttqre_workspace *w)

{
  size_t N = w->A->size1;
  size_t i;
  gsl_complex z;

  w->info = 0;
  w->lwork = w->work->size;

  /* balance */
  dgebal_((char *) &w->jobg,
          (int *) &N,
          w->A->data,
          (int *) &w->A->tda,
          &w->ilo,
          &w->ihi,
          w->wr->data,
          (int *) &w->info);

  /* convert to hessenberg */
  dgehrd_((int *) &N,
          &w->ilo,
          &w->ihi,
          w->A->data,
          (int *) &w->A->tda,
          w->work2->data,
          w->work->data,
          (int *) &w->lwork,
          (int *) &w->info);

  ttqre_((char *) &w->job,
         (char *) &w->compz,
         (int *) &N,
         &w->ilo,
         &w->ihi,
         w->A->data,
         (int *) &w->A->tda,
         w->wr->data,
         w->wi->data,
         (double *) NULL,
         (int *) &w->ldz,
         w->work->data,
         (int *) &w->lwork,
         (int *) &w->info);

  if (w->info != 0)
    printf("ttqre_proc: info = %d\n", w->info);

  for (i = 0; i < N; ++i)
    {
      GSL_SET_COMPLEX(&z, w->wr->data[i], w->wi->data[i]);
      gsl_vector_complex_set(w->eval, i, z);
    }

  return (w->info);
} /* ttqre_proc() */

#endif /* TTQRE */

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
        printf("%10.9f%s",
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
print_evals(gsl_vector_complex *eval)

{
  size_t N = eval->size;
  size_t i;
  gsl_complex z;

  printf("EV = [\n");

  for (i = 0; i < N; ++i)
    {
      z = gsl_vector_complex_get(eval, i);
      printf("%.18e %.18e;\n", GSL_REAL(z), GSL_IMAG(z));
    }

  printf("\n]\n");
} /* print_evals() */

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
           gsl_matrix *A, const char *obsname)

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

        printf("=== eval - lapack dgees ===\n");
        print_evals(expected);

        printf("=== eval - %s ===\n", obsname);
        print_evals(obs);

        printf("max abserr = %g  max relerr = %g\n", max_abserr, max_relerr);

        printf("=========================================\n\n");
      }

    return k;
} /* test_evals() */

/* test if A = ZTZ^t (or AZ = ZT) */

int
test_Z(gsl_matrix *A, gsl_matrix *Z, gsl_matrix *T, gsl_matrix *T1,
       gsl_matrix *T2)

{
  size_t N = A->size1;
  size_t i, j, k;
  double rhs, lhs;
  double abserr;

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

  return k;
} /* test_Z() */

void
my_error_handler(const char *reason, const char *file, int line,
                 int err)

{
  printf("[caught: %s:%d: errno=%d %s]\n", file, line, err, reason);
} /* my_error_handler() */

int
main(int argc, char *argv[])

{
  unsymm_workspace *unsymm_workspace_p;
  lapack_workspace *lapack_workspace_p;
  size_t N;
  int c;
  gsl_matrix *A;
  gsl_rng *r;
  int incremental; /* incremental/random matrices */
  int lower;  /* lower bound */
  int upper;  /* upper bound */
  unsigned int nmat;   /* number of matrices to solve */
  int tint;   /* interval between time comparisons */
  int status;
  int compute_z;
  struct tms lapack_start, lapack_end;
  struct tms unsymm_start, unsymm_end;
  double unsymm_sum, lapack_sum;
  double clock_ticks_per_second;
  gsl_complex z_zero;
#ifdef TTQRE
  ttqre_workspace *ttqre_workspace_p;
  struct tms ttqre_start, ttqre_end;
  double ttqre_sum = 0.0;
#endif

  gsl_ieee_env_setup();
  gsl_rng_env_setup();

  /*gsl_set_error_handler(&my_error_handler);*/

  clock_ticks_per_second = (double) sysconf(_SC_CLK_TCK);
  lapack_sum = unsymm_sum = 0.0;

  N = 30;
  incremental = 0;
  compute_z = 0;
  lower = -10;
  upper = 10;
  nmat = 0;
  tint = 0;

  while ((c = getopt(argc, argv, "izc:n:l:u:t:")) != (-1))
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

          case 't':
            tint = strtol(optarg, NULL, 0);
            break;

          case 'z':
            compute_z = 1;
            break;

          case '?':
          default:
            printf("usage: %s [-i] [-z] [-n size] [-l lower-bound] [-u upper-bound] [-t count] [-c num]\n", argv[0]);
            exit(1);
            break;
        } /* switch (c) */
    }

  A = gsl_matrix_alloc(N, N);
  unsymm_workspace_p = unsymm_alloc(N, compute_z);
  lapack_workspace_p = lapack_alloc(N, compute_z);
#ifdef TTQRE
  ttqre_workspace_p = ttqre_alloc(N);
#endif

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

  fprintf(stderr, " on element range [%d, %d]\n", lower, upper);

  GSL_SET_COMPLEX(&z_zero, 0.0, 0.0);

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

      /*if (count != 1437)
        continue;*/

      /* make copies of the matrix */
      gsl_matrix_memcpy(lapack_workspace_p->A, A);
      gsl_matrix_memcpy(unsymm_workspace_p->A, A);
#ifdef TTQRE
      gsl_matrix_memcpy(ttqre_workspace_p->A, A);
#endif

      /* compute eigenvalues with LAPACK */

      times(&lapack_start);
      status = lapack_proc(lapack_workspace_p);
      times(&lapack_end);

      lapack_sum += (lapack_end.tms_utime - lapack_start.tms_utime) /
                    clock_ticks_per_second;

      if (status)
        continue;

      sort_complex_vector(lapack_workspace_p->eval);

#ifdef TTQRE
      /* compute eigenvalues with TTQRE */

      times(&ttqre_start);
      status = ttqre_proc(ttqre_workspace_p);
      times(&ttqre_end);

      ttqre_sum += (ttqre_end.tms_utime - ttqre_start.tms_utime) /
                   clock_ticks_per_second;

      sort_complex_vector(ttqre_workspace_p->eval);

      test_evals(ttqre_workspace_p->eval,
                 lapack_workspace_p->eval,
                 A,
                 "ttqre");
#endif /* TTQRE */

      /* compute eigenvalues with my GSL code */

      gsl_vector_complex_set_all(unsymm_workspace_p->eval, z_zero);

      times(&unsymm_start);
      status = unsymm_proc(unsymm_workspace_p);
      times(&unsymm_end);

      /*printf("eigenvalues found = %d\n", status);*/
      /*printf("tot its = %u\n", unsymm_workspace_p->unsymm_p->sbaed_workspace_p->state.tot_it);
      printf("qr its = %u\n", unsymm_workspace_p->unsymm_p->sbaed_workspace_p->state.qr_it);*/

      if (status != (int) N)
        {
          printf("=========== CASE %lu ============\n", count);
          printf("Failed to converge: found %d eigenvalues\n", status);
          print_matrix(A, "A");
          print_matrix(unsymm_workspace_p->A, "T");
          if (compute_z)
            print_matrix(unsymm_workspace_p->Z, "Z");
        }

      unsymm_sum += (unsymm_end.tms_utime - unsymm_start.tms_utime) /
                    clock_ticks_per_second;

      sort_complex_vector(unsymm_workspace_p->eval);

      c = test_evals(unsymm_workspace_p->eval,
                     lapack_workspace_p->eval,
                     A,
                     "alken");
      /*if (c)
        exit(1);*/

      if (compute_z)
        {
          test_Z(A,
                 unsymm_workspace_p->Z,
                 unsymm_workspace_p->A,
                 unsymm_workspace_p->T1,
                 unsymm_workspace_p->T2);
        }

      if (tint && (count % tint == 0))
        {
          printf("=========== CASE %lu ============\n", count);
          printf("avg lapack time = %g sec\n", lapack_sum / count);
#ifdef TTQRE
          printf("avg ttqre time = %g sec\n", ttqre_sum / count);
#endif
          printf("avg alken time = %g sec\n", unsymm_sum / count);
#ifdef TTQRE
          printf("ttqre/lapack = %g\n", ttqre_sum / lapack_sum);
          printf("ttqre/alken = %g\n", ttqre_sum / unsymm_sum);
#endif
          printf("lapack/alken = %g\n", lapack_sum / unsymm_sum);
        }
    }

  gsl_matrix_free(A);
  unsymm_free(unsymm_workspace_p);
  lapack_free(lapack_workspace_p);
#ifdef TTQRE
  ttqre_free(ttqre_workspace_p);
#endif

  if (r)
    gsl_rng_free(r);

  return 0;
} /* main() */
