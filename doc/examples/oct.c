#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "oct.h"

void
print_octave(const gsl_matrix *m, const char *str)
{
  FILE *fp;
  size_t i, j;
  const size_t N = m->size1;
  const size_t M = m->size2;

  if (str == NULL)
    fp = stdout;
  else
    fp = fopen(str, "w");

  if (!fp)
    return;

  fprintf(fp, "%% Created by Octave 2.1.73, Tue Aug 01 15:00:27 2006 MDT <blah@blah>\n");
  fprintf(fp, "%% name: %s\n", str);
  fprintf(fp, "%% type: matrix\n");
  fprintf(fp, "%% rows: %zu\n", N);
  fprintf(fp, "%% columns: %zu\n", M);

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < M; ++j)
        {
          fprintf(fp,
                  "%10.12e%s",
                  gsl_matrix_get(m, i, j),
                  (j < M - 1) ? " " : "\n");
        }
    }

  if (str != NULL)
    fclose(fp);
}

void
printv_octave(const gsl_vector *v, const char *str)
{
  FILE *fp;
  size_t i;
  const size_t N = v->size;

  fp = fopen(str, "w");

  if (!fp)
    return;

  fprintf(fp, "%% Created by Octave 2.1.73, Tue Aug 01 15:00:27 2006 MDT <blah@blah>\n");
  fprintf(fp, "%% name: %s\n", str);
  fprintf(fp, "%% type: matrix\n");
  fprintf(fp, "%% rows: %zu\n", N);
  fprintf(fp, "%% columns: %u\n", 1);

  for (i = 0; i < N; ++i)
    {
      fprintf(fp, "%10.12e\n", gsl_vector_get(v, i));
    }

  fclose(fp);
}

void
printc_octave(const gsl_matrix_complex *m, const char *str)
{
  FILE *fp;
  size_t i, j;
  const size_t N = m->size1;
  const size_t M = m->size2;

  fp = fopen(str, "w");

  if (!fp)
    return;

  fprintf(fp, "%% Created by Octave 2.1.73, Tue Aug 01 15:00:27 2006 MDT <blah@blah>\n");
  fprintf(fp, "%% name: %s\n", str);
  fprintf(fp, "%% type: complex matrix\n");
  fprintf(fp, "%% rows: %zu\n", N);
  fprintf(fp, "%% columns: %zu\n", M);

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < M; ++j)
        {
          gsl_complex z = gsl_matrix_complex_get(m, i, j);
          fprintf(fp,
                  "(%.12e,%.12e)%s",
                  GSL_REAL(z),
                  GSL_IMAG(z),
                  (j < M - 1) ? " " : "\n");
        }
    }

  fclose(fp);
}

/* print symmetric matrix, using upper triangle */
void
printsym_octave(const gsl_matrix *m, const char *str)
{
  FILE *fp;
  size_t i, j;
  const size_t N = m->size1;
  const size_t M = m->size2;

  fp = fopen(str, "w");

  if (!fp)
    return;

  fprintf(fp, "%% Created by Octave 2.1.73, Tue Aug 01 15:00:27 2006 MDT <blah@blah>\n");
  fprintf(fp, "%% name: %s\n", str);
  fprintf(fp, "%% type: matrix\n");
  fprintf(fp, "%% rows: %zu\n", N);
  fprintf(fp, "%% columns: %zu\n", M);

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < M; ++j)
        {
          double z;

          if (j >= i)
            z = gsl_matrix_get(m, i, j);
          else
            z = gsl_matrix_get(m, j, i);

          fprintf(fp,
                  "%.12e%s",
                  z,
                  (j < M - 1) ? " " : "\n");
        }
    }

  fclose(fp);
}

/* print Hermitian matrix, using upper triangle */
void
printherm_octave(const gsl_matrix_complex *m, const char *str)
{
  FILE *fp;
  size_t i, j;
  const size_t N = m->size1;
  const size_t M = m->size2;

  fp = fopen(str, "w");

  if (!fp)
    return;

  fprintf(fp, "%% Created by Octave 2.1.73, Tue Aug 01 15:00:27 2006 MDT <blah@blah>\n");
  fprintf(fp, "%% name: %s\n", str);
  fprintf(fp, "%% type: complex matrix\n");
  fprintf(fp, "%% rows: %zu\n", N);
  fprintf(fp, "%% columns: %zu\n", M);

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < M; ++j)
        {
          gsl_complex z;

          if (j >= i)
            z = gsl_matrix_complex_get(m, i, j);
          else
            z = gsl_complex_conjugate(gsl_matrix_complex_get(m, j, i));

          fprintf(fp,
                  "(%.12e,%.12e)%s",
                  GSL_REAL(z),
                  GSL_IMAG(z),
                  (j < M - 1) ? " " : "\n");
        }
    }

  fclose(fp);
}

void
printcv_octave(const gsl_vector_complex *v, const char *str)
{
  FILE *fp;
  size_t i;
  const size_t N = v->size;

  fp = fopen(str, "w");

  if (!fp)
    return;

  fprintf(fp, "%% Created by Octave 2.1.73, Tue Aug 01 15:00:27 2006 MDT <blah@blah>\n");
  fprintf(fp, "%% name: %s\n", str);
  fprintf(fp, "%% type: complex matrix\n");
  fprintf(fp, "%% rows: %zu\n", N);
  fprintf(fp, "%% columns: %u\n", 1);

  for (i = 0; i < N; ++i)
    {
      gsl_complex z = gsl_vector_complex_get(v, i);
      fprintf(fp, "(%.12e,%.12e)\n", GSL_REAL(z), GSL_IMAG(z));
    }

  fclose(fp);
}

/* print triangular matrix, using upper triangle */
void
printtri_octave(const gsl_matrix *m, const char *str)
{
  FILE *fp;
  size_t i, j;
  const size_t N = m->size1;
  const size_t M = m->size2;

  if (str == NULL)
    fp = stdout;
  else
    fp = fopen(str, "w");

  if (!fp)
    return;

  fprintf(fp, "%% Created by Octave 2.1.73, Tue Aug 01 15:00:27 2006 MDT <blah@blah>\n");
  fprintf(fp, "%% name: %s\n", str);
  fprintf(fp, "%% type: matrix\n");
  fprintf(fp, "%% rows: %zu\n", N);
  fprintf(fp, "%% columns: %zu\n", M);

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < M; ++j)
        {
          double mij;

          if (j >= i)
            mij = gsl_matrix_get(m, i, j);
          else
            mij = 0.0;

          fprintf(fp,
                  "%10.12e%s",
                  mij,
                  (j < M - 1) ? " " : "\n");
        }
    }

  if (str != NULL)
    fclose(fp);
}
