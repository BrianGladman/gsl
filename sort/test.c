#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

int cmp_dbl (const void *a, const void *b);
void test_sort (size_t N);
void initialize (double *data, size_t N);
void cpy (double *dest, double *src, size_t N);
void randomize (double *data, size_t n);
void reverse (double *data, size_t N);
int check (double *data, double *orig, size_t N);
int pcheck (gsl_permutation * p, double *data, double *orig, size_t N);
double urand (void);

#define BASE_LONG_DOUBLE
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_CHAR

int
main (void)
{
  size_t i, s;

  /* Test for lengths of 1 ... 31, then 32, 64, 128, 256, ... */

  for (i = 1; i < 1024; i = (i < 32) ? i + 1 : 2 * i)
    test_sort (i);

  for (i = 1; i < 1024; i = (i < 32) ? i + 1 : 2 * i)
    {
      for (s = 1; s < 4; s++)
	{
	  test_sort_vector (i, s);
	  test_sort_vector_float (i, s);
	  test_sort_vector_long_double (i, s);
	  test_sort_vector_ulong (i, s);
	  test_sort_vector_long (i, s);
	  test_sort_vector_uint (i, s);
	  test_sort_vector_int (i, s);
	  test_sort_vector_ushort (i, s);
	  test_sort_vector_short (i, s);
	  test_sort_vector_uchar (i, s);
	  test_sort_vector_char (i, s);
	}
    }

  return gsl_test_summary ();
}

void
test_sort (size_t N)
{
  int status;

  double *orig = (double *) malloc (N * sizeof (double));
  double *data = (double *) malloc (N * sizeof (double));
  gsl_permutation *p = gsl_permutation_alloc (N);

  initialize (orig, N);

  /* Already sorted */
  cpy (data, orig, N);

  status = gsl_sort_index (p, data, N, sizeof (double), (gsl_comparison_fn_t) & cmp_dbl);
  status |= pcheck (p, data, orig, N);
  gsl_test (status, "indexing array, n = %u, ordered", N);

  gsl_sort (data, N, sizeof (double), (gsl_comparison_fn_t) & cmp_dbl);
  status = check (data, orig, N);

  gsl_test (status, "sorting, array, n = %u, ordered", N);

  /* Reverse the data */

  cpy (data, orig, N);
  reverse (data, N);

  status = gsl_sort_index (p, data, N, sizeof (double), (gsl_comparison_fn_t) & cmp_dbl);
  status |= pcheck (p, data, orig, N);
  gsl_test (status, "indexing array, n = %u, reversed", N);

  gsl_sort (data, N, sizeof (double), (gsl_comparison_fn_t) & cmp_dbl);
  status = check (data, orig, N);

  gsl_test (status, "sorting, array, n = %u, reversed", N);

  /* Perform some shuffling */

  cpy (data, orig, N);
  randomize (data, N);

  status = gsl_sort_index (p, data, N, sizeof (double), (gsl_comparison_fn_t) & cmp_dbl);
  status |= pcheck (p, data, orig, N);
  gsl_test (status, "indexing array, n = %u, randomized", N);

  gsl_sort (data, N, sizeof (double), (gsl_comparison_fn_t) & cmp_dbl);
  status = check (data, orig, N);

  gsl_test (status, "sorting, array, n = %u, randomized", N);

  free (orig);
  free (data);
  gsl_permutation_free (p);
}

void
initialize (double *data, size_t N)
{
  size_t i;

  for (i = 0; i < N; i++)
    {
      data[i] = i;
    }
}

void
cpy (double *dest, double *src, size_t N)
{
  size_t i;

  for (i = 0; i < N; i++)
    {
      dest[i] = src[i];
    }
}

void
randomize (double *data, size_t N)
{
  size_t i;

  for (i = 0; i < N; i++)
    {
      size_t j = urand () * N;
      double tmp = data[i];
      data[i] = data[j];
      data[j] = tmp;
    }
}

void
reverse (double *data, size_t N)
{
  size_t i;

  for (i = 0; i < N / 2; i++)
    {
      size_t j = N - i - 1;

      {
	double tmp = data[i];
	data[i] = data[j];
	data[j] = tmp;
      }
    }
}

int
check (double *data, double *orig, size_t N)
{
  size_t i;

  for (i = 0; i < N; i++)
    {
      if (data[i] != orig[i])
	{
	  return GSL_FAILURE;
	}
    }

  return GSL_SUCCESS;
}

int
pcheck (gsl_permutation * p, double *data, double *orig, size_t N)
{
  size_t i;

  for (i = 0; i < N; i++)
    {
      if (data[p->data[i]] != orig[i])
	{
	  return GSL_FAILURE;
	}
    }

  return GSL_SUCCESS;
}


double 
urand (void)
{
  static unsigned long int x = 1;
  x = (1103515245 * x + 12345) & 0x7fffffffUL;
  return x / 2147483648.0;
}

int
cmp_dbl (const void *a, const void *b)
{
  const double x = *(const double *) a;
  const double y = *(const double *) b;
  if (x > y)
    return 1;
  if (x == y)
    return 0;
  else
    return -1;
}
