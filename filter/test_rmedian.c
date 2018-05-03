/* filter/test_rmedian.c
 * 
 * Copyright (C) 2018 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifdef _MSC_VER
#  include <gettimeofday.h>
#else
#  include <sys/time.h>
#endif
#define TIMEDIFF(a, b)    ((double) ((b).tv_sec - (a).tv_sec) + 1.0e-6 * ((b).tv_usec - (a).tv_usec))

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* slow/dumb rmedian which constructs actual window for each sample, sorts
 * it and finds median */
static int
slow_rmedian(const gsl_vector * x, gsl_vector * y, const size_t k)
{
  const size_t n = x->size;
  const size_t window_size = (k % 2 == 0) ? k + 1 : k;
  const int H = (int) (window_size / 2);
  double *window = malloc(window_size * sizeof(double));
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double yi;
      int idx;
      size_t j = 0;

      /* fill window with previous filter output values */
      for (idx = (int) i - H; idx < (int) i; ++idx)
        {
          if (idx < 0)
            window[j++] = 0.0; /* zero pad */
          else
            window[j++] = gsl_vector_get(y, idx);
        }

      /* fill remainder of window with input samples */
      for (idx = (int) i; idx <= (int) i + H; ++idx)
        {
          if (idx < (int) n)
            window[j++] = gsl_vector_get(x, idx);
          else
            window[j++] = 0.0; /* zero pad */
        }

      yi = median_find((long) window_size, window);
      gsl_vector_set(y, i, yi);
    }

  free(window);

  return GSL_SUCCESS;
}

/* test square wave input (root signal) */
static void
test_rmedian_root(const gsl_filter_end_t etype, const size_t n, const size_t k)
{
  const double tol = 1.0e-12;
  gsl_filter_rmedian_workspace *w = gsl_filter_rmedian_alloc(k);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  char buf[2048];
  size_t i;

  /* test a root sequence (square input): x = [zero one zero] */
  gsl_vector_set_all(x, 0.0);

  for (i = n / 3; i <= n / 2; ++i)
    gsl_vector_set(x, i, 1.0);

  /* compute y = rmedian(x) and test y = x */
  gsl_filter_rmedian(etype, x, y, w);

  sprintf(buf, "n=%zu k=%zu RMF square wave root sequence", n, k);
  compare_vectors(tol, y, x, buf);

  /* compute y = rmedian(x) with second algorithm and test y = x */
  gsl_filter_rmedian2(x, y, w);

  sprintf(buf, "n=%zu k=%zu RMF square wave root sequence 2nd algo", n, k);
  compare_vectors(tol, y, x, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_filter_rmedian_free(w);
}

/* test random input and in-place */
static void
test_rmedian_random(const gsl_filter_end_t etype, const size_t n, const size_t k, gsl_rng * r)
{
  const double tol = 1.0e-12;
  gsl_filter_rmedian_workspace *w = gsl_filter_rmedian_alloc(k);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *z = gsl_vector_alloc(n);
  char buf[2048];

  /* test filter with random input against slow algorithm */
  random_vector(x, r);

  /* y = rmedian(x) */
  gsl_filter_rmedian(etype, x, y, w);

  /* y = rmedian(x) with slow algorithm */
  slow_rmedian(x, z, w->K);

  /* test y = z */
  sprintf(buf, "n=%zu k=%zu RMF symmetric random slow test", n, k);
  compare_vectors(tol, y, z, buf);

  /* test in-place filter */

  /* z = rmedian(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_filter_rmedian(etype, z, z, w);

  sprintf(buf, "n=%zu k=%zu RMF symmetric random in-place", n, k);
  compare_vectors(tol, z, y, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_filter_rmedian_free(w);
}

void
test_rmedian_sine(const gsl_filter_end_t etype, const size_t n, const size_t k, gsl_rng * r)
{
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y_rmedian = gsl_vector_alloc(n);
  gsl_vector *z_rmedian = gsl_vector_alloc(n);
  gsl_filter_rmedian_workspace *w = gsl_filter_rmedian_alloc(k);
  size_t i;
  struct timeval tv0, tv1;
  double t1, t2;

  for (i = 0; i < n; ++i)
    {
      double ti = (double) i / (n - 1.0);
      double xi = sin(2.0 * M_PI * ti) +
                  sin(2.0 * M_PI * 40.0 * ti);
      double ei = gsl_ran_gaussian(r, 0.2);

      gsl_vector_set(x, i, xi + ei);
    }

  /* y = rmedian(x) */
  gettimeofday(&tv0, NULL);
  gsl_filter_rmedian(etype, x, y_rmedian, w);
  gettimeofday(&tv1, NULL);
  t1 = TIMEDIFF(tv0, tv1);

  /* y = rmedian(x) with 2nd algorithm */
  gettimeofday(&tv0, NULL);
  gsl_filter_rmedian2(x, z_rmedian, w);
  gettimeofday(&tv1, NULL);
  t2 = TIMEDIFF(tv0, tv1);
  compare_vectors(GSL_DBL_EPSILON, z_rmedian, y_rmedian, "test_rmedian_sine 2nd algo");

  fprintf(stderr, "1st algo = %g [sec], 2nd algo = %g [sec]\n", t1, t2);

  /* z = rmedian(y) and check y = z */
  gsl_filter_rmedian(etype, y_rmedian, z_rmedian, w);
  compare_vectors(GSL_DBL_EPSILON, z_rmedian, y_rmedian, "test_rmedian_sine root sequence");

  gsl_filter_rmedian_free(w);
  gsl_vector_free(x);
  gsl_vector_free(y_rmedian);
  gsl_vector_free(z_rmedian);
}

void
test_rmedian(gsl_rng * rng_p)
{
  /* test root sequences */

  test_rmedian_root(GSL_FILTER_END_PADZERO, 1000, 3);
  test_rmedian_root(GSL_FILTER_END_PADZERO, 2000, 101);

  test_rmedian_root(GSL_FILTER_END_PADVALUE, 1000, 3);
  test_rmedian_root(GSL_FILTER_END_PADVALUE, 2000, 101);

  /* test random input */

  test_rmedian_random(GSL_FILTER_END_PADZERO, 10, 1, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADZERO, 100, 3, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADZERO, 1000, 3, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADZERO, 100, 1001, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADZERO, 5, 7, rng_p);

  test_rmedian_random(GSL_FILTER_END_PADVALUE, 10, 1, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADVALUE, 100, 3, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADVALUE, 1000, 3, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADVALUE, 100, 1001, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADVALUE, 5, 7, rng_p);

#if 0
  test_rmedian_sine(1000, 5, rng_p);
  test_rmedian_sine(5000, 71, rng_p);
  test_rmedian_sine(5000, 201, rng_p);
#elif 0
  test_rmedian_sine(GSL_FILTER_END_PADZERO, 500000, 6001, rng_p);
#endif
}
