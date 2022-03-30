#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>

#include <shtns.h>

#include "shtools.h"

#define TIME_DIFF(a,b) (((b).tv_sec + (b).tv_usec * 1.0e-6) - ((a).tv_sec + (a).tv_usec * 1.0e-6))

/* GSL */
double
proc_P_gsl(const size_t flags, const size_t lmax, const size_t n, const double x, double * Plm)
{
  struct timeval tv0, tv1;
  double dt;
  size_t i;

  gettimeofday(&tv0, NULL);

  gsl_sf_legendre_precompute(GSL_SF_LEGENDRE_SPHARM, lmax, flags, Plm);

  for (i = 0; i < n; ++i)
    gsl_sf_legendre_arrayx(GSL_SF_LEGENDRE_SPHARM, lmax, x, Plm);

  gettimeofday(&tv1, NULL);

  dt = TIME_DIFF(tv0, tv1);

  return dt;
}

/* SHTOOLS, l indexing */
double
proc_P_shtools(const size_t lmax, const size_t n, const double x, double * Plm)
{
  struct timeval tv0, tv1;
  size_t i;
  int lmaxi = (int) lmax;
  int csphase = 1;
  int cnorm = 1;
  int status;
  double dt;

  gettimeofday(&tv0, NULL);

  for (i = 0; i < n; ++i)
    shtools::PlmON(&Plm[0], lmaxi, x, &csphase, &cnorm, &status);

  gettimeofday(&tv1, NULL);

  dt = TIME_DIFF(tv0, tv1);

  /* deallocate memory */
  shtools::PlmON(&Plm[0], -1, x, &csphase, &cnorm, &status);

  return dt;
}

double
proc_P_shtns(shtns_cfg shtns, const size_t lmax, const size_t n, const double x, double * Plm)
{
  struct timeval tv0, tv1;
  size_t i;
  double dt;

  gettimeofday(&tv0, NULL);

  for (i = 0; i < n; ++i)
    {
      size_t idx = 0;
      int m;

      for (m = 0; m <= (int) lmax; ++m)
        {
          legendre_sphPlm_array(shtns, lmax, m, x, &Plm[idx]);
          idx += lmax + 1 - m;
        }
    }

  gettimeofday(&tv1, NULL);

  dt = TIME_DIFF(tv0, tv1);

  return dt;
}

int
main(int argc, char * argv[])
{
  const size_t lmax = 1500;
  const double x = -0.75;
  size_t eval_lmin = 2;
  size_t eval_lmax = 400;
  size_t n = 2000;
  const size_t plm_size = gsl_sf_legendre_array_n(lmax);
  double * Plm = (double *) malloc(plm_size * sizeof(double));
  shtns_cfg shtns = shtns_create((int) lmax, (int) lmax, 1, (shtns_norm) (sht_orthonormal | SHT_NO_CS_PHASE));
  size_t l;

  if (argc > 1)
    eval_lmin = (size_t) atoi(argv[1]);
  if (argc > 2)
    eval_lmax = (size_t) atoi(argv[2]);
  if (argc > 3)
    n = (size_t) atoi(argv[3]);

  l = 1;
  printf("Field %zu: ALF degree\n", l++);
  printf("Field %zu: GSL (l indexing) (seconds)\n", l++);
  printf("Field %zu: GSL (m indexing) (seconds)\n", l++);
  printf("Field %zu: SHTOOLS (seconds)\n", l++);
  printf("Field %zu: SHTns (seconds)\n", l++);

  for (l = eval_lmin; l <= eval_lmax; ++l)
    {
      double dt_gsll = proc_P_gsl(GSL_SF_LEGENDRE_FLG_INDEXL, l, n, x, Plm);
      double dt_gslm = proc_P_gsl(0, l, n, x, Plm);
      double dt_shtools = proc_P_shtools(l, n, x, Plm);
      double dt_shtns = proc_P_shtns(shtns, l, n, x, Plm);

      printf("%zu %.12e %.12e %.12e %.12e\n",
             l,
             dt_gsll / n,
             dt_gslm / n,
             dt_shtools / n,
             dt_shtns / n);
      fflush(stdout);
    }

  free(Plm);
  shtns_destroy(shtns);

  return 0;
}
