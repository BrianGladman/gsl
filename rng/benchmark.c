#include <config.h>
#include <time.h>
#include <stdio.h>
#include <gsl_rng.h>
#include <gsl_errno.h>

void benchmark (const gsl_rng_type * T);

#define N  1000000
int sum;

int
main (void)
{
  benchmark (gsl_rng_bsdrand);
  benchmark (gsl_rng_randu);
  benchmark (gsl_rng_cmrg);
  benchmark (gsl_rng_minstd);
  benchmark (gsl_rng_mrg);
  benchmark (gsl_rng_mt19937);
  benchmark (gsl_rng_rand);
  benchmark (gsl_rng_ranlux);
  benchmark (gsl_rng_ranlux389);
  benchmark (gsl_rng_ranmar);
  benchmark (gsl_rng_r250);
  benchmark (gsl_rng_taus);
  benchmark (gsl_rng_transputer);
  benchmark (gsl_rng_tt800);
  benchmark (gsl_rng_uni);
  benchmark (gsl_rng_uni32);
  benchmark (gsl_rng_vax);
  benchmark (gsl_rng_zuf);
  return 0;
}

void
benchmark (const gsl_rng_type * T)
{
  int start, end;
  int i = 0, d = 0 ;
  double t1, t2;

  gsl_rng *r = gsl_rng_alloc (T);

  start = clock ();
  do
    {
      int j;
      for (j = 0; j < N; j++)
	sum += gsl_rng_get (r);

      i += N;
      end = clock ();
    }
  while (end < start + CLOCKS_PER_SEC/10);

  t1 = (end - start) / (double) CLOCKS_PER_SEC;

  start = clock ();
  do
    {
      int j;
      for (j = 0; j < N; j++)
	sum += gsl_rng_uniform (r);

      d += N;
      end = clock ();
    }
  while (end < start + CLOCKS_PER_SEC/10);

  t2 = (end - start) / (double) CLOCKS_PER_SEC;


  printf ("%15s %.2e random ints (%.2f%% doubles) per second\n",
	  gsl_rng_name (r),
	  i / t1, 100*(d/t2)/(i/t1));

  gsl_rng_free (r);
}
