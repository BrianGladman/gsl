#include <config.h>
#include <stdio.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_test.h>

void rng_test (const gsl_rng_type * T, unsigned long int seed, unsigned int n,
	       unsigned long int result);
void rng_float_test (const gsl_rng_type * T);
void generic_rng_test (const gsl_rng_type * T);
void rng_state_test (const gsl_rng_type * T);
void rng_parallel_state_test (const gsl_rng_type * T);
int rng_max_test (gsl_rng * r, unsigned long int *kmax, unsigned long int ran_max) ;
int rng_min_test (gsl_rng * r, unsigned long int *kmin, unsigned long int ran_min, unsigned long int ran_max) ;
int rng_sum_test (gsl_rng * r, double *sigma);
int rng_bin_test (gsl_rng * r, double *sigma);

#define N  10000
#define N2 200000

int
main (void)
{
  gsl_rng_env_setup ();

  /* specific tests of known results for 10000 iterations with seed = 1 */

  rng_test (gsl_rng_rand, 1, 10000, 1910041713);
  rng_test (gsl_rng_randu, 1, 10000, 1623524161);
  rng_test (gsl_rng_cmrg, 1, 10000, 719452880);
  rng_test (gsl_rng_minstd, 1, 10000, 1043618065);
  rng_test (gsl_rng_mrg, 1, 10000, 2064828650);
  rng_test (gsl_rng_taus, 1, 10000, 2733957125UL);
  rng_test (gsl_rng_transputer, 1, 10000, 1244127297UL);
  rng_test (gsl_rng_vax, 1, 10000, 3051034865UL);

  /* FIXME: the ranlux tests below were made by running the fortran code and
     getting the expected value from that. An analytic calculation
     would be preferable. */

  rng_test (gsl_rng_ranlux, 314159265, 10000, 12077992);
  rng_test (gsl_rng_ranlux389, 314159265, 10000, 165942);

  /* FIXME: the tests below were made by running the original code in
     the ../random directory and getting the expected value from
     that. An analytic calculation would be preferable. */

  rng_test (gsl_rng_slatec, 1, 10000, 45776);
  rng_test (gsl_rng_uni, 1, 10000, 9214);
  rng_test (gsl_rng_uni32, 1, 10000, 1155229825);
  rng_test (gsl_rng_zuf, 1, 10000, 3970);

  /* The tests below were made by running the original code and
     getting the expected value from that. An analytic calculation
     would be preferable. */

  rng_test (gsl_rng_r250, 1, 10000, 1100653588);
  rng_test (gsl_rng_mt19937, 4357, 1000, 1309179303);
  rng_test (gsl_rng_tt800, 0, 10000, 2856609219UL);

  rng_test (gsl_rng_ran0, 0, 10000, 1115320064);
  rng_test (gsl_rng_ran1, 0, 10000, 1491066076);
  rng_test (gsl_rng_ran2, 0, 10000, 1701364455);
  rng_test (gsl_rng_ran3, 0, 10000, 186340785);

  rng_test (gsl_rng_ranmar, 1, 10000, 14428370);

  rng_test (gsl_rng_rand48, 0, 10000, 0xDE095043UL);
  rng_test (gsl_rng_rand48, 1, 10000, 0xEDA54977UL);

  rng_test (gsl_rng_random_glibc2, 0, 10000, 1908609430);
  rng_test (gsl_rng_random8_glibc2, 0, 10000, 1910041713);
  rng_test (gsl_rng_random32_glibc2, 0, 10000, 1587395585);
  rng_test (gsl_rng_random64_glibc2, 0, 10000, 52848624);
  rng_test (gsl_rng_random128_glibc2, 0, 10000, 1908609430);
  rng_test (gsl_rng_random256_glibc2, 0, 10000, 179943260);

  rng_test (gsl_rng_random_bsd, 0, 10000, 1457025928);
  rng_test (gsl_rng_random8_bsd, 0, 10000, 1910041713);
  rng_test (gsl_rng_random32_bsd, 0, 10000, 1663114331);
  rng_test (gsl_rng_random64_bsd, 0, 10000, 864469165);
  rng_test (gsl_rng_random128_bsd, 0, 10000, 1457025928);
  rng_test (gsl_rng_random256_bsd, 0, 10000, 1216357476);

  rng_test (gsl_rng_random_libc5, 0, 10000, 428084942);
  rng_test (gsl_rng_random8_libc5, 0, 10000, 1910041713);
  rng_test (gsl_rng_random32_libc5, 0, 10000, 1967452027);
  rng_test (gsl_rng_random64_libc5, 0, 10000, 2106639801);
  rng_test (gsl_rng_random128_libc5, 0, 10000, 428084942);
  rng_test (gsl_rng_random256_libc5, 0, 10000, 116367984);

  rng_test (gsl_rng_ranf, 0, 10000, 2152890433UL);
  rng_test (gsl_rng_ranf, 2, 10000, 339327233);

  /* Test constant relationship between int and double functions */

  rng_float_test (gsl_rng_slatec);
  rng_float_test (gsl_rng_gfsr4);
  rng_float_test (gsl_rng_cmrg);
  rng_float_test (gsl_rng_minstd);
  rng_float_test (gsl_rng_mrg);
  rng_float_test (gsl_rng_mt19937);
  rng_float_test (gsl_rng_r250);
  rng_float_test (gsl_rng_ran0);
  rng_float_test (gsl_rng_ran1);
  rng_float_test (gsl_rng_ran2);
  rng_float_test (gsl_rng_ran3);
  rng_float_test (gsl_rng_rand);
  rng_float_test (gsl_rng_rand48);
  rng_float_test (gsl_rng_random8_bsd);
  rng_float_test (gsl_rng_random8_glibc2);
  rng_float_test (gsl_rng_random8_libc5);
  rng_float_test (gsl_rng_random128_bsd);
  rng_float_test (gsl_rng_random128_glibc2);
  rng_float_test (gsl_rng_random128_libc5);
  rng_float_test (gsl_rng_random256_bsd);
  rng_float_test (gsl_rng_random256_glibc2);
  rng_float_test (gsl_rng_random256_libc5);
  rng_float_test (gsl_rng_random32_bsd);
  rng_float_test (gsl_rng_random32_glibc2);
  rng_float_test (gsl_rng_random32_libc5);
  rng_float_test (gsl_rng_random64_bsd);
  rng_float_test (gsl_rng_random64_glibc2);
  rng_float_test (gsl_rng_random64_libc5);
  rng_float_test (gsl_rng_random_bsd);
  rng_float_test (gsl_rng_random_glibc2);
  rng_float_test (gsl_rng_random_libc5);
  rng_float_test (gsl_rng_randu);
  rng_float_test (gsl_rng_ranf);
  rng_float_test (gsl_rng_ranlux);
  rng_float_test (gsl_rng_ranlux389);
  rng_float_test (gsl_rng_ranmar);
  rng_float_test (gsl_rng_taus);
  rng_float_test (gsl_rng_transputer);
  rng_float_test (gsl_rng_tt800);
  rng_float_test (gsl_rng_uni);
  rng_float_test (gsl_rng_uni32);
  rng_float_test (gsl_rng_vax);
  rng_float_test (gsl_rng_zuf);

  /* Test save/restore functions */

  rng_state_test (gsl_rng_slatec);
  rng_state_test (gsl_rng_cmrg);
  rng_state_test (gsl_rng_gfsr4);
  rng_state_test (gsl_rng_minstd);
  rng_state_test (gsl_rng_mrg);
  rng_state_test (gsl_rng_mt19937);
  rng_state_test (gsl_rng_r250);
  rng_state_test (gsl_rng_ran0);
  rng_state_test (gsl_rng_ran1);
  rng_state_test (gsl_rng_ran2);
  rng_state_test (gsl_rng_ran3);
  rng_state_test (gsl_rng_rand);
  rng_state_test (gsl_rng_rand48);
  rng_state_test (gsl_rng_random8_bsd);
  rng_state_test (gsl_rng_random8_glibc2);
  rng_state_test (gsl_rng_random8_libc5);
  rng_state_test (gsl_rng_random128_bsd);
  rng_state_test (gsl_rng_random128_glibc2);
  rng_state_test (gsl_rng_random128_libc5);
  rng_state_test (gsl_rng_random256_bsd);
  rng_state_test (gsl_rng_random256_glibc2);
  rng_state_test (gsl_rng_random256_libc5);
  rng_state_test (gsl_rng_random32_bsd);
  rng_state_test (gsl_rng_random32_glibc2);
  rng_state_test (gsl_rng_random32_libc5);
  rng_state_test (gsl_rng_random64_bsd);
  rng_state_test (gsl_rng_random64_glibc2);
  rng_state_test (gsl_rng_random64_libc5);
  rng_state_test (gsl_rng_random_bsd);
  rng_state_test (gsl_rng_random_glibc2);
  rng_state_test (gsl_rng_random_libc5);
  rng_state_test (gsl_rng_randu);
  rng_state_test (gsl_rng_ranf);
  rng_state_test (gsl_rng_ranlux);
  rng_state_test (gsl_rng_ranlux389);
  rng_state_test (gsl_rng_ranmar);
  rng_state_test (gsl_rng_taus);
  rng_state_test (gsl_rng_transputer);
  rng_state_test (gsl_rng_tt800);
  rng_state_test (gsl_rng_uni);
  rng_state_test (gsl_rng_uni32);
  rng_state_test (gsl_rng_vax);
  rng_state_test (gsl_rng_zuf);

  rng_parallel_state_test (gsl_rng_slatec);
  rng_parallel_state_test (gsl_rng_cmrg);
  rng_parallel_state_test (gsl_rng_gfsr4);
  rng_parallel_state_test (gsl_rng_minstd);
  rng_parallel_state_test (gsl_rng_mrg);
  rng_parallel_state_test (gsl_rng_mt19937);
  rng_parallel_state_test (gsl_rng_r250);
  rng_parallel_state_test (gsl_rng_ran0);
  rng_parallel_state_test (gsl_rng_ran1);
  rng_parallel_state_test (gsl_rng_ran2);
  rng_parallel_state_test (gsl_rng_ran3);
  rng_parallel_state_test (gsl_rng_rand);
  rng_parallel_state_test (gsl_rng_rand48);
  rng_parallel_state_test (gsl_rng_random8_bsd);
  rng_parallel_state_test (gsl_rng_random8_glibc2);
  rng_parallel_state_test (gsl_rng_random8_libc5);
  rng_parallel_state_test (gsl_rng_random128_bsd);
  rng_parallel_state_test (gsl_rng_random128_glibc2);
  rng_parallel_state_test (gsl_rng_random128_libc5);
  rng_parallel_state_test (gsl_rng_random256_bsd);
  rng_parallel_state_test (gsl_rng_random256_glibc2);
  rng_parallel_state_test (gsl_rng_random256_libc5);
  rng_parallel_state_test (gsl_rng_random32_bsd);
  rng_parallel_state_test (gsl_rng_random32_glibc2);
  rng_parallel_state_test (gsl_rng_random32_libc5);
  rng_parallel_state_test (gsl_rng_random64_bsd);
  rng_parallel_state_test (gsl_rng_random64_glibc2);
  rng_parallel_state_test (gsl_rng_random64_libc5);
  rng_parallel_state_test (gsl_rng_random_bsd);
  rng_parallel_state_test (gsl_rng_random_glibc2);
  rng_parallel_state_test (gsl_rng_random_libc5);
  rng_parallel_state_test (gsl_rng_randu);
  rng_parallel_state_test (gsl_rng_ranf);
  rng_parallel_state_test (gsl_rng_ranlux);
  rng_parallel_state_test (gsl_rng_ranlux389);
  rng_parallel_state_test (gsl_rng_ranmar);
  rng_parallel_state_test (gsl_rng_taus);
  rng_parallel_state_test (gsl_rng_transputer);
  rng_parallel_state_test (gsl_rng_tt800);
  rng_parallel_state_test (gsl_rng_uni);
  rng_parallel_state_test (gsl_rng_uni32);
  rng_parallel_state_test (gsl_rng_vax);
  rng_parallel_state_test (gsl_rng_zuf);

  /* generic statistical tests (these are just to make sure that we
     don't get any crazy results back from the generator, i.e. they
     aren't a test of the algorithm, just the implementation) */

  generic_rng_test (gsl_rng_slatec);
  generic_rng_test (gsl_rng_cmrg);
  generic_rng_test (gsl_rng_gfsr4);
  generic_rng_test (gsl_rng_minstd);
  generic_rng_test (gsl_rng_mrg);
  generic_rng_test (gsl_rng_mt19937);
  generic_rng_test (gsl_rng_r250);
  generic_rng_test (gsl_rng_ran0);
  generic_rng_test (gsl_rng_ran1);
  generic_rng_test (gsl_rng_ran2);
  generic_rng_test (gsl_rng_ran3);
  generic_rng_test (gsl_rng_rand);
  generic_rng_test (gsl_rng_rand48);
  generic_rng_test (gsl_rng_random8_bsd);
  generic_rng_test (gsl_rng_random8_glibc2);
  generic_rng_test (gsl_rng_random8_libc5);
  generic_rng_test (gsl_rng_random128_bsd);
  generic_rng_test (gsl_rng_random128_glibc2);
  generic_rng_test (gsl_rng_random128_libc5);
  generic_rng_test (gsl_rng_random256_bsd);
  generic_rng_test (gsl_rng_random256_glibc2);
  generic_rng_test (gsl_rng_random256_libc5);
  generic_rng_test (gsl_rng_random32_bsd);
  generic_rng_test (gsl_rng_random32_glibc2);
  generic_rng_test (gsl_rng_random32_libc5);
  generic_rng_test (gsl_rng_random64_bsd);
  generic_rng_test (gsl_rng_random64_glibc2);
  generic_rng_test (gsl_rng_random64_libc5);
  generic_rng_test (gsl_rng_random_bsd);
  generic_rng_test (gsl_rng_random_glibc2);
  generic_rng_test (gsl_rng_random_libc5);
  generic_rng_test (gsl_rng_randu);
  generic_rng_test (gsl_rng_ranf);
  generic_rng_test (gsl_rng_ranlux);
  generic_rng_test (gsl_rng_ranlux389);
  generic_rng_test (gsl_rng_ranmar);
  generic_rng_test (gsl_rng_taus);
  generic_rng_test (gsl_rng_transputer);
  generic_rng_test (gsl_rng_tt800);
  generic_rng_test (gsl_rng_uni);
  generic_rng_test (gsl_rng_uni32);
  generic_rng_test (gsl_rng_vax);
  generic_rng_test (gsl_rng_zuf);

  return gsl_test_summary ();
}


void
rng_test (const gsl_rng_type * T, unsigned long int seed, unsigned int n,
	  unsigned long int result)
{
  gsl_rng *r = gsl_rng_alloc (T);
  unsigned int i;
  unsigned long int k = 0;
  int status;

  if (seed != 0)
    {
      gsl_rng_set (r, seed);
    }

  for (i = 0; i < n; i++)
    {
      k = gsl_rng_get (r);
    }

  status = (k != result);
  gsl_test (status, "%s, %u steps (%u observed vs %u expected)",
	    gsl_rng_name (r), n, k, result);

  gsl_rng_free (r);
}

void
rng_float_test (const gsl_rng_type * T)
{
  gsl_rng *ri = gsl_rng_alloc (T);
  gsl_rng *rf = gsl_rng_alloc (T);

  double u, c ; 
  unsigned int i;
  unsigned long int k = 0;
  int status = 0 ;

  do 
    {
      k = gsl_rng_get (ri);
      u = gsl_rng_get (rf);
    } 
  while (k == 0) ;

  c = k / u ;
  for (i = 0; i < N2; i++)
    {
      k = gsl_rng_get (ri);
      u = gsl_rng_get (rf);
      if (c*k != u)
	{
	  status = 1 ;
	  break ;
	}
    }

  gsl_test (status, "%s, ratio of int to double (%g observed vs %g expected)",
	    gsl_rng_name (ri), c, k/u);

  gsl_rng_free (ri);
  gsl_rng_free (rf);
}


void
rng_state_test (const gsl_rng_type * T)
{
  unsigned long int test_a[N], test_b[N];

  int i;

  gsl_rng *r = gsl_rng_alloc (T);
  gsl_rng *r_save = gsl_rng_alloc (T);

  for (i = 0; i < N; ++i)
    {
      gsl_rng_get (r);	/* throw away N iterations */
    }

  gsl_rng_cpy (r_save, r);	/* save the intermediate state */

  for (i = 0; i < N; ++i)
    {
      test_a[i] = gsl_rng_get (r);
    }

  gsl_rng_cpy (r, r_save);	/* restore the intermediate state */
  gsl_rng_free (r_save);

  for (i = 0; i < N; ++i)
    {
      test_b[i] = gsl_rng_get (r);
    }

  {
    int status = 0;
    for (i = 0; i < N; ++i)
      {
	status |= (test_b[i] != test_a[i]);
      }
    gsl_test (status, "%s, random number state consistency",
	      gsl_rng_name (r));
  }

  gsl_rng_free (r);
}


void
rng_parallel_state_test (const gsl_rng_type * T)
{
  unsigned long int test_a[N], test_b[N];
  unsigned long int test_c[N], test_d[N];
  double test_e[N], test_f[N];

  int i;

  gsl_rng *r1 = gsl_rng_alloc (T);
  gsl_rng *r2 = gsl_rng_alloc (T);

  for (i = 0; i < N; ++i)
    {
      gsl_rng_get (r1);		/* throw away N iterations */
    }

  gsl_rng_cpy (r2, r1);		/* save the intermediate state */

  for (i = 0; i < N; ++i)
    {
      /* check that there is no hidden state intermixed between r1 and r2 */
      test_a[i] = gsl_rng_get (r1);	
      test_b[i] = gsl_rng_get (r2);
      test_c[i] = gsl_rng_uniform_int (r1, 1234);	
      test_d[i] = gsl_rng_uniform_int (r2, 1234);
      test_e[i] = gsl_rng_uniform (r1);	
      test_f[i] = gsl_rng_uniform (r2);
    }

  {
    int status = 0;
    for (i = 0; i < N; ++i)
      {
	status |= (test_b[i] != test_a[i]);
	status |= (test_c[i] != test_d[i]);
	status |= (test_e[i] != test_f[i]);
      }
    gsl_test (status, "%s, parallel random number state consistency",
	      gsl_rng_name (r1));
  }

  gsl_rng_free (r1);
  gsl_rng_free (r2);

}

void
generic_rng_test (const gsl_rng_type * T)
{
  gsl_rng *r = gsl_rng_alloc (T);
  const char *name = gsl_rng_name (r);
  unsigned long int kmax = 0, kmin = 1000;
  double sigma = 0;
  const unsigned long int ran_max = gsl_rng_max (r);
  const unsigned long int ran_min = gsl_rng_min (r);

  int status = rng_max_test (r, &kmax, ran_max);

  gsl_test (status,
	    "%s, observed vs theoretical maximum (%lu vs %lu)",
	    name, kmax, ran_max);

  status = rng_min_test (r, &kmin, ran_min, ran_max);

  gsl_test (status,
	    "%s, observed vs theoretical minimum (%lu vs %lu)",
	    name, kmin, ran_min);

  status = rng_sum_test (r, &sigma);

  gsl_test (status,
	    "%s, sum test within acceptable sigma (observed %.2g sigma)",
	    name, sigma);

  status = rng_bin_test (r, &sigma);

  gsl_test (status,
	    "%s, bin test within acceptable chisq (observed %.2g sigma)",
	    name, sigma);

  gsl_rng_set (r, 1);	/* set seed to 1 */
  status = rng_max_test (r, &kmax, ran_max);

  gsl_rng_set (r, 1);	/* set seed to 1 */
  status |= rng_min_test (r, &kmin, ran_min, ran_max);

  gsl_rng_set (r, 1);	/* set seed to 1 */
  status |= rng_sum_test (r, &sigma);

  gsl_rng_set (r, 12345);	/* set seed to a "typical" value */
  status |= rng_max_test (r, &kmax, ran_max);

  gsl_rng_set (r, 12345);	/* set seed to a "typical" value */
  status |= rng_min_test (r, &kmin, ran_min, ran_max);

  gsl_rng_set (r, 12345);	/* set seed to a "typical" value */
  status |= rng_sum_test (r, &sigma);

  gsl_test (status, "%s, maximum and sum tests for non-default seeds", name);

  gsl_rng_free (r);
}

int
rng_max_test (gsl_rng * r, unsigned long int *kmax, unsigned long int ran_max)
{
  unsigned long int actual_uncovered;
  double expect_uncovered;
  int status;
  unsigned long int max = 0;
  int i;

  for (i = 0; i < N2; ++i)
    {
      unsigned long int k = gsl_rng_get (r);
      if (k > max)
	max = k;
    }

  *kmax = max;

  actual_uncovered = ran_max - max;
  expect_uncovered = (double) ran_max / (double) N2;

  status = (max > ran_max) || (actual_uncovered > 5 * expect_uncovered) ;

  return status;
}

int
rng_min_test (gsl_rng * r, unsigned long int *kmin, 
	      unsigned long int ran_min, unsigned long int ran_max)
{
  unsigned long int actual_uncovered;
  double expect_uncovered;
  int status;
  unsigned long int min = 1e9;
  int i;

  for (i = 0; i < N2; ++i)
    {
      unsigned long int k = gsl_rng_get (r);
      if (k < min)
	min = k;
    }

  *kmin = min;

  actual_uncovered = min - ran_min;
  expect_uncovered = (double) ran_max / (double) N2;

  status = (min < ran_min) || (actual_uncovered > 5 * expect_uncovered);

  return status;
}

int
rng_sum_test (gsl_rng * r, double *sigma)
{
  double sum = 0;
  int i, status;

  for (i = 0; i < N2; ++i)
    {
      double x = gsl_rng_uniform (r) - 0.5;
      sum += x;
    }

  sum /= N2;

  /* expect the average to have a variance of 1/(12 n) */

  *sigma = sum * sqrt (12.0 * N2);

  /* more than 3 sigma is an error */
  
  status = (fabs (*sigma) > 3 || fabs(*sigma) < 0.003);

  if (status) {
      fprintf(stderr,"sum=%g, sigma=%g\n",sum,*sigma);
  }

  return status;
}

#define BINS 17
#define EXTRA 10
int
rng_bin_test (gsl_rng * r, double *sigma)
{
  int count[BINS+EXTRA];
  double chisq = 0;
  int i, status;

  for (i = 0; i < BINS+EXTRA; i++)
      count[i] = 0 ;


  for (i = 0; i < N2; i++)
    {
      int j = gsl_rng_uniform_int (r, BINS);
      count[j]++ ;
    }

  chisq = 0 ;
  for (i = 0; i < BINS; i++)
    {
      double x = (double)N2/(double)BINS ;
      double d = (count[i] - x) ;
      chisq += (d*d) / x;
    }

  *sigma = sqrt(chisq/BINS) ;

  /* more than 3 sigma is an error */
  
  status = (fabs (*sigma) > 3 || fabs(*sigma) < 0.003);

  for (i = BINS; i < BINS+EXTRA; i++)
    {
      if (count[i] != 0)
	{
	  status = 1 ;
	  gsl_test (status, 
		    "%s, wrote outside range in bin test "
		    "(%d observed vs %d expected)",
		    gsl_rng_name(r), i, BINS - 1);
	}
    }

  return status;
}

