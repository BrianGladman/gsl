#include <stdio.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_test.h>

void rng_test (const gsl_rng_type * T, unsigned long int seed, unsigned int n,
	       unsigned long int result);
void generic_rng_test (const gsl_rng_type * T);
void rng_state_test (const gsl_rng_type * T);
void rng_parallel_state_test (const gsl_rng_type * T);
int rng_max_test (gsl_rng * r, unsigned long int *kmax, unsigned long int ran_max) ;
int rng_min_test (gsl_rng * r, unsigned long int *kmin, unsigned long int ran_min, unsigned long int ran_max) ;
int rng_sum_test (gsl_rng * r, double *sigma);

#define N  10000
#define N2 200000

int
main (void)
{
  gsl_rng_env_setup ();

  /* specific tests of known results for 10000 iterations with seed = 1 */

  rng_test (gsl_rng_bsdrand, 1, 10000, 1910041713);
  rng_test (gsl_rng_randu, 1, 10000, 1623524161);
  rng_test (gsl_rng_cmrg, 1, 10000, 719452880);
  rng_test (gsl_rng_minstd, 1, 10000, 1043618065);
  rng_test (gsl_rng_mrg, 1, 10000, 2064828650);
  rng_test (gsl_rng_taus, 1, 10000, 2733957125UL);
  rng_test (gsl_rng_tds, 1, 10000, 1244127297UL);
  rng_test (gsl_rng_vax, 1, 10000, 3051034865UL);

  /* FIXME: the ranlux tests below were made by running the fortran code and
     getting the expected value from that. An analytic calculation
     would be preferable. */

  rng_test (gsl_rng_ranlux, 314159265, 10000, 12077992);
  rng_test (gsl_rng_ranlux389, 314159265, 10000, 165942);

  /* FIXME: the tests below were made by running the original code in
     the ../random directory and getting the expected value from
     that. An analytic calculation would be preferable. */

  rng_test (gsl_rng_rand, 1, 10000, 45776);
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

  /* Test save/restore functions */

  rng_state_test (gsl_rng_bsdrand);
  rng_state_test (gsl_rng_cmrg);
  rng_state_test (gsl_rng_minstd);
  rng_state_test (gsl_rng_mrg);
  rng_state_test (gsl_rng_mt19937);
  rng_state_test (gsl_rng_r250);
  rng_state_test (gsl_rng_ran0);
  rng_state_test (gsl_rng_ran1);
  rng_state_test (gsl_rng_ran2);
  rng_state_test (gsl_rng_ran3);
  rng_state_test (gsl_rng_rand);
  rng_state_test (gsl_rng_randu);
  rng_state_test (gsl_rng_ranlux);
  rng_state_test (gsl_rng_ranlux389);
  rng_state_test (gsl_rng_ranmar);
  rng_state_test (gsl_rng_taus);
  rng_state_test (gsl_rng_tds);
  rng_state_test (gsl_rng_tt800);
  rng_state_test (gsl_rng_uni);
  rng_state_test (gsl_rng_uni32);
  rng_state_test (gsl_rng_vax);
  rng_state_test (gsl_rng_zuf);

  rng_parallel_state_test (gsl_rng_bsdrand);
  rng_parallel_state_test (gsl_rng_cmrg);
  rng_parallel_state_test (gsl_rng_minstd);
  rng_parallel_state_test (gsl_rng_mrg);
  rng_parallel_state_test (gsl_rng_mt19937);
  rng_parallel_state_test (gsl_rng_r250);
  rng_parallel_state_test (gsl_rng_ran0);
  rng_parallel_state_test (gsl_rng_ran1);
  rng_parallel_state_test (gsl_rng_ran2);
  rng_parallel_state_test (gsl_rng_ran3);
  rng_parallel_state_test (gsl_rng_rand);
  rng_parallel_state_test (gsl_rng_randu);
  rng_parallel_state_test (gsl_rng_ranlux);
  rng_parallel_state_test (gsl_rng_ranlux389);
  rng_parallel_state_test (gsl_rng_ranmar);
  rng_parallel_state_test (gsl_rng_taus);
  rng_parallel_state_test (gsl_rng_tds);
  rng_parallel_state_test (gsl_rng_tt800);
  rng_parallel_state_test (gsl_rng_uni);
  rng_parallel_state_test (gsl_rng_uni32);
  rng_parallel_state_test (gsl_rng_vax);
  rng_parallel_state_test (gsl_rng_zuf);

  /* generic statistical tests (these are just to make sure that we
     don't get any crazy results back from the generator, i.e. they
     aren't a test of the algorithm, just the implementation) */

  generic_rng_test (gsl_rng_bsdrand);
  generic_rng_test (gsl_rng_cmrg);
  generic_rng_test (gsl_rng_minstd);
  generic_rng_test (gsl_rng_mrg);
  generic_rng_test (gsl_rng_mt19937);
  generic_rng_test (gsl_rng_r250);
  generic_rng_test (gsl_rng_ran0);
  generic_rng_test (gsl_rng_ran1);
  generic_rng_test (gsl_rng_ran2);
  generic_rng_test (gsl_rng_ran3);
  generic_rng_test (gsl_rng_rand);
  generic_rng_test (gsl_rng_randu);
  generic_rng_test (gsl_rng_ranlux);
  generic_rng_test (gsl_rng_ranlux389);
  generic_rng_test (gsl_rng_ranmar);
  generic_rng_test (gsl_rng_taus);
  generic_rng_test (gsl_rng_tds);
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
  gsl_test (status, "%s, %u iterations (%u observed vs %u expected)",
	    gsl_rng_name (r), n, k, result);

  gsl_rng_free (r);
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
    }

  {
    int status = 0;
    for (i = 0; i < N; ++i)
      {
	status |= (test_b[i] != test_a[i]);
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
  status = (fabs (*sigma) > 3);		/* more than 3 sigma is an error */

  return status;
}
