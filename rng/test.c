#include <stdio.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_test.h>

void generic_rng_test (const gsl_rng_type * T);

int
main (void)
{
  gsl_rng_env_setup() ;
  generic_rng_test (gsl_rng_cmrg);
  generic_rng_test (gsl_rng_mrg);
  generic_rng_test (gsl_rng_rand);
  generic_rng_test (gsl_rng_taus);
  generic_rng_test (gsl_rng_uni);
  generic_rng_test (gsl_rng_uni32);
  generic_rng_test (gsl_rng_zuf);

  return gsl_test_summary ();
}

void
generic_rng_test (const gsl_rng_type * T)
{
  int n = 1000000;

  gsl_rng *r = gsl_rng_alloc (T);
  const unsigned long int ran_max = gsl_rng_max (r);
  const char *name = gsl_rng_name (r);

  {
    unsigned long int actual_uncovered;
    double expected_uncovered;
    int status;
    unsigned long int kmax = 0;
    int i;

    for (i = 0; i < n; ++i)
      {
	unsigned long int k = gsl_rng_get (r);
	if (k > kmax)
	  kmax = k;
      }

    actual_uncovered = ran_max - kmax;
    expected_uncovered = (double) ran_max / (double) n;

    /* The uni generator never actually reaches its ran_max in practice,
       due to the way the initial state is generated from the seed.
       Thus it only hits 32766 instead of 32767. 

       We'll let it pass by checking if the observed max is just 1 below
       the theoretical max.  */

    status = (kmax > ran_max)
      || (expected_uncovered < actual_uncovered && actual_uncovered > 1);

    gsl_test (status,
	      "%s, observed vs theoretical maximum (%lu vs %lu)",
	      name, kmax, ran_max, n);
  };

  {
    double sum = 0, sigma;
    int i, status;

    for (i = 0; i < n; ++i)
      sum += gsl_rng_get_uni (r);

    sum /= n;

    /* expect sum to have variance == n*(1/12), so average should have
       variance == 1/(12*n) */

    sigma = (sum - 0.5) * sqrt (12.0 * n);

    status = (fabs (sigma) > 3);	/* more than 3 sigma is an error */
    gsl_test (status,
	      "%s, sum test within acceptable sigma (observed %g sigma)",
	      name, sigma);
  }


}

