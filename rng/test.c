#include <stdio.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_test.h>

void rng_test (const gsl_rng_type * T, unsigned int seed, unsigned int n, 
	       unsigned long int result);
void generic_rng_test (const gsl_rng_type * T);

int
main (void)
{
  gsl_rng_env_setup() ;

  /* specific tests of known results for 10000 iterations with seed = 1*/

  rng_test (gsl_rng_minstd,1,10000,1043618065); 
  rng_test (gsl_rng_bad_rand,1,10000,1910041713); 
  rng_test (gsl_rng_bad_randu,1,10000,1623524161); 
  rng_test (gsl_rng_cmrg,1,10000,1477798470); 
  rng_test (gsl_rng_mrg,1,10000,1711374253); 
  rng_test (gsl_rng_taus,1,10000,676146779);
  rng_test (gsl_rng_vax,1,10000,3051034865UL); 


  rng_test (gsl_rng_ranlux,314159265,10000,12077992); 
  rng_test (gsl_rng_ranlux389,314159265,10000,165942); 

  /* FIXME: the ranlux test was made by running the fortran code and
     getting the expected value from that. An analytic calculation
     would be preferable. */

  /* generic statistical tests */

  generic_rng_test (gsl_rng_cmrg);
  generic_rng_test (gsl_rng_mrg);
  generic_rng_test (gsl_rng_minstd);
  generic_rng_test (gsl_rng_rand);
  generic_rng_test (gsl_rng_taus);
  generic_rng_test (gsl_rng_uni);
  generic_rng_test (gsl_rng_uni32);
  generic_rng_test (gsl_rng_zuf);
  generic_rng_test (gsl_rng_ranlux);

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
      || (2 * expected_uncovered < actual_uncovered && actual_uncovered > 1);

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


void 
rng_test (const gsl_rng_type * T, unsigned int seed, unsigned int n, 
	  unsigned long int result)
{
  gsl_rng * r = gsl_rng_alloc (T);
  unsigned int i ;
  unsigned long int k = 0;
  int status;

  if (seed != 1) {
    gsl_rng_set(r,seed) ;
  }

  for (i = 0; i < n; i++)
    {
      k = gsl_rng_get(r) ;
    }

  status = (k != result) ;
  gsl_test(status, "%s, %u iterations (%u observed vs %u expected)",
	   gsl_rng_name(r), n, k, result) ;
 
  gsl_rng_free(r) ;
}
