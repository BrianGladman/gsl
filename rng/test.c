#include <stdio.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_test.h>

void rng_test (const gsl_rng_type * T, unsigned long int seed, unsigned int n, 
	       unsigned long int result);
void generic_rng_test (const gsl_rng_type * T);
void rng_state_test (const gsl_rng_type * T);
void rng_parallel_state_test (const gsl_rng_type * T);

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

  /* FIXME: the ranlux tests below were made by running the fortran code and
     getting the expected value from that. An analytic calculation
     would be preferable. */

  rng_test (gsl_rng_ranlux,314159265,10000,12077992); 
  rng_test (gsl_rng_ranlux389,314159265,10000,165942); 

  /* FIXME: the tests below were made by running the original code in
     the ../random directory and getting the expected value from
     that. An analytic calculation would be preferable. */

  rng_test (gsl_rng_rand,1,10000,45776); 
  rng_test (gsl_rng_uni,1,10000,9214); 
  rng_test (gsl_rng_uni32,1,10000,1155229825); 
  rng_test (gsl_rng_zuf,1,10000,3970); 

  rng_test (gsl_rng_mt19937,4357,1000,1309179303); 
  rng_test (gsl_rng_tt800,0,10000,2856609219UL); 

  /* Test save/restore functions */

  rng_state_test (gsl_rng_cmrg);
  rng_state_test (gsl_rng_mrg);
  rng_state_test (gsl_rng_minstd);
  rng_state_test (gsl_rng_rand);
  rng_state_test (gsl_rng_taus);
  rng_state_test (gsl_rng_uni);
  rng_state_test (gsl_rng_uni32);
  rng_state_test (gsl_rng_zuf);
  rng_state_test (gsl_rng_ranlux);
  rng_state_test (gsl_rng_mt19937);
  rng_state_test (gsl_rng_tt800);

  rng_parallel_state_test (gsl_rng_cmrg);
  rng_parallel_state_test (gsl_rng_mrg);
  rng_parallel_state_test (gsl_rng_minstd);
  rng_parallel_state_test (gsl_rng_rand);
  rng_parallel_state_test (gsl_rng_taus);
  rng_parallel_state_test (gsl_rng_uni);
  rng_parallel_state_test (gsl_rng_uni32);
  rng_parallel_state_test (gsl_rng_zuf);
  rng_parallel_state_test (gsl_rng_ranlux);
  rng_parallel_state_test (gsl_rng_mt19937);
  rng_parallel_state_test (gsl_rng_tt800);

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
  generic_rng_test (gsl_rng_mt19937);
  generic_rng_test (gsl_rng_tt800);

  return gsl_test_summary ();
}


void 
rng_test (const gsl_rng_type * T, unsigned long int seed, unsigned int n, 
	  unsigned long int result)
{
  gsl_rng * r = gsl_rng_alloc (T);
  unsigned int i ;
  unsigned long int k = 0;
  int status;

  if (seed != 0) {
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


void
rng_state_test (const gsl_rng_type * T)
{
  unsigned long int test_a[1000], test_b[1000] ;

  int i ;

  gsl_rng * r = gsl_rng_alloc (T);
  gsl_rng * r_save = gsl_rng_alloc (T) ;

  for (i = 0; i < 1000; ++i)
    {
      gsl_rng_get (r) ;   /* throw away 1000 iterations */
    }

  gsl_rng_cpy(r_save, r) ;  /* save the intermediate state */

  for (i = 0; i < 1000; ++i)
    {
      test_a[i] = gsl_rng_get (r) ;
    }
    
  gsl_rng_cpy(r, r_save) ;  /* restore the intermediate state */
  gsl_rng_free(r_save) ;

  for (i = 0; i < 1000; ++i)
    {
      test_b[i] = gsl_rng_get (r) ;
    }

  { 
    int status = 0 ;
    for (i = 0; i < 1000; ++i)
      {
	status |= (test_b[i] != test_a[i]) ;
      }
    gsl_test (status, "%s, random number state consistency, 1000 iterations",
	      gsl_rng_name(r)) ;
  }
  
  gsl_rng_free(r) ;
}


void
rng_parallel_state_test (const gsl_rng_type * T)
{
  unsigned long int test_a[1000], test_b[1000] ;

  int i ;

  gsl_rng * r1 = gsl_rng_alloc (T);
  gsl_rng * r2 = gsl_rng_alloc (T) ;

  for (i = 0; i < 1000; ++i)
    {
      gsl_rng_get (r1) ;   /* throw away 1000 iterations */
    }

  gsl_rng_cpy(r2, r1) ;  /* save the intermediate state */

  for (i = 0; i < 1000; ++i)
    {
      test_a[i] = gsl_rng_get (r1) ; /* check that there is no hidden state */
      test_b[i] = gsl_rng_get (r2) ;
    }

  { 
    int status = 0 ;
    for (i = 0; i < 1000; ++i)
      {
	status |= (test_b[i] != test_a[i]) ;
      }
    gsl_test (status, "%s, parallel random number state consistency, "
	      "1000 iterations", gsl_rng_name(r1)) ;
  }

  gsl_rng_free (r1) ;
  gsl_rng_free (r2) ;

}

void
generic_rng_test (const gsl_rng_type * T)
{
  long int n = 1000000;

  gsl_rng * r = gsl_rng_alloc (T);
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
      || (3 * expected_uncovered < actual_uncovered && actual_uncovered > 1);

    gsl_test (status,
	      "%s, observed vs theoretical maximum (%lu vs %lu)",
	      name, kmax, ran_max, n);
  };

  {
    double sum = 0, sigma;
    int i, status;

    for (i = 0; i < n; ++i)
      {
	double x = gsl_rng_get_uni (r) - 0.5;
	sum += x ;
      }
	
    sum /= n;

    /* expect the average to have a variance of 1/(12 n) */

    sigma = sum * sqrt (12.0 * n);
    status = (fabs (sigma) > 3);	/* more than 3 sigma is an error */

    gsl_test (status,
	      "%s, sum test within acceptable sigma (observed %g sigma)",
	      name, sigma);

  }

  gsl_rng_free (r) ;
}

