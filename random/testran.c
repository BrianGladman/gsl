/* Test routine for the uniform random number generators */
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl_test.h>
#include "gsl_ran.h"

/* Usage: test [n] [seed] 
   Exercise a random number routine. By default 1000000 random numbers are
   tested. The exit status indicates success or failure. */

int verbose = 0 ;

int
main(int argc, char *argv[])
{
  int i ;
  int n = 1000000; /* default test samples 1 million random numbers */
  int status ;
  unsigned long r,rmax;
  double ran_max;
  double sum,sigma;
  double actual_uncovered, expected_uncovered ;
  double test_a[5], test_b[5], test_c[5] ;
  int randseed=17;
  void *tmpState;

  if (argc >= 2) 
    n = strtol (argv[1], NULL, 0);
  
  if (argc == 3) 
    {
      randseed = strtol (argv[2], NULL, 0);
      gsl_ran_seed(randseed);
    }

  if (verbose)
    printf("n=%d seed=%d\n", n, randseed) ;

  if (verbose) 
    {
      printf("Random State:\n");
      gsl_ran_printState(gsl_ran_getRandomState());
    }
  
  ran_max = gsl_ran_max() ;

  if (verbose)
    printf("Nominal Maximum random  = %g\n", ran_max);

  rmax=0;
  for (i=0; i<n; ++i) {
    r = gsl_ran_random();
    if (rmax < r) rmax = r;
  }

  if (verbose)
    {
      printf("Largest Observed random = %lu\n",rmax);
      printf("                        = %g\n",(double)rmax);
    }
  
  actual_uncovered = ran_max - rmax ;


  expected_uncovered = ran_max / n ;   /* approximately */

  status = (rmax > ran_max) 
    || (expected_uncovered < actual_uncovered && actual_uncovered > 1) ;

  /* 
     the uni generator never actually reaches its ran_max in practice,
     due to the way the initial state is generated from the seed.
     Thus it only hits 32766 instead of 32767. 
     
     We'll let it pass by checking if the observed max is just 1 below
     the theoretical max.

     */
  

  
  gsl_test (status, 
	    "observed vs theoretical maximum (%lu vs %g)",
	    rmax, ran_max, n) ;
  
  sum=0;
  for (i=0; i<n; ++i)
    sum += gsl_ran_uniform();
  sum /= n;
  /* expect sum to have variance == n*(1/12) */
  /* so average should have variance == 1/(12*n) */
  sigma = (sum - 0.5)*sqrt(12.0*n);

  if (verbose)
    printf("Sum test: %.2f sigmas\n",sigma);

  status = (fabs(sigma) > 3) ;  /* more than 3 sigma is an error */
  gsl_test (status,
	    "sum test within acceptable sigma (observed %g sigma)",
	    sigma) ;

  if (verbose)
    {
      printf("Testing getRandomState/setRandomState:\n");
      printf("The following sets of numbers should be identical.\n");
    }

  tmpState = gsl_ran_getRandomState();
  for (i=0; i<5; ++i) 
    {
      test_a[i] = gsl_ran_uniform() ;
      if (verbose)
	printf("%c %.6f",(i==0 ? '(' : ','), test_a[i]);
    }
  if (verbose)
  printf(" )\n");

  gsl_ran_setRandomState(tmpState); 
  for (i=0; i<5; ++i)
    {
      test_b[i] = gsl_ran_uniform() ;
      if (verbose) 
	printf("%c %.6f",(i==0 ? '(' : ','), test_b[i]);
    }
  if (verbose) 
    printf(" )\n");

  for (i = 0; i < 5; ++i) {
    status = (test_b[i] != test_a[i]) ;
    gsl_test (status, "random number state consistency, iteration %d", i) ;
  }


  gsl_ran_setRandomState(tmpState);
  for (i=0; i<5; ++i) {
    test_c[i] = gsl_ran_uniform_wstate(tmpState) ;
    if (verbose)
      printf("%c %.6f",(i==0 ? '(' : ','), test_c[i]);
    /* we have two parallel streams now, gsl_ran_uniform() */
    gsl_ran_uniform();
    gsl_ran_uniform();
    gsl_ran_uniform();
  }
  if (verbose)
    printf(" )\n");
  
  for (i = 0; i < 5; ++i) {
    status = (test_c[i] != test_b[i]) ;
    gsl_test (status, 
	      "parallel random number state consistency, iteration %d",
	      i) ;
  }

  return gsl_test_summary () ;
}
    




