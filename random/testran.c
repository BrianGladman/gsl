
/* $Id$ */
/* Test routine for the uniform random number generators */
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl_test.h>
#include "gsl_ran.h"

#include <getopt.h>

void
usage()
{
  printf("Usage: test [OPTION]\n"
"Exercise a random number routine. By default 1000000 random numbers are\n"
"tested.\n"
"\n"
"  -n, --number=NUM       sample NUM random numbers\n"
"  -s, --seed=NUM         set random number seed (optional)\n"
"  -v, --verbose          verbosely list tests\n"
"\n"
"Without the -v option the test is quiet. The exit status indicates\n"
"success or failure.\n"
) ; 
  exit(0) ;
}

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
  
  while (1) {
    
    static struct option long_options[] = 
    {
      {"verbose", 0, 0, 'v'},
      {"number", 1, 0, 'n'},
      {"seed", 1, 0, 's'},
      {"help", 0, 0, 'h'},
      {0, 0, 0, 0}
    } ;
    
    int option_index = 0 ;
    
    int c = getopt_long (argc, argv, "hn:s:v",
			 long_options, &option_index) ;
    
    if (c == -1)   /* end of options */
      break ;   
    
    if (c == 0 && long_options[option_index].flag == 0)
      c = long_options[option_index].val;
    
    switch (c) 
      {
      case 'v':
	verbose = 1 ;
	break ;

      case 'n':
	if (optarg) 
	  n = strtol (optarg, NULL, 0);
	else 
	  usage () ;
	break ;

      case 's':
	if (optarg) 
	  {
	    randseed = strtol (optarg, NULL, 0);
	    gsl_ran_seed(randseed);
	  }
	else 
	  usage () ;
	break ;

      case 'h':
      default:
	usage () ;
      }
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
  expected_uncovered = ran_max / n ;
  
  status = (rmax > ran_max) || (expected_uncovered < actual_uncovered) ;
  
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
    




