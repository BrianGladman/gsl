#include <stdio.h>
#include <gsl_rng.h>

int main (void)
{
  gsl_rng * r1 = gsl_rng_alloc (gsl_rng_cmrg) ;
  gsl_rng * r2 = gsl_rng_alloc (gsl_rng_cmrg) ;
  
  int result = gsl_rng_get (r1) ;
  printf("result = %d\n", result) ;

  result = gsl_rng_get (r2) ;
  printf("result = %d\n", result) ;

  result = gsl_rng_get (r1) ;
  printf("result = %d\n", result) ;

  result = gsl_rng_get (r2) ;
  printf("result = %d\n", result) ;

  gsl_rng_set (r1,100) ;
  gsl_rng_set (r1,101) ;

  result = gsl_rng_get (r1) ;
  printf("result = %d\n", result) ;

  result = gsl_rng_get (r2) ;
  printf("result = %d\n", result) ;

  return 0 ;
}
