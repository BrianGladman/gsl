#include <stdio.h>
#include <gsl_rng.h>

int main (void)
{
  gsl_rng * r = gsl_rng_alloc (gsl_rng_blah) ;
  
  int result = gsl_rng_get (r) ;
  printf("result = %d\n", result) ;

  result = gsl_rng_get (r) ;
  printf("result = %d\n", result) ;

  gsl_rng_set (r,100) ;

  result = gsl_rng_get (r) ;
  printf("result = %d\n", result) ;

  return 0 ;
}
