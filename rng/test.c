#include <stdio.h>
#include <gsl_rng.h>

int main (void)
{
  gsl_rng * r1 = gsl_rng_alloc (gsl_rng_cmrg) ;
  gsl_rng * r2 = gsl_rng_alloc (gsl_rng_mrg) ;
  
  int result = gsl_rng_get (r1) ;
  printf("result = %d\n", result) ;

  printf("r1 name is %s\n", r1->name) ;
  printf("r1 state is:\n") ; gsl_rng_print_state (r1) ; printf("\n") ;

  result = gsl_rng_get (r2) ;
  printf("result = %d\n", result) ;

  printf("r2 name is %s\n", r2->name) ;

  result = gsl_rng_get (r1) ;
  printf("result = %d\n", result) ;
  printf("r1 state is:\n") ; gsl_rng_print_state (r1) ; printf("\n") ;

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
