
int main (void)
{
  gsl_rng * r = gsl_rng_alloc() ;
  
  int result = gsl_rng_get(r) ;
  printf("result = %d\n", result) ;

  result = gsl_rng_get(r) ;
  printf("result = %d\n", result) ;

  gsl_rng_set(r,100.0) ;

  result = gsl_rng_get(r) ;
  printf("result = %d\n", result) ;

  return 0 ;
}
