#include <stdio.h>
#include <time.h>
/* #define GSL_RANGE_CHECK */
#include <gsl_vector.h>

#define N 100000

int main (void) 
{
  gsl_vector * v ;
  size_t i,k;
  double *p ;
  double start, end, tot=0 ;

  v = gsl_vector_alloc(N) ;

  start = clock() / (double)CLOCKS_PER_SEC ;
 
  p = v->data ;

  for (k = 0 ; k < 100 ; k++) {
  for (i = 0 ; i< N ; i++) {
    p[i] = i ;
  } ;
  }

  end = clock() / (double)CLOCKS_PER_SEC ;

  printf("time = %g\n", end - start) ;
  
  start = clock() / (double)CLOCKS_PER_SEC ;
 
  for (k = 0 ; k < 100 ; k++) {
  for (i = 0 ; i< N ; i++) {
    gsl_vector_set(v,i,(double)i) ;
  } ;
  }

  end = clock() / (double)CLOCKS_PER_SEC ;

  printf("time = %g\n", end - start) ;

  start = clock() / (double)CLOCKS_PER_SEC ;

  for (k = 0 ; k < 100 ; k++) {
  for (i = 0 ; i< N ; i++) {
    tot += gsl_vector_get(v,i) ;
  } ;
  }
  end = clock() / (double)CLOCKS_PER_SEC ;

  printf("tot=%g\n",tot) ;
  printf("time = %g\n", end - start) ;
  return 0 ;
}
