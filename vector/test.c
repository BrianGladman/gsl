#include <stdio.h>
#include <time.h>
/* #define GSL_RANGE_CHECK */
#include <gsl_vector.h>

#define N 100000
main () 
{
  gsl_vector v ;
  int i,k;
  double *p ;
  double start, end, tot=0 ;

  gsl_vector_alloc(&v, N) ;

  start = clock() / (double)CLOCKS_PER_SEC ;
 
  p = v.data ;

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
    *gsl_vector_set(v,i) = i ;
  } ;
  }

  end = clock() / (double)CLOCKS_PER_SEC ;

  printf("time = %g\n", end - start) ;

  start = clock() / (double)CLOCKS_PER_SEC ;
  
  for (k = 0 ; k < 100 ; k++) {
  for (i = 0 ; i< N ; i++) {
    gsl_vector_set_directly(v,i,i) ;
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

}
