#include <config.h>
#include <fcntl.h>
#include <stdio.h>
#include <gsl_vector.h>
#include <gsl_test.h>

#define N 10000

int main (void) 
{
  gsl_vector * v, * w;
  size_t i;

  v = gsl_vector_alloc(N) ;

  gsl_test(v->data == 0, "gsl_vector_alloc returns valid pointer") ;
  gsl_test(v->size != N, "gsl_vector_alloc returns valid size") ;

  for (i = 0 ; i < N ; i++) {
    gsl_vector_set(v,i,(double)i) ;
  } ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      if(v->data[i] != i) 
	status = 1 ;
    } ;
    
    gsl_test(status, "gsl_vector_set writes into array correctly") ;
  }

  {
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      if(gsl_vector_get(v,i) != i) 
	status = 1 ;
    } ;
    gsl_test(status, "gsl_vector_get reads from array correctly") ;
  }

  gsl_vector_free(v);  /* free whatever is in v */

  v = gsl_vector_calloc(N) ;

  gsl_test(v->data == 0, "gsl_vector_calloc returns valid pointer") ;
  gsl_test(v->size != N, "gsl_vector_calloc returns valid size") ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      if(v->data[i] != 0.0) 
	status = 1 ;
    } ;
    
    gsl_test(status, "gsl_vector_calloc initializes array to zero") ;
  }

  {
    FILE * f = fopen("test.dat","w") ;
    
    for (i = 0 ; i < N ; i++) {
      gsl_vector_set(v,i,(double)i) ;
    } ;
    
    gsl_vector_fprintf(f, v, "%.4f") ;

    fclose(f) ;
  }

  w = gsl_vector_calloc(N) ;

  {
    int status = 0 ;
    FILE * f = fopen("test.dat","r") ;
    
    gsl_vector_fscanf(f, w) ;

    for (i = 0 ; i < N ; i++) {
      if (w->data[i] != (double)i)
	status = 1 ;
    } ;

    gsl_test(status, "gsl_vector_fprintf and fscanf work correctly") ;

    fclose(f) ;
  }


  {
    int f = open("test.dat",O_WRONLY|O_TRUNC) ;
    
    for (i = 0 ; i < N ; i++) {
      gsl_vector_set(v,i,(double)(N-i)) ;
    } ;
    
    gsl_vector_write(f, v) ;

    close(f) ;
  }

  {
    int status = 0 ;
    int  f = open("test.dat",O_RDONLY) ;
    
    gsl_vector_read(f, w) ;

    for (i = 0 ; i < N ; i++) {
      if (w->data[i] != (double)(N-i))
	status = 1 ;
    } ;

    gsl_test(status, "gsl_vector_write and read work correctly") ;

    close(f) ;
  }

  return gsl_test_summary ();
}
