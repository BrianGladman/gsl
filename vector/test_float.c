#include <config.h>
#include <fcntl.h>
#include <stdio.h>
#include <gsl_vector_float.h>
#include <gsl_test.h>

#define N 10000

int main (void) 
{
  gsl_vector_float * v, * w;
  size_t i;

  v = gsl_vector_float_alloc(N) ;

  gsl_test(v->data == 0, "gsl_vector_float_alloc returns valid pointer") ;
  gsl_test(v->size != N, "gsl_vector_float_alloc returns valid size") ;

  for (i = 0 ; i < N ; i++) {
    float x = i ;
    gsl_vector_float_set(v,i,x) ;
  } ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      if(v->data[i] != i) 
	status = 1 ;
    } ;
    
    gsl_test(status, "gsl_vector_float_set writes into array correctly") ;
  }

  {
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      if(gsl_vector_float_get(v,i) != i) 
	status = 1 ;
    } ;
    gsl_test(status, "gsl_vector_float_get reads from array correctly") ;
  }

  gsl_vector_float_free(v);  /* free whatever is in v */

  v = gsl_vector_float_calloc(N) ;

  gsl_test(v->data == 0, "gsl_vector_float_calloc returns valid pointer") ;
  gsl_test(v->size != N, "gsl_vector_float_calloc returns valid size") ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      if(v->data[i] != 0.0) 
	status = 1 ;
    } ;
    
    gsl_test(status, "gsl_vector_float_calloc initializes array to zero") ;
  }

  {
    FILE * f = fopen("test.txt","w") ;
    
    for (i = 0 ; i < N ; i++) {
      gsl_vector_float_set(v,i,(float)i) ;
    } ;
    
    gsl_vector_float_fprintf(f, v, "%.4f") ;

    fclose(f) ;
  }

  w = gsl_vector_float_calloc(N) ;

  {
    int status = 0 ;
    FILE * f = fopen("test.txt","r") ;
    
    gsl_vector_float_fscanf(f, w) ;

    for (i = 0 ; i < N ; i++) {
      if (w->data[i] != (float)i)
	status = 1 ;
    } ;

    gsl_test(status, "gsl_vector_float_fprintf and fscanf work correctly") ;

    fclose(f) ;
  }


  {
    FILE * f = fopen("test.dat", "w") ;
    
    for (i = 0 ; i < N ; i++) {
      gsl_vector_float_set(v,i,(float)(N-i)) ;
    } ;
    
    gsl_vector_float_fwrite(f, v) ;

    fclose(f) ;
  }

  {
    int status = 0 ;
    FILE * f = fopen("test.dat","r") ;
    
    gsl_vector_float_fread(f, w) ;

    for (i = 0 ; i < N ; i++) {
      if (w->data[i] != (float)(N-i))
	status = 1 ;
    } ;

    gsl_test(status, "gsl_vector_float_write and read work correctly") ;

    fclose(f) ;
  }

  return gsl_test_summary ();
}
