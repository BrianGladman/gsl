#include <config.h>
#include <fcntl.h>
#include <stdio.h>
#include <gsl_vector_complex.h>
#include <gsl_test.h>

#define N 10000

int main (void) 
{
  gsl_vector_complex * v, * w;
  size_t i;

  v = gsl_vector_complex_alloc(N) ;

  gsl_test(v->data == 0, "gsl_vector_complex_alloc returns valid pointer") ;
  gsl_test(v->size != N, "gsl_vector_complex_alloc returns valid size") ;

  for (i = 0 ; i < N ; i++) {
    gsl_complex x ;
    x.real = i ; x.imag = i+1 ;
    gsl_vector_complex_set(v,i,x) ;
  } ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      gsl_complex x, y ;
      x.real = i ; x.imag = i+1 ;
      if(v->data[i].real != x.real || v->data[i].imag != x.imag) 
	status = 1 ;
    } ;
    
    gsl_test(status, "gsl_vector_complex_set writes into array correctly") ;
  }

  {
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      gsl_complex x, y ;
      x.real = i ; x.imag = i+1 ;
      y = gsl_vector_complex_get(v,i) ;
      if(y.real != x.real || y.imag != x.imag) 
	status = 1 ;
    } ;
    gsl_test(status, "gsl_vector_complex_get reads from array correctly") ;
  }

  gsl_vector_complex_free(v);  /* free whatever is in v */

  v = gsl_vector_complex_calloc(N) ;

  gsl_test(v->data == 0, "gsl_vector_complex_calloc returns valid pointer") ;
  gsl_test(v->size != N, "gsl_vector_complex_calloc returns valid size") ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      if(v->data[i].real != 0 || v->data[i].imag != 0) 
	status = 1 ;
    } ;
    
    gsl_test(status, "gsl_vector_complex_calloc initializes array to zero") ;
  }

  {
    FILE * f = fopen("test.txt","w") ;
    
    for (i = 0 ; i < N ; i++) {
      gsl_complex x ;
      x.real = i ; x.imag = i+1 ;
      gsl_vector_complex_set(v,i,x) ;
    } ;
    
    gsl_vector_complex_fprintf(f, v, "%.4f") ;

    fclose(f) ;
  }

  w = gsl_vector_complex_calloc(N) ;

  {
    int status = 0 ;
    FILE * f = fopen("test.txt","r") ;
    
    gsl_vector_complex_fscanf(f, w) ;

    for (i = 0 ; i < N ; i++) {
      if (w->data[i].real != i || w->data[i].imag != i+1 )
	status = 1 ;
    } ;

    gsl_test(status, "gsl_vector_complex_fprintf and fscanf work correctly") ;

    fclose(f) ;
  }


  {
    FILE * f = fopen("test.dat", "w") ;
    
    for (i = 0 ; i < N ; i++) {
      gsl_complex x;
      x.real = (N-i) ; x.imag = (N-i)+1 ;
      gsl_vector_complex_set(v,i,x) ;
    } ;
    
    gsl_vector_complex_fwrite(f, v) ;

    fclose(f) ;
  }

  {
    int status = 0 ;
    FILE * f = fopen("test.dat","r") ;
    
    gsl_vector_complex_fread(f, w) ;

    for (i = 0 ; i < N ; i++) {
      if (w->data[i].real != (N-i) || w->data[i].imag != (N-i)+1 )
	status = 1 ;
    } ;

    gsl_test(status, "gsl_vector_complex_write and read work correctly") ;

    fclose(f) ;
  }

  return gsl_test_summary ();
}
