#include <config.h>
#include <fcntl.h>
#include <stdio.h>
#include <gsl_vector.h>
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
    gsl_complex x = {{i,i+1}};
    gsl_vector_complex_set(v,i,x) ;
  } ;
  
  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      double xr = i, xi = i+1;
      if (v->data[2*i] != xr || v->data[2*i+1] != xi)
	status = 1 ;
    } ;

    gsl_test(status, "gsl_vector_complex_set writes into array correctly") ;
  }

  {
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      gsl_complex x = {{i,i+1}} ;
      gsl_complex y = GSL_VECTOR_COMPLEX(v,i) ;
      if(!GSL_COMPLEX_EQ(x,y)) 
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
      if(v->data[2*i] != 0 || v->data[2*i+1] != 0) 
	status = 1 ;
    } ;
    
    gsl_test(status, "gsl_vector_complex_calloc initializes array to zero") ;
  }

  {
    FILE * f = fopen("test.txt","w") ;
    
    for (i = 0 ; i < N ; i++) {
      gsl_complex x = {{i,i+1}} ;
      GSL_VECTOR_COMPLEX(v,i) = x ;
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
      if (w->data[2*i] != i || w->data[2*i+1] != i+1 )
	status = 1 ;
    } ;

    gsl_test(status, "gsl_vector_complex_fprintf and fscanf work correctly") ;

    fclose(f) ;
  }


  {
    FILE * f = fopen("test.dat", "w") ;
    
    for (i = 0 ; i < N ; i++) {
      gsl_complex x = {{(N-i),(N-i)+1}} ;
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
      if (w->data[2*i] != (N-i) || w->data[2*i+1] != (N-i)+1 )
	status = 1 ;
    } ;

    gsl_test(status, "gsl_vector_complex_write and read work correctly") ;

    fclose(f) ;
  }

  return gsl_test_summary ();
}
