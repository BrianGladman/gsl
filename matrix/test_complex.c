#include <config.h>
#include <stdio.h>
#include <gsl_matrix_complex.h>
#include <gsl_test.h>

#define N 107
#define M 239

int main (void) 
{
  gsl_matrix_complex * m;
  size_t i, j;
  int k = 0;

  m = gsl_matrix_complex_alloc(N,M) ;

  gsl_test(m->data == 0, "gsl_matrix_complex_alloc returns valid pointer") ;
  gsl_test(m->size1 != N, "gsl_matrix_complex_alloc returns valid size1") ;
  gsl_test(m->size2 != M, "gsl_matrix_complex_alloc returns valid size2") ;

  for (i = 0 ; i < N ; i++) {
    for (j = 0 ; j < M ; j++) {
      gsl_complex z ;
      k++ ; z.real = k ; z.imag = k + 1000 ;
      gsl_matrix_complex_set(m,i,j,z) ;
    }
  }

  { 
    int status = 0 ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ;
	if(m->data[i*M + j].real != k || m->data[i*M + j].imag != k + 1000) 
	  status = 1 ;
      } ;
    } ;
    
    gsl_test(status, "gsl_matrix_complex_set writes into array correctly") ;
  }

  {
    int status = 0 ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	gsl_complex z = gsl_matrix_complex_get(m,i,j) ;
	k++ ;
	if(z.real != k || z.imag != k + 1000) 
	  status = 1 ;
      } ;
    } ;
    gsl_test(status, "gsl_matrix_complex_get reads from array correctly") ;
  }

  gsl_matrix_complex_free(m);  /* free whatever is in m */

  m = gsl_matrix_complex_calloc(N,M) ;

  gsl_test(m->data == 0, "gsl_matrix_complex_calloc returns valid pointer") ;
  gsl_test(m->size1 != N, "gsl_matrix_complex_calloc returns valid size1") ;
  gsl_test(m->size2 != M, "gsl_matrix_complex_calloc returns valid size2") ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	if(m->data[i*M + j].real != 0.0 || m->data[i*M + j].imag != 0.0) 
	  status = 1 ;
      } 
    }
    gsl_test(status, "gsl_matrix_complex_calloc initializes array to zero") ;
  }

  {
    FILE * f = fopen("test.txt","w") ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	gsl_complex z ;
	k++ ; z.real = k ; z.imag = k + 1000 ;
	gsl_matrix_complex_set(m,i,j,z) ;
      }
    }
    
    gsl_matrix_complex_fprintf(f, m, "%.4f") ;

    fclose(f) ;
  }

  {
    int status = 0 ;
    FILE * f = fopen("test.txt","r") ;
    gsl_matrix_complex * mm = gsl_matrix_complex_calloc(N,M) ;

    gsl_matrix_complex_fscanf(f, mm) ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; 
	if (mm->data[i*M + j].real != k || mm->data[i*M + j].imag != k + 1000)
	  status = 1 ;
      }
    }
    
    gsl_test(status, "gsl_matrix_complex_fprintf and fscanf work correctly") ;
    
    fclose (f) ;
    gsl_matrix_complex_free (mm) ;
  }

  {
    FILE * f = fopen("test.dat", "w") ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	gsl_complex z ;
	k++ ; z.real = k ; z.imag = k + 1000 ;
	gsl_matrix_complex_set(m,i,j,z) ;
      }
    }
    
    gsl_matrix_complex_fwrite(f, m) ;

    fclose (f) ;
  }

  {
    int status = 0 ;
    FILE * f = fopen("test.dat","r") ;
    gsl_matrix_complex * mm = gsl_matrix_complex_calloc(N,M) ;

    gsl_matrix_complex_fread(f, mm) ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; 
	if (mm->data[i*M + j].real != k || mm->data[i*M + j].imag != k + 1000 )
	    status = 1 ;
      }
    }
    
    gsl_test(status, "gsl_matrix_complex_write and read work correctly") ;

    fclose (f) ;
    gsl_matrix_complex_free (mm) ;
  }

  return gsl_test_summary ();
}
