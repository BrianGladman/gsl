#include <config.h>
#include <stdio.h>
#include <gsl_matrix_int.h>
#include <gsl_test.h>

#define N 107
#define M 239

int main (void) 
{
  gsl_matrix_int * m, * mm;
  size_t i, j;
  int k = 0;

  m = gsl_matrix_int_alloc(N,M) ;

  gsl_test(m->data == 0, "gsl_matrix_int_alloc returns valid pointer") ;
  gsl_test(m->size1 != N, "gsl_matrix_int_alloc returns valid size1") ;
  gsl_test(m->size2 != M, "gsl_matrix_int_alloc returns valid size2") ;

  for (i = 0 ; i < N ; i++) {
    for (j = 0 ; j < M ; j++) {
      k++ ; gsl_matrix_int_set(m,i,j,k) ;
    }
  }

  { 
    int status = 0 ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ;
	if(m->data[i*M + j] != k) 
	  status = 1 ;
      } ;
    } ;
    
    gsl_test(status, "gsl_matrix_int_set writes into array correctly") ;
  }

  {
    int status = 0 ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ;
	if(gsl_matrix_int_get(m,i,j) != k) 
	  status = 1 ;
      } ;
    } ;
    gsl_test(status, "gsl_matrix_int_get reads from array correctly") ;
  }

  gsl_matrix_int_free(m);  /* free whatever is in m */

  m = gsl_matrix_int_calloc(N,M) ;

  gsl_test(m->data == 0, "gsl_matrix_int_calloc returns valid pointer") ;
  gsl_test(m->size1 != N, "gsl_matrix_int_calloc returns valid size1") ;
  gsl_test(m->size2 != M, "gsl_matrix_int_calloc returns valid size2") ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	if(m->data[i*M + j] != 0.0) 
	  status = 1 ;
      } 
    }
    gsl_test(status, "gsl_matrix_int_calloc initializes array to zero") ;
  }

  {
    FILE * f = fopen("test.txt","w") ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; gsl_matrix_int_set(m,i,j,k) ;
      }
    }
    
    gsl_matrix_int_fprintf(f, m, "%d") ;

    fclose(f) ;
  }

  mm = gsl_matrix_int_calloc(N,M) ;

  {
    int status = 0 ;
    FILE * f = fopen("test.txt","r") ;
    
    gsl_matrix_int_fscanf(f, mm) ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; if (mm->data[i*M + j] != k)
	status = 1 ;
      }
    }

    gsl_test(status, "gsl_matrix_int_fprintf and fscanf work correctly") ;

    fclose(f) ;
  }

  {
    FILE * f = fopen("test.dat", "w") ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; gsl_matrix_int_set(m,i,j,k) ;
      }
    }
    
    gsl_matrix_int_fwrite(f, m) ;

    fclose(f) ;
  }

  {
    int status = 0 ;
    FILE * f = fopen("test.dat","r") ;
    
    gsl_matrix_int_fread(f, m) ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; if (m->data[i*M + j] != k)
	  status = 1 ;
      }
    }

    gsl_test(status, "gsl_matrix_int_write and read work correctly") ;

    fclose(f) ;
  }

  return gsl_test_summary ();
}
