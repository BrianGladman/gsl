#include <gsl_matrix.h>
#include <gsl_test.h>

#define N 107
#define M 239

int main (void) 
{
  gsl_matrix * m;
  size_t i, j;
  size_t k = 0;

  m = gsl_matrix_alloc(N,M) ;

  gsl_test(m->data == 0, "gsl_matrix_alloc returns valid pointer") ;
  gsl_test(m->size1 != N, "gsl_matrix_alloc returns valid size1") ;
  gsl_test(m->size2 != M, "gsl_matrix_alloc returns valid size2") ;

  for (i = 0 ; i < N ; i++) {
    for (j = 0 ; j < M ; j++) {
      k++ ; gsl_matrix_set(m,i,j,k) ;
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
    
    gsl_test(status, "gsl_matrix_set writes into array correctly") ;
  }

  {
    int status = 0 ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ;
	if(gsl_matrix_get(m,i,j) != k) 
	  status = 1 ;
      } ;
    } ;
    gsl_test(status, "gsl_matrix_get reads from array correctly") ;
  }

  gsl_matrix_free(m);  /* free whatever is in m */

  m = gsl_matrix_calloc(N,M) ;

  gsl_test(m->data == 0, "gsl_matrix_calloc returns valid pointer") ;
  gsl_test(m->size1 != N, "gsl_matrix_calloc returns valid size1") ;
  gsl_test(m->size2 != M, "gsl_matrix_calloc returns valid size2") ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	if(m->data[i*M + j] != 0.0) 
	  status = 1 ;
      } 
    } ;
    gsl_test(status, "gsl_matrix_calloc initializes array to zero") ;
  }

  return gsl_test_summary ();
}
