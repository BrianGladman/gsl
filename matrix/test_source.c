void test (void);
 
void 
test (void) 
{
  BASE * m;
  size_t i, j;
  int k = 0;

  m = FUNCTION(gsl_matrix,alloc) (N,M) ;

  gsl_test(m->data == 0, NAME(gsl_matrix) "_alloc returns valid pointer") ;
  gsl_test(m->size1 != N, NAME(gsl_matrix) "_alloc returns valid size1") ;
  gsl_test(m->size2 != M, NAME(gsl_matrix) "_alloc returns valid size2") ;

  for (i = 0 ; i < N ; i++) {
    for (j = 0 ; j < M ; j++) {
      k++ ; FUNCTION(gsl_matrix,set)(m,i,j,(BASE)k) ;
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
    
    gsl_test(status, NAME(gsl_matrix) "_set writes into array correctly") ;
  }

  {
    int status = 0 ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ;
	if(FUNCTION(gsl_matrix,get)(m,i,j) != k) 
	  status = 1 ;
      } ;
    } ;
    gsl_test(status, NAME(gsl_matrix) "_get reads from array correctly") ;
  }

  FUNCTION(gsl_matrix,free)(m);  /* free whatever is in m */
  
  m = FUNCTION(gsl_matrix,calloc)(N,M) ;

  gsl_test(m->data == 0, NAME(gsl_matrix) "_calloc returns valid pointer") ;
  gsl_test(m->size1 != N, NAME(gsl_matrix) "_calloc returns valid size1") ;
  gsl_test(m->size2 != M, NAME(gsl_matrix) "_calloc returns valid size2") ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	if(m->data[i*M + j] != 0.0) 
	  status = 1 ;
      } 
    }
    gsl_test(status, NAME(gsl_matrix) "_calloc initializes array to zero") ;
  }

  {
    FILE * f = fopen("test.txt","w") ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; FUNCTION(gsl_matrix,set)(m,i,j,(double)k) ;
      }
    }
    
    FUNCTION(gsl_matrix,fprintf)(f, m, "%.4f") ;

    fclose(f) ;
  }


  {
    int status = 0 ;
    FILE * f = fopen("test.txt","r") ;
    BASE * mm = FUNCTION(gsl_matrix,calloc)(N,M) ;

    FUNCTION(gsl_matrix,fscanf)(f, mm) ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; 
	if (mm->data[i*M + j] != k)
	  status = 1 ;
      }
    }
    
    gsl_test(status, NAME(gsl_matrix) "_fprintf and fscanf work correctly") ;
    
    fclose (f) ;
    gsl_matrix_free (mm) ;
  }


  {
    FILE * f = fopen("test.dat", "w") ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; FUNCTION(gsl_matrix,set)(m,i,j,(double)k) ;
      }
    }

    FUNCTION(gsl_matrix,fwrite)(f, m) ;
    fclose (f) ;
  }

  {
    int status = 0 ;
    FILE * f = fopen("test.dat","r") ;
    BASE * mm = FUNCTION(gsl_matrix,calloc)(N,M) ;

    FUNCTION(gsl_matrix,fread)(f, mm) ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; 
	if (mm->data[i*M + j] != k)
	    status = 1 ;
      }
    }

    gsl_test(status, NAME(gsl_matrix) "_write and read work correctly") ;

    fclose (f) ;
    FUNCTION(gsl_matrix,free) (mm) ;
  }

}

#undef GSL_RANGE_CHECK_OFF
#define GSL_WARNINGS_OFF

#include <gsl_matrix.h>
#include <gsl_test.h>

#define N 107
#define M 239
void my_error_handler (const char *reason, const char *file, 
		       int line, int err);

int status = 0 ;

int main (void) 
{
  gsl_matrix * m = gsl_matrix_alloc(N,M) ;
  gsl_matrix_int * mi = gsl_matrix_int_alloc(N,M) ;
  gsl_matrix_float * mf = gsl_matrix_float_alloc(N,M) ;
  gsl_matrix_complex * mc = gsl_matrix_complex_alloc(N,M) ;

  size_t i = 0, j = 0;
  double x; 
  gsl_complex z = {{1.2, 3.4}} ;
  gsl_complex z1 ;
  gsl_complex zero = {{0.0, 0.0}} ;
  
  gsl_warnings_off = 1 ;

  gsl_set_error_handler (&my_error_handler);

  /* double */

  status = 0 ;
  gsl_matrix_set(m, N+1, 0, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_set traps first index above upper array bound") ;

  status = 0 ;
  gsl_matrix_set(m, 0, M+1, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_set traps second index above upper array bound") ;

  status = 0 ;
  gsl_matrix_set(m, N, 0, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_set traps first index at upper array bound") ;

  status = 0 ;
  gsl_matrix_set(m, 0, M, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_set traps second index at upper array bound") ;

  status = 0 ;
  x = gsl_matrix_get(m, i-1, 0) ;
  gsl_test(!status,
	   "gsl_matrix_get traps first index below lower array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_get returns zero for first index below lower array bound") ;

  status = 0 ;
  x = gsl_matrix_get(m, 0, j-1) ;
  gsl_test(!status,
	   "gsl_matrix_get traps second index below lower array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_get returns zero for second index below lower array bound") ;

  status = 0 ;
  x = gsl_matrix_get(m, N+1, 0) ;
  gsl_test(!status, 
	   "gsl_matrix_get traps first index above upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_get returns zero for first index above upper array bound") ;

  status = 0 ;
  x = gsl_matrix_get(m, 0, M+1) ;
  gsl_test(!status, 
	   "gsl_matrix_get traps second index above upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_get returns zero for second index above upper array bound") ;

  status = 0 ;
  x = gsl_matrix_get(m, N, 0) ;
  gsl_test(!status, 
	   "gsl_matrix_get traps first index at upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_get returns zero for first index at upper array bound") ;

  status = 0 ;
  x = gsl_matrix_get(m, 0, M) ;
  gsl_test(!status, 
	   "gsl_matrix_get traps second index at upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_get returns zero for second index at upper array bound") ;

  /* float */

  status = 0 ;
  gsl_matrix_float_set(mf, i-1, j, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_float_set traps first index below lower array bound") ;

  status = 0 ;
  gsl_matrix_float_set(mf, i, j-1, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_float_set traps second index below lower array bound") ;

  status = 0 ;
  gsl_matrix_float_set(mf, N+1, 0, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_float_set traps first index above upper array bound") ;

  status = 0 ;
  gsl_matrix_float_set(mf, 0, M+1, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_float_set traps second index above upper array bound") ;

  status = 0 ;
  gsl_matrix_float_set(mf, N, 0, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_float_set traps first index at upper array bound") ;

  status = 0 ;
  gsl_matrix_float_set(mf, 0, M, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_float_set traps second index at upper array bound") ;

  status = 0 ;
  x = gsl_matrix_float_get(mf, i-1, 0) ;
  gsl_test(!status,
	   "gsl_matrix_float_get traps first index below lower array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_float_get returns zero for first index below lower array bound") ;

  status = 0 ;
  x = gsl_matrix_float_get(mf, 0, j-1) ;
  gsl_test(!status,
	   "gsl_matrix_float_get traps second index below lower array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_float_get returns zero for second index below lower array bound") ;

  status = 0 ;
  x = gsl_matrix_float_get(mf, N+1, 0) ;
  gsl_test(!status, 
	   "gsl_matrix_float_get traps first index above upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_float_get returns zero for first index above upper array bound") ;

  status = 0 ;
  x = gsl_matrix_float_get(mf, 0, M+1) ;
  gsl_test(!status, 
	   "gsl_matrix_float_get traps second index above upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_float_get returns zero for second index above upper array bound") ;

  status = 0 ;
  x = gsl_matrix_float_get(mf, N, 0) ;
  gsl_test(!status, 
	   "gsl_matrix_float_get traps first index at upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_float_get returns zero for first index at upper array bound") ;

  status = 0 ;
  x = gsl_matrix_float_get(mf, 0, M) ;
  gsl_test(!status, 
	   "gsl_matrix_float_get traps second index at upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_float_get returns zero for second index at upper array bound") ;



  /* int */

  status = 0 ;
  gsl_matrix_int_set(mi, i-1, j, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_int_set traps first index below lower array bound") ;

  status = 0 ;
  gsl_matrix_int_set(mi, i, j-1, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_int_set traps second index below lower array bound") ;

  status = 0 ;
  gsl_matrix_int_set(mi, N+1, 0, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_int_set traps first index above upper array bound") ;

  status = 0 ;
  gsl_matrix_int_set(mi, 0, M+1, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_int_set traps second index above upper array bound") ;

  status = 0 ;
  gsl_matrix_int_set(mi, N, 0, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_int_set traps first index at upper array bound") ;

  status = 0 ;
  gsl_matrix_int_set(mi, 0, M, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_int_set traps second index at upper array bound") ;

  status = 0 ;
  x = gsl_matrix_int_get(mi, i-1, 0) ;
  gsl_test(!status,
	   "gsl_matrix_int_get traps first index below lower array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_int_get returns zero for first index below lower array bound") ;

  status = 0 ;
  x = gsl_matrix_int_get(mi, 0, j-1) ;
  gsl_test(!status,
	   "gsl_matrix_int_get traps second index below lower array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_int_get returns zero for second index below lower array bound") ;

  status = 0 ;
  x = gsl_matrix_int_get(mi, N+1, 0) ;
  gsl_test(!status, 
	   "gsl_matrix_int_get traps first index above upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_int_get returns zero for first index above upper array bound") ;

  status = 0 ;
  x = gsl_matrix_int_get(mi, 0, M+1) ;
  gsl_test(!status, 
	   "gsl_matrix_int_get traps second index above upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_int_get returns zero for second index above upper array bound") ;

  status = 0 ;
  x = gsl_matrix_int_get(mi, N, 0) ;
  gsl_test(!status, 
	   "gsl_matrix_int_get traps first index at upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_int_get returns zero for first index at upper array bound") ;

  status = 0 ;
  x = gsl_matrix_int_get(mi, 0, M) ;
  gsl_test(!status, 
	   "gsl_matrix_int_get traps second index at upper array bound") ;
  gsl_test(x != 0, 
	   "gsl_matrix_int_get returns zero for second index at upper array bound") ;


  /* complex */

  status = 0 ;
  gsl_matrix_complex_set(mc, i-1, j, z) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_set traps first index below lower array bound") ;

  status = 0 ;
  gsl_matrix_complex_set(mc, i, j-1, z) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_set traps second index below lower array bound") ;

  status = 0 ;
  gsl_matrix_complex_set(mc, N+1, 0, z) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_set traps first index above upper array bound") ;

  status = 0 ;
  gsl_matrix_complex_set(mc, 0, M+1, z) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_set traps second index above upper array bound") ;

  status = 0 ;
  gsl_matrix_complex_set(mc, N, 0, z) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_set traps first index at upper array bound") ;

  status = 0 ;
  gsl_matrix_complex_set(mc, 0, M, z) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_set traps second index at upper array bound") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, i-1, 0) ;
  gsl_test(!status,
	   "gsl_matrix_complex_get traps first index below lower bound") ;
  gsl_test(z1.dat[0] != zero.dat[0], 
	   "gsl_matrix_complex_get returns zero real for first index below l.b.") ;
  gsl_test(z1.dat[1] != zero.dat[1], 
	   "gsl_matrix_complex_get returns zero imag for first index below l.b.") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, 0, j-1) ;
  gsl_test(!status,
	   "gsl_matrix_complex_get traps second index below lower bound") ;
  gsl_test(z1.dat[0] != zero.dat[0], 
	   "gsl_matrix_complex_get returns zero real for second index below l.b.") ;
  gsl_test(z1.dat[1] != zero.dat[1], 
	   "gsl_matrix_complex_get returns zero imag for second index below l.b.") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, N+1, 0) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_get traps first index above upper bound") ;
  gsl_test(z1.dat[0] != zero.dat[0], 
	   "gsl_matrix_complex_get returns zero real for first index above u.b.") ;
  gsl_test(z1.dat[1] != zero.dat[1], 
	   "gsl_matrix_complex_get returns zero imag for first index above u.b.") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, 0, M+1) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_get traps second index above upper bound") ;
  gsl_test(z1.dat[0] != zero.dat[0], 
	   "gsl_matrix_complex_get returns zero real for second index above u.b.") ;
  gsl_test(z1.dat[1] != zero.dat[1], 
	   "gsl_matrix_complex_get returns zero imag for second index above u.b.") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, N, 0) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_get traps first index at upper bound") ;
  gsl_test(z1.dat[0] != zero.dat[0], 
	   "gsl_matrix_complex_get returns zero real for first index at u.b.") ;
  gsl_test(z1.dat[1] != zero.dat[1], 
	   "gsl_matrix_complex_get returns zero imag for first index at u.b.") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, 0, M) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_get traps second index at upper bound") ;
  gsl_test(z1.dat[0] != zero.dat[0], 
	   "gsl_matrix_complex_get returns zero real for second index at u.b.") ;
  gsl_test(z1.dat[1] != zero.dat[1], 
	   "gsl_matrix_complex_get returns zero imag for second index at u.b.") ;

  return gsl_test_summary ();
}

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
  status = 1 ;
}
