#include <config.h>
#include <stdio.h>

#undef GSL_RANGE_CHECK_OFF
#define GSL_WARNINGS_OFF

#include <gsl_matrix.h>
#include <gsl_matrix_int.h>
#include <gsl_matrix_float.h>
#include <gsl_matrix_complex.h>
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
  gsl_complex z = {1.2, 3.4} ;
  gsl_complex z1 ;
  gsl_complex zero = {0.0, 0.0} ;
  
  gsl_warnings_off = 1 ;

  gsl_set_error_handler (&my_error_handler);

  /* double */

  status = 0 ;
  gsl_matrix_set(m, i-1, j, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_set traps first index below lower array bound") ;

  status = 0 ;
  gsl_matrix_set(m, i, j-1, 1.2) ;
  gsl_test(!status, 
	   "gsl_matrix_set traps second index below lower array bound") ;

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
  gsl_test(z1.real != zero.real, 
	   "gsl_matrix_complex_get returns zero real for first index below l.b.") ;
  gsl_test(z1.imag != zero.imag, 
	   "gsl_matrix_complex_get returns zero imag for first index below l.b.") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, 0, j-1) ;
  gsl_test(!status,
	   "gsl_matrix_complex_get traps second index below lower bound") ;
  gsl_test(z1.real != zero.real, 
	   "gsl_matrix_complex_get returns zero real for second index below l.b.") ;
  gsl_test(z1.imag != zero.imag, 
	   "gsl_matrix_complex_get returns zero imag for second index below l.b.") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, N+1, 0) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_get traps first index above upper bound") ;
  gsl_test(z1.real != zero.real, 
	   "gsl_matrix_complex_get returns zero real for first index above u.b.") ;
  gsl_test(z1.imag != zero.imag, 
	   "gsl_matrix_complex_get returns zero imag for first index above u.b.") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, 0, M+1) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_get traps second index above upper bound") ;
  gsl_test(z1.real != zero.real, 
	   "gsl_matrix_complex_get returns zero real for second index above u.b.") ;
  gsl_test(z1.imag != zero.imag, 
	   "gsl_matrix_complex_get returns zero imag for second index above u.b.") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, N, 0) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_get traps first index at upper bound") ;
  gsl_test(z1.real != zero.real, 
	   "gsl_matrix_complex_get returns zero real for first index at u.b.") ;
  gsl_test(z1.imag != zero.imag, 
	   "gsl_matrix_complex_get returns zero imag for first index at u.b.") ;

  status = 0 ;
  z1 = gsl_matrix_complex_get(mc, 0, M) ;
  gsl_test(!status, 
	   "gsl_matrix_complex_get traps second index at upper bound") ;
  gsl_test(z1.real != zero.real, 
	   "gsl_matrix_complex_get returns zero real for second index at u.b.") ;
  gsl_test(z1.imag != zero.imag, 
	   "gsl_matrix_complex_get returns zero imag for second index at u.b.") ;

  return gsl_test_summary ();
}

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
  status = 1 ;
}
