#include <config.h>
#include <stdio.h>

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
  gsl_matrix * m;
  size_t i = 0, j = 0;
  double x; 
  
  gsl_warnings_off = 1 ;

  gsl_set_error_handler (&my_error_handler);

  m = gsl_matrix_alloc(N,M) ;

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

  return gsl_test_summary ();
}

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
  status = 1 ;
}
