#include <config.h>
#include <fcntl.h>
#include <stdio.h>
#include <gsl_histogram2d.h>
#include <gsl_test.h>
#include <gsl_errno.h>

#define N 107
#define M 239

void my_error_handler (const char *reason, const char *file, 
		       int line, int err);
int status = 0 ;

int main (void) 
{
  gsl_histogram2d * h;
  double result ;
  size_t i, j;

  gsl_set_error_handler (&my_error_handler);

  status = 0 ;
  h = gsl_histogram2d_calloc(0,10) ;
  gsl_test(!status, "gsl_histogram_calloc traps zero-width histogram") ;
  gsl_test(h != 0, 
	   "gsl_histogram2d_calloc returns NULL for zero-width histogram") ;


  status = 0 ;
  h = gsl_histogram2d_calloc(10,0) ;
  gsl_test(!status, "gsl_histogram_calloc traps zero-length histogram") ;
  gsl_test(h != 0, 
	   "gsl_histogram2d_calloc returns NULL for zero-length histogram") ;


  status = 0 ;
  h = gsl_histogram2d_calloc_uniform (0, 10, 0.0, 1.0, 0.0, 1.0) ;
  gsl_test(!status, 
	   "gsl_histogram2d_calloc_uniform traps zero-width histogram") ;
  gsl_test(h != 0, 
	   "gsl_histogram2d_calloc_uniform returns NULL for zero-width histogram") ;

  status = 0 ;
  h = gsl_histogram2d_calloc_uniform (10, 0, 0.0, 1.0, 0.0, 1.0) ;
  gsl_test(!status, 
	   "gsl_histogram2d_calloc_uniform traps zero-length histogram") ;
  gsl_test(h != 0, 
	   "gsl_histogram2d_calloc_uniform returns NULL for zero-length histogram") ;

  status = 0 ;
  h = gsl_histogram2d_calloc_uniform (10, 10, 0.0, 1.0, 1.0, 1.0) ;
  gsl_test(!status, 
	   "gsl_histogram2d_calloc_uniform traps equal endpoints") ;
  gsl_test(h != 0, 
	   "gsl_histogram2d_calloc_uniform returns NULL for equal endpoints") ;


  status = 0 ;
  h = gsl_histogram2d_calloc_uniform (10, 10, 1.0, 1.0, 0.0, 1.0) ;
  gsl_test(!status, 
	   "gsl_histogram2d_calloc_uniform traps equal endpoints") ;
  gsl_test(h != 0, 
	   "gsl_histogram2d_calloc_uniform returns NULL for equal endpoints") ;


  status = 0 ;
  h = gsl_histogram2d_calloc_uniform (10, 10, 0.0, 1.0, 2.0, 1.0) ;
  gsl_test(!status, 
	   "gsl_histogram2d_calloc_uniform traps invalid range") ;
  gsl_test(h != 0, 
	   "gsl_histogram2d_calloc_uniform returns NULL for invalid range") ;


  status = 0 ;
  h = gsl_histogram2d_calloc_uniform (10, 10, 2.0, 1.0, 0.0, 1.0) ;
  gsl_test(!status, 
	   "gsl_histogram2d_calloc_uniform traps invalid range") ;
  gsl_test(h != 0, 
	   "gsl_histogram2d_calloc_uniform returns NULL for invalid range") ;


  h = gsl_histogram2d_calloc_uniform (N, M, 0.0, 1.0, 0.0, 1.0) ;

  status = gsl_histogram2d_accumulate (h, 1.0, 0.0, 10.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_accumulate traps x at xmax") ;

  status = gsl_histogram2d_accumulate (h, 2.0, 0.0, 100.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_accumulate traps x above xmax") ;

  status = gsl_histogram2d_accumulate (h, -1.0, 0.0, 1000.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_accumulate traps x below xmin") ;

  status = gsl_histogram2d_accumulate (h, 0.0, 1.0, 10.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_accumulate traps y at ymax") ;

  status = gsl_histogram2d_accumulate (h, 0.0, 2.0, 100.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_accumulate traps y above ymax") ;

  status = gsl_histogram2d_accumulate (h, 0.0, -1.0, 1000.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_accumulate traps y below ymin") ;


  status = gsl_histogram2d_increment (h, 1.0, 0.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_increment traps x at xmax") ;

  status = gsl_histogram2d_increment (h, 2.0, 0.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_increment traps x above xmax") ;

  status = gsl_histogram2d_increment (h, -1.0, 0.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_increment traps x below xmin") ;


  status = gsl_histogram2d_increment (h, 0.0, 1.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_increment traps y at ymax") ;

  status = gsl_histogram2d_increment (h, 0.0, 2.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_increment traps y above ymax") ;

  status = gsl_histogram2d_increment (h, 0.0, -1.0) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_increment traps y below ymin") ;


  result = gsl_histogram2d_get (h, N, 0) ;
  gsl_test(result != 0, "gsl_histogram2d_get traps x index at nx") ;

  result = gsl_histogram2d_get (h, N+1, 0) ;
  gsl_test(result != 0, "gsl_histogram2d_get traps x index above nx") ;

  result = gsl_histogram2d_get (h, 0, M) ;
  gsl_test(result != 0, "gsl_histogram2d_get traps y index at ny") ;

  result = gsl_histogram2d_get (h, 0, M+1) ;
  gsl_test(result != 0, "gsl_histogram2d_get traps y index above ny") ;


  result = gsl_histogram2d_get_xlowerlimit (h, N) ;
  gsl_test(result != 0, "gsl_histogram2d_get_xlowerlimit traps index at nx") ;

  result = gsl_histogram2d_get_xupperlimit (h, N) ;
  gsl_test(result != 0, "gsl_histogram2d_get_xupperlimit traps index at nx") ;

  result = gsl_histogram2d_get_xlowerlimit (h, N + 1) ;
  gsl_test(result != 0, "gsl_histogram2d_get_xlowerlimit traps index above nx") ;

  result = gsl_histogram2d_get_xupperlimit (h, N + 1) ;
  gsl_test(result != 0, "gsl_histogram2d_get_xupperlimit traps index above nx") ;


  result = gsl_histogram2d_get_ylowerlimit (h, M) ;
  gsl_test(result != 0, "gsl_histogram2d_get_ylowerlimit traps index at ny") ;

  result = gsl_histogram2d_get_yupperlimit (h, M) ;
  gsl_test(result != 0, "gsl_histogram2d_get_yupperlimit traps index at ny") ;

  result = gsl_histogram2d_get_ylowerlimit (h, M + 1) ;
  gsl_test(result != 0, "gsl_histogram2d_get_ylowerlimit traps index above ny") ;

  result = gsl_histogram2d_get_yupperlimit (h, M + 1) ;
  gsl_test(result != 0, "gsl_histogram2d_get_yupperlimit traps index above ny") ;

  status = 0;
  gsl_histogram2d_find (h, -0.01, 0.0, &i, &j) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_find traps x below xmin") ;

  status = 0;
  gsl_histogram2d_find (h, 1.0, 0.0, &i, &j) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_find traps x at xmax") ;

  status = 0;
  gsl_histogram2d_find (h, 1.1, 0.0, &i, &j) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_find traps x above xmax") ;


  status = 0;
  gsl_histogram2d_find (h, 0.0, -0.01, &i, &j) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_find traps y below ymin") ;

  status = 0;
  gsl_histogram2d_find (h, 0.0, 1.0, &i, &j) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_find traps y at ymax") ;

  status = 0;
  gsl_histogram2d_find (h, 0.0, 1.1, &i, &j) ;
  gsl_test(status != GSL_EDOM, "gsl_histogram2d_find traps y above ymax") ;

  return gsl_test_summary ();
}


void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
  status = 1 ;
}
