#include <config.h>
#include <fcntl.h>
#include <stdio.h>
#include <gsl_histogram.h>
#include <gsl_test.h>

#define N 397

void my_error_handler (const char *reason, const char *file, 
		       int line, int err);
int status = 0 ;

int main (void) 
{
  gsl_histogram * h;
  double result ;
  size_t i, j;

  gsl_set_error_handler (&my_error_handler);

  status = 0 ;
  h = gsl_histogram_alloc(0) ;
  gsl_test(!status, "gsl_histogram_alloc traps zero-length histogram") ;
  gsl_test(h != 0, 
	   "gsl_histogram_alloc returns NULL for zero-length histogram") ;

  status = 0 ;
  h = gsl_histogram_alloc_uniform (0, 0.0, 1.0) ;
  gsl_test(!status, 
	   "gsl_histogram_alloc_uniform traps zero-length histogram") ;
  gsl_test(h != 0, 
	   "gsl_histogram_alloc_uniform returns NULL for zero-length histogram") ;

  status = 0 ;
  h = gsl_histogram_alloc_uniform (10, 1.0, 1.0) ;
  gsl_test(!status, 
	   "gsl_histogram_alloc_uniform traps equal endpoints") ;
  gsl_test(h != 0, 
	   "gsl_histogram_alloc_uniform returns NULL for equal endpoints") ;

  status = 0 ;
  h = gsl_histogram_alloc_uniform (10, 2.0, 1.0) ;
  gsl_test(!status, 
	   "gsl_histogram_alloc_uniform traps invalid range") ;
  gsl_test(h != 0, 
	   "gsl_histogram_alloc_uniform returns NULL for invalid range") ;


  h = gsl_histogram_alloc_uniform (N, 0.0, 1.0) ;

  status = gsl_histogram_accumulate (h, 1.0, 10.0) ;
  gsl_test(status != 1, "gsl_histogram_accumulate traps x at xmax") ;

  status = gsl_histogram_accumulate (h, 2.0, 100.0) ;
  gsl_test(status != 1, "gsl_histogram_accumulate traps x above xmax") ;

  status = gsl_histogram_accumulate (h, -1.0, 1000.0) ;
  gsl_test(status != -1, "gsl_histogram_accumulate traps x below xmin") ;

  status = gsl_histogram_increment (h, 1.0) ;
  gsl_test(status != 1, "gsl_histogram_increment traps x at xmax") ;

  status = gsl_histogram_increment (h, 2.0) ;
  gsl_test(status != 1, "gsl_histogram_increment traps x above xmax") ;

  status = gsl_histogram_increment (h, -1.0) ;
  gsl_test(status != -1, "gsl_histogram_increment traps x below xmin") ;


  result = gsl_histogram_get (h, N) ;
  gsl_test(result != 0, "gsl_histogram_get traps index at nbins") ;

  result = gsl_histogram_get (h, N+1) ;
  gsl_test(result != 0, "gsl_histogram_get traps index above nbins") ;

  result = gsl_histogram_get_lowerlimit (h, N) ;
  gsl_test(result != 0, "gsl_histogram_get_lowerlimit traps index at nbins") ;

  result = gsl_histogram_get_upperlimit (h, N) ;
  gsl_test(result != 0, "gsl_histogram_get_upperlimit traps index at nbins") ;

  result = gsl_histogram_get_lowerlimit (h, N + 1) ;
  gsl_test(result != 0, "gsl_histogram_get_lowerlimit traps index above nbins") ;

  result = gsl_histogram_get_upperlimit (h, N + 1) ;
  gsl_test(result != 0, "gsl_histogram_get_upperlimit traps index above nbins") ;

  status = 0;
  status = gsl_histogram_find (h->nbins, h->range, -0.01, &i) ;
  gsl_test(status != -1, "gsl_histogram_find traps x below xmin") ;

  status = 0;
  status = gsl_histogram_find (h->nbins, h->range, 1.0, &i) ;
  gsl_test(status != 1, "gsl_histogram_find traps x at xmax") ;

  status = 0;
  status = gsl_histogram_find (h->nbins, h->range, 1.1, &i) ;
  gsl_test(status != 1, "gsl_histogram_find traps x above xmax") ;

  return gsl_test_summary ();
}


void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
  status = 1 ;
}
