#include <config.h>
#include <math.h>
#include <gsl_histogram.h>
#include <gsl_test.h>

int main (void)
{
  size_t i;
  int status = 0 ;

  gsl_histogram * h = gsl_histogram_alloc_uniform (10, 0.0, 1.0) ;
  
  gsl_histogram_increment (h, 0.1) ;
  gsl_histogram_increment (h, 0.2) ;
  gsl_histogram_increment (h, 0.2) ;
  gsl_histogram_increment (h, 0.3) ;
  
  { 
    gsl_histogram_pdf * p = gsl_histogram_pdf_alloc (h) ;

    gsl_histogram * hh = gsl_histogram_alloc_uniform (100, 0.0, 1.0) ;

    for (i = 0; i < 100000 ; i++)
      {
	double u = ((double) rand ()) / RAND_MAX;
	double x = gsl_histogram_pdf_sample (p, u) ;
	gsl_histogram_increment (hh, x) ;
      }

    for (i = 0 ; i < 100 ; i++)
      {
	double y = gsl_histogram_get (hh, i) /  2500  ;
	double x = gsl_histogram_get_lowerlimit (hh, i) ;
	size_t k ; double ya ;

	gsl_histogram_find (h, x, &k) ;
	ya = gsl_histogram_get (h, k) ;

	if (ya == 0) {
	  if (y != 0) 
	    status = 1 ;
	} else {
	  double err = 1/sqrt(gsl_histogram_get (hh, i)) ;
	  double sigma = fabs((y-ya)/(y*err)) ;
	  if (sigma > 3) 
	    status = 1 ;
	}
      }
    gsl_test(status, "gsl_histogram_pdf_sample works") ;
  }
  
  return gsl_test_summary () ;
}

