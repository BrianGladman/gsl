#include <config.h>
#include <math.h>
#include <gsl_histogram2d.h>
#include <gsl_test.h>

int main (void)
{
  size_t i, j;
  int status = 0 ;

  gsl_histogram2d * h = gsl_histogram2d_alloc_uniform (10, 10,
						       0.0, 1.0, 
						       0.0, 1.0) ;
  
  gsl_histogram2d_increment (h, 0.1, 0.1) ;
  gsl_histogram2d_increment (h, 0.2, 0.2) ;
  gsl_histogram2d_increment (h, 0.2, 0.2) ;
  gsl_histogram2d_increment (h, 0.3, 0.3) ;
  
  { 
    gsl_histogram2d_pdf * p = gsl_histogram2d_pdf_alloc (h) ;

    gsl_histogram2d * hh = gsl_histogram2d_alloc_uniform (20, 20, 
							  0.0, 1.0,
							  0.0, 1.0) ;
    
    for (i = 0; i < 100000 ; i++)
      {
	double u = ((double) rand ()) / (1+RAND_MAX);
	double v = ((double) rand ()) / (1+RAND_MAX);
	double x, y ;
	status = gsl_histogram2d_pdf_sample (p, u, v, &x, &y) ;
	status = gsl_histogram2d_increment (hh, x, y) ;
      }

    status = 0 ;
    for (i = 0 ; i < 20 ; i++)
      {
	for (j = 0; j < 20; j++)
	  {
	    double z = gsl_histogram2d_get (hh, i, j) / (100000.0 / 16.0)  ;
	    double x = gsl_histogram2d_get_xlowerlimit (hh, i) ;
	    double y = gsl_histogram2d_get_ylowerlimit (hh, j) ;
	    size_t k1, k2 ; double ya ;
	    
	    gsl_histogram2d_find (h, x, y, &k1, &k2) ;
	    ya = gsl_histogram2d_get (h, k1, k2) ;

	    if (ya == 0) {
	      if (z != 0) 
		{
		  status = 1 ;
		  printf("(%d,%d): %g vs %g\n", i,j, z, ya) ;
		}
	    } else {
	      double err = 1/sqrt(gsl_histogram2d_get (hh, i, j)) ;
	      double sigma = fabs((z-ya)/(ya*err)) ;
	      if (sigma > 3)
		{
		  status = 1 ;
		  printf("%g vs %g err=%g sigma=%g\n", z, ya, err, sigma) ;
		}
	    }
	  }
      }
    gsl_test(status, "gsl_histogram2d_pdf_sample works") ;
  }
  return gsl_test_summary () ;
}


