#include <config.h>
#include <gsl_errno.h>
#include <gsl_histogram.h>

int
gsl_histogram_get (const gsl_histogram * h, size_t i, double * y)
{
  const size_t n = h->nbins ;

  if (i > n) 
    {
      GSL_ERROR("index lies outside valid range of 0 .. nbins - 1", EDOM) ;
    }
  
  *y = h->bin[i] ;

  return 0 ;
}

int
gsl_histogram_get_binrange (const gsl_histogram * h, size_t i,
			    double * x0, double * x1)
{
  const size_t n = h->nbins ;

  if (i > n) 
    {
      GSL_ERROR("index lies outside valid range of 0 .. nbins - 1", EDOM) ;
    }
  
  *x0 = h->range[i] ;
  *x1 = h->range[i + 1] ;

  return 0 ;
}

