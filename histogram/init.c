#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_histogram.h>

gsl_histogram *
gsl_histogram_alloc_uniform (const size_t n, const double xmin,
			     const double xmax)
{
  gsl_histogram * h ;
  
  if (xmin >= xmax) 
    {
      GSL_ERROR_RETURN ("xmin must be less than xmax", GSL_EINVAL, 0) ;
    }
  
  h = gsl_histogram_alloc (n) ;
  
  if (h == 0) 
    {
      return h ;
    }
    
  {
    size_t i ;
    
    for (i = 0; i < n + 1; i++) 
      {
	h->range[i] = xmin + ((double)i/(double)n) * (xmax - xmin) ;
      }
  }

  return h ;
}

gsl_histogram *
gsl_histogram_alloc (size_t n)
{
  gsl_histogram * h ;
 
  if (n == 0)
    {
      GSL_ERROR_RETURN ("histogram length n must be positive integer", 
			GSL_EDOM, 0) ;
    }

  h = (gsl_histogram *) malloc(sizeof(gsl_histogram)) ;

  if (h == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for histogram struct",
			GSL_ENOMEM, 0);
    }

  h->range = (double *) malloc((n + 1) * sizeof(double)) ;

  if (h->range == 0) 
    {
      free(h) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for histogram ranges", 
			GSL_ENOMEM, 0);
    }

  h->bin = (double *) malloc(n * sizeof(double)) ;

  if (h->bin == 0) 
    {
      free(h->range) ;
      free(h) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for histogram bins", 
			GSL_ENOMEM, 0);
    }

  {
    size_t i ;

    for (i = 0; i < n + 1; i++)
      {
	h->range[i] = i ;
      }
    
    for (i = 0; i < n; i++)
      {
	h->bin[i] = 0 ;
      }
  }

  h->nbins = n;

  return h ;
}


void
gsl_histogram_free (gsl_histogram * h)
{
  free (h->range);
  free (h->bin);
  free (h);
}
