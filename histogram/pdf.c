#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_histogram.h>

double
gsl_histogram_pdf_sample (const gsl_histogram_pdf * p, const double r)
{
  size_t i ;
  int status = gsl_histogram_find (p->nbins, p->sum, r, &i) ;
  if (status) {
    return 0 ;
  } else {
    double delta = (r - p->sum[i])/(p->sum[i+1] - p->sum[i]) ;
    double x = p->range[i] + delta * (p->range[i+1] - p->range[i]) ;
    return x ;
  }
}

gsl_histogram_pdf *
gsl_histogram_pdf_alloc (const gsl_histogram * h)
{
  size_t i;
  size_t n = h->nbins ;
  gsl_histogram_pdf * p;
  
  if (n == 0)
    {
      GSL_ERROR_RETURN ("histogram length n must be positive integer", 
			GSL_EDOM, 0) ;
    }

  p = (gsl_histogram_pdf *) malloc(sizeof(gsl_histogram_pdf)) ;

  if (p == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for histogram pdf struct",
			GSL_ENOMEM, 0);
    }

  p->range = (double *) malloc((n + 1) * sizeof(double)) ;

  if (p->range == 0) 
    {
      free(p) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for histogram pdf ranges", 
			GSL_ENOMEM, 0);
    }

  p->sum = (double *) malloc((n + 1) * sizeof(double)) ;

  if (p->sum == 0) 
    {
      free(p->range) ;
      free(p) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for histogram pdf sums", 
			GSL_ENOMEM, 0);
    }

  for (i = 0; i < n + 1; i++)
    {
      p->range[i] = h->range[i] ;
    }

  {
    double total = 0, sum = 0;
    
    for (i = 0; i < n; i++)
      {
	const double b = h->bin[i] ;
	if (b >= 0) 
	  {
	    total += b ;
	  }
      }

    p->sum[0] = 0 ;

    for (i = 0; i < n; i++)
      {
	double b = h->bin[i] ;
	if (b >= 0) 
	  {
	    sum += b / total ;
	  }
	p->sum[i+1] = sum  ;
      }
  }

  p->nbins = n;

  return p ;
}

void
gsl_histogram_pdf_free (gsl_histogram_pdf * p)
{
  free (p->range);
  free (p->sum);
  free (p);
}
