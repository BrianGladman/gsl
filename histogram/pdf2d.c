#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_histogram.h>
#include <gsl_histogram2d.h>

int
gsl_histogram2d_pdf_sample (const gsl_histogram2d_pdf * p, 
			    const double r1, const double r2, 
			    double * x, double * y)
{
  size_t k ;
  int status = gsl_histogram_find_impl (p->nx * p->ny, p->sum, r1, &k) ;

  if (status) {
    return status ;
  } else {
    size_t i = k / p->ny ;
    size_t j = k - (i * p->ny) ;
    double delta = (r1 - p->sum[k])/(p->sum[k+1] - p->sum[k]) ;
    *x = p->xrange[i] + delta * (p->xrange[i+1] - p->xrange[i]) ;
    *y = p->yrange[j] + r2 * (p->yrange[j+1] - p->yrange[j]) ;
    return 0 ;
  }
}

gsl_histogram2d_pdf *
gsl_histogram2d_pdf_alloc (const gsl_histogram2d * h)
{
  size_t i;
  const size_t nx = h->nx ;
  const size_t ny = h->ny ;
  const size_t n = nx * ny ;

  gsl_histogram2d_pdf * p;
  
  if (n == 0)
    {
      GSL_ERROR_RETURN ("histogram2d length n must be positive integer", 
			GSL_EDOM, 0) ;
    }

  p = (gsl_histogram2d_pdf *) malloc(sizeof(gsl_histogram2d_pdf)) ;

  if (p == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for histogram2d pdf struct",
			GSL_ENOMEM, 0);
    }

  p->xrange = (double *) malloc((nx + 1) * sizeof(double)) ;

  if (p->xrange == 0) 
    {
      free(p) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for histogram2d pdf xranges", 
			GSL_ENOMEM, 0);
    }

  p->yrange = (double *) malloc((ny + 1) * sizeof(double)) ;

  if (p->yrange == 0) 
    {
      free(p->xrange) ;
      free(p) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for histogram2d pdf yranges", 
			GSL_ENOMEM, 0);
    }

  p->sum = (double *) malloc((n + 1) * sizeof(double)) ;

  if (p->sum == 0) 
    {
      free(p->yrange) ;
      free(p->xrange) ;
      free(p) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for histogram2d pdf sums", 
			GSL_ENOMEM, 0);
    }

  for (i = 0; i < nx + 1; i++)
    {
      p->xrange[i] = h->xrange[i] ;
    }

  for (i = 0; i < ny + 1; i++)
    {
      p->yrange[i] = h->yrange[i] ;
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

  p->nx = nx;
  p->ny = ny ;

  return p ;
}

void
gsl_histogram2d_pdf_free (gsl_histogram2d_pdf * p)
{
  free (p->xrange);
  free (p->yrange);
  free (p->sum);
  free (p);
}
