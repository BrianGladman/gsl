#include <config.h>
#include <gsl_errno.h>
#include <gsl_histogram2d.h>

double
gsl_histogram2d_get (const gsl_histogram2d * h, const size_t i, const size_t j)
{
  const size_t nx = h->nx ;
  const size_t ny = h->ny ;

  if (i >= nx) 
    {
      GSL_ERROR_RETURN ("index i lies outside valid range of 0 .. nx - 1",
			GSL_EDOM, 0) ;
    }

  if (j >= ny) 
    {
      GSL_ERROR_RETURN ("index j lies outside valid range of 0 .. ny - 1",
			GSL_EDOM, 0) ;
    }
  
  return h->bin[i * ny + j] ;
}

double
gsl_histogram2d_get_xlowerlimit (const gsl_histogram2d * h, const size_t i)
{
  const size_t nx = h->nx ;

  if (i >= nx) 
    {
      GSL_ERROR_RETURN ("index i lies outside valid range of 0 .. nx - 1",
			GSL_EDOM, 0) ;
    }

  return h->xrange[i] ;
}

double
gsl_histogram2d_get_xupperlimit (const gsl_histogram2d * h, const size_t i)
{
  const size_t nx = h->nx ;

  if (i >= nx) 
    {
      GSL_ERROR_RETURN ("index i lies outside valid range of 0 .. nx - 1",
			GSL_EDOM, 0) ;
    }

  return h->xrange[i + 1] ;
}



double
gsl_histogram2d_get_ylowerlimit (const gsl_histogram2d * h, const size_t j)
{
  const size_t ny = h->ny ;

  if (j >= ny) 
    {
      GSL_ERROR_RETURN ("index j lies outside valid range of 0 .. ny - 1",
			GSL_EDOM, 0) ;
    }

  return h->yrange[j] ;
}

double
gsl_histogram2d_get_yupperlimit (const gsl_histogram2d * h, const size_t j)
{
  const size_t ny = h->ny ;

  if (j >= ny) 
    {
      GSL_ERROR_RETURN ("index j lies outside valid range of 0 .. ny - 1",
			GSL_EDOM, 0) ;
    }

  return h->yrange[j + 1] ;
}

