#include <stdlib.h>
#include <gsl_errno.h>

#include "source.h"

BASE * 
FUNCTION(gsl_monte_vector,alloc) (const size_t n)
{
  BASE * v ;

  if (n == 0)
    {
      GSL_ERROR_RETURN ("vector length n must be positive integer", 
			GSL_EDOM, 0) ;
    }

  v = (BASE *) malloc(n * sizeof(BASE)) ;

  if (v == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for vector data", 
			GSL_ENOMEM, 0);
    }
  return v ;
}

BASE *
FUNCTION(gsl_monte_vector,calloc) (const size_t n)
{
  size_t i ;

  BASE * v = (BASE *) FUNCTION(gsl_vector,alloc) (n) ;
  
  if (v == 0) 
    return 0 ;

  for (i = 0 ; i < n; i++)  /* initialize vector to zero */
    {
      v[i] = ZERO ;
    }

  return v ;
}


void
FUNCTION(gsl_monte_vector,free) (BASE * v)
{
  if ( v == (BASE *) NULL) {
    GSL_ERROR_RETURN_NOTHING("Attempt to free null pointer", EFAULT);
  }
  free(v) ;
}

