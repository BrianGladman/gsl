#include <stdlib.h>

#include "source.h"

TYPE(gsl_vector) * 
FUNCTION(gsl_vector,alloc) (const size_t n)
{
  TYPE(gsl_vector) * v ;

  if (n == 0)
    {
      GSL_ERROR_RETURN ("vector length n must be positive integer", 
			GSL_EDOM, 0);
    }

  v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector))) ;
  
  if (v == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
    }

  v->data = malloc(n * sizeof(BASE)) ;

  if (v->data == 0) 
    {
      free(v) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for vector data", 
			GSL_ENOMEM, 0);
    }
  
  v->size = n ;

  return v ;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector,calloc) (const size_t n)
{
  size_t i ;

  TYPE(gsl_vector) * v = FUNCTION(gsl_vector,alloc) (n) ;
  
  if (v == 0) 
    return 0 ;

  for (i = 0 ; i < n; i++)  /* initialize vector to zero */
    {
      v->data[i] = 0 ;
    }

  return v ;
}


void
FUNCTION(gsl_vector,free) (TYPE(gsl_vector) * v)
{
  free(v->data) ;
  free(v) ;
}

