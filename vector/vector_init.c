#include <stdlib.h>
#include <gsl_vector.h>

gsl_vector * 
gsl_vector_alloc (size_t n)
{
  gsl_vector * v ;

  if (n == 0)
    {
      GSL_ERROR_RETURN ("vector length n must be positive integer", 
			GSL_EDOM, 0);
    }

  v = (gsl_vector *) malloc(sizeof(gsl_vector)) ;
  
  if (v == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
    }
    
  v->data = malloc(n * sizeof(double)) ;

  if (v->data == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for vector data", 
			GSL_ENOMEM, 0);
    }
  
  v->size = n ;

  return v ;
}

gsl_vector *
gsl_vector_calloc (size_t n)
{
  size_t i ;

  gsl_vector * v = gsl_vector_alloc (n) ;
  
  if (v == 0) 
    return 0 ;

  for (i = 0 ; i < n; i++)  /* initialize vector to zero */
    {
      v->data[i] = 0.0 ;
    }

  return v ;
}


void
gsl_vector_free (gsl_vector * v)
{
  free(v->data) ;
  free(v) ;
}

