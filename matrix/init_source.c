#include <stdlib.h>
#include <gsl_matrix.h>


TYPE(gsl_matrix) * 
FUNCTION(gsl_matrix,alloc) (const size_t n1, const size_t n2)
{
  TYPE(gsl_matrix) * m ;

  if (n1 == 0)
    {
      GSL_ERROR_RETURN ("matrix dimension n1 must be positive integer", 
			GSL_EDOM, 0);
    }
  else if (n2 == 0)
    {
      GSL_ERROR_RETURN ("matrix dimension n2 must be positive integer", 
			GSL_EDOM, 0);
    }

  m = (TYPE(gsl_matrix) *) malloc(sizeof(TYPE(gsl_matrix))) ;
  
  if (m == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for matrix struct",
			GSL_ENOMEM, 0);
    }

  m->data = (ATOMIC *) malloc(MULTIPLICITY * n1 * n2 * sizeof(ATOMIC)) ;

  if (m->data == 0) 
    {
      free(m) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for matrix data", 
			GSL_ENOMEM, 0);
    }
  
  m->size1 = n1 ;
  m->size2 = n2 ;

  return m ;
}

TYPE(gsl_matrix) *
FUNCTION(gsl_matrix,calloc) (const size_t n1, const size_t n2)
{
  size_t i ;

  TYPE(gsl_matrix) * m = FUNCTION(gsl_matrix,alloc) (n1, n2) ;
  
  if (m == 0) 
    return 0 ;

  /* initialize matrix to zero */

  for (i = 0 ; i < MULTIPLICITY * n1 * n2; i++)  
    {
      m->data[i] = 0.0 ;
    }

  return m ;
}

void
FUNCTION(gsl_matrix,free) (TYPE(gsl_matrix) * m)
{
  free(m->data) ;
  free(m) ;
}

