#include <gsl_errno.h>

#include "source.h"

int gsl_check_range = 0; 

BASE
FUNCTION(gsl_matrix,get)(const TYPE(gsl_matrix) * m, 
			 const size_t i, const size_t j)
{
#ifdef GSL_CHECK_RANGE
  if (gsl_check_range) 
    {
      if (i >= m->size1)  /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN("first index out of range", GSL_EINVAL, 0) ;
	}
      else if (j >= m->size2) /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN("second index out of range", GSL_EINVAL, 0) ;
	}
    }
#endif
  return m->data[i * m->size2 + j] ;
} 

void
FUNCTION(gsl_matrix,set)(TYPE(gsl_matrix) * m,
			 const size_t i, const size_t j, 
			 const BASE x)
{
#ifdef GSL_CHECK_RANGE
  if (gsl_check_range) 
    {
      if (i >= m->size) /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN("first index out of range", GSL_EINVAL, /* void */);
	}
      else if (j >= m->size2) /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN("second index out of range", GSL_EINVAL,/* void */);
	}
    }
#endif
  m->data[i * m->size2 + j] = x ;
}

