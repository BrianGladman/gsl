#include <gsl_errno.h>

#include "source.h"

inline BASE
FUNCTION(gsl_matrix,get)(const TYPE(gsl_matrix) * m, 
			 const size_t i, const size_t j)
{
  if (gsl_check_range) 
    {
      if (i >= m->size1)  /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_CONTINUE("first index out of range", GSL_EINVAL) ;
	}
      else if (j >= m->size2) /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_CONTINUE("second index out of range", GSL_EINVAL) ;
	}
    }
  return *(BASE *)(m->data + MULTIPLICITY*(i * m->size2 + j)) ;
} 

inline void
FUNCTION(gsl_matrix,set)(TYPE(gsl_matrix) * m,
			 const size_t i, const size_t j, 
			 const BASE x)
{
  if (gsl_check_range) 
    {
      if (i >= m->size1) /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN_NOTHING("first index out of range", GSL_EINVAL);
	}
      else if (j >= m->size2) /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN_NOTHING("second index out of range", GSL_EINVAL);
	}
    }
  *(BASE *)(m->data + MULTIPLICITY*(i * m->size2 + j)) = x ;
}

