#include <config.h>
#include <gsl_errno.h>

#include "source.h"

inline BASE
FUNCTION(gsl_vector,get)(const TYPE(gsl_vector) * v, const size_t i)
{
  if (gsl_check_range) 
    {
      if (i >= v->size) /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
	}
    }

  return v->data[i] ;
}

inline void
FUNCTION(gsl_vector,set)(TYPE(gsl_vector) * v, const size_t i, const BASE x)
{
  if (gsl_check_range) 
    {
      if (i >= v->size) /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN_NOTHING("index out of range", GSL_EINVAL) ;
	}
    }
  v->data[i] = x ;
}

