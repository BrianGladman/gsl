#include <gsl_errno.h>
#include "source.h"

inline BASE *
FUNCTION(gsl_vector,ptr)(const TYPE(gsl_vector) * v, const size_t i)
{
  if (gsl_check_range) 
    {
      if (i >= v->size) /* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
	}
    }

  return (BASE *) (v->data + i) ;
}


#if WANT_GET_PROTO==1

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

  return v->data[i];
}

#endif /* WANT_GET_PROTO */
