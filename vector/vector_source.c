#include <gsl_errno.h>

#include "source.h"

BASE
FUNCTION(gsl_vector,get)(const TYPE(gsl_vector) * v, const size_t i)
{
#ifdef GSL_CHECK_RANGE
  if (i < 0 || i >= v->size)  /* if size_t is unsigned i<0 is impossible! */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return v->data[i] ;
}

void
FUNCTION(gsl_vector,set)(TYPE(gsl_vector) * v, const size_t i, const BASE x)
{
#ifdef GSL_CHECK_RANGE
  if (i < 0 || i >= v->size) /* if size_t is unsigned  i<0 is impossible! */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, /* nothing */) ;
    }
#endif
  v->data[i] = x ;
}

