#include <gsl_errno.h>
#include <gsl_vector.h>

double
gsl_vector_get(const gsl_vector * v, const size_t i)
{
#ifdef GSL_CHECK_RANGE
  if (i < 0 || i >= v->size)  /* if size_t is unsigned i<0 is impossible! */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0.0) ;
    }
#endif
  return v->data[i] ;
}

void
gsl_vector_set(gsl_vector * v, const size_t i, const double x)
{
#ifdef GSL_CHECK_RANGE
  if (i < 0 || i >= v->size) /* if size_t is unsigned  i<0 is impossible! */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, /* nothing */) ;
    }
#endif
  v->data[i] = x ;
}

