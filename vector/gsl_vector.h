#include <stdlib.h>
#include <gsl_errno.h>

typedef struct
{
  size_t size;
  double * data;
} gsl_vector ;

gsl_vector * gsl_vector_alloc (size_t n);
gsl_vector * gsl_vector_calloc (size_t n);
void gsl_vector_free (gsl_vector * v);

double gsl_vector_get(const gsl_vector * v, const size_t i);
void gsl_vector_set(gsl_vector * v, const size_t i, const double x);

/* inline functions if you are using GCC */

#ifndef __STRICT_ANSI__
extern inline 
double
gsl_vector_get(const gsl_vector * v, const size_t i)
{
#ifdef GSL_CHECK_RANGE
  if (i < 0 || i >= v->size)  /* if size_t is unsigned i<0 is impossible! */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return v->data[i] ;
} 

extern inline 
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
#endif

