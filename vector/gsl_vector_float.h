#ifndef GSL_VECTOR_FLOAT_H 
#define GSL_VECTOR_FLOAT_H 

#include <stdlib.h>
#include <gsl_errno.h>

typedef struct
{
  size_t size;
  float * data;
} gsl_vector_float ;

gsl_vector_float * gsl_vector_float_alloc (size_t n);
gsl_vector_float * gsl_vector_float_calloc (size_t n);
void gsl_vector_float_free (gsl_vector_float * v);

float gsl_vector_float_get(const gsl_vector_float * v, const size_t i);
void gsl_vector_float_set(gsl_vector_float * v, const size_t i, const float x);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifndef __STRICT_ANSI__
extern inline 
float
gsl_vector_float_get(const gsl_vector_float * v, const size_t i)
{
#ifdef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return v->data[i] ;
} 

extern inline 
void
gsl_vector_float_set(gsl_vector_float * v, const size_t i, const float x)
{
#ifdef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, /* nothing */) ;
    }
#endif
  v->data[i] = x ;
}
#endif

#endif /* GSL_VECTOR_FLOAT_H */
