#ifndef GSL_VECTOR_INT_H 
#define GSL_VECTOR_INT_H 

#include <stdlib.h>
#include <gsl_errno.h>

typedef struct
{
  size_t size;
  int * data;
} gsl_vector_int ;

gsl_vector_int * gsl_vector_int_alloc (size_t n);
gsl_vector_int * gsl_vector_int_calloc (size_t n);
void gsl_vector_int_free (gsl_vector_int * v);

int gsl_vector_int_get(const gsl_vector_int * v, const size_t i);
void gsl_vector_int_set(gsl_vector_int * v, const size_t i, const int x);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifndef __STRICT_ANSI__
extern inline 
int
gsl_vector_int_get(const gsl_vector_int * v, const size_t i)
{
#ifdef GSL_CHECK_RANGE
  if (i >= v->size)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return v->data[i] ;
} 

extern inline 
void
gsl_vector_int_set(gsl_vector_int * v, const size_t i, const int x)
{
#ifdef GSL_CHECK_RANGE
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, /* nothing */) ;
    }
#endif
  v->data[i] = x ;
}
#endif

#endif /* GSL_VECTOR_INT_H */
