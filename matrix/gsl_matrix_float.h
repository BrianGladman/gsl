#include <stdlib.h>
#include <gsl_errno.h>

typedef struct
{
  size_t size1;
  size_t size2;
  float * data;
  
} gsl_matrix_float ;

gsl_matrix_float * gsl_matrix_float_alloc (size_t n1, size_t n2);
gsl_matrix_float * gsl_matrix_float_calloc (size_t n1, size_t n2);
void gsl_matrix_float_free (gsl_matrix_float * m);

float gsl_matrix_float_get(const gsl_matrix_float * m, size_t i, size_t j);
void gsl_matrix_float_set(gsl_matrix_float * m, size_t i, size_t j, float x);

/* inline functions if you are using GCC */

#ifndef __STRICT_ANSI__
extern inline 
float
gsl_matrix_float_get(const gsl_matrix_float * m, 
		     const size_t i, const size_t j)
{
#ifdef GSL_CHECK_RANGE
  if (i >= m->size1)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("first index out of range", GSL_EINVAL, 0) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("second index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return m->data[i * m->size2 + j] ;
} 

extern inline 
void
gsl_matrix_float_set(gsl_matrix_float * m, 
		     const size_t i, const size_t j, const float x)
{
#ifdef GSL_CHECK_RANGE
  if (i >= m->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING("first index out of range", GSL_EINVAL) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING("second index out of range", GSL_EINVAL) ;
    }
#endif
  m->data[i * m->size2 + j] = x ;
}
#endif

