#include <stdlib.h>
#include <gsl_errno.h>

typedef struct
{
  size_t size1;
  size_t size2;
  double * data;
  
} gsl_matrix ;

gsl_matrix * gsl_matrix_alloc (size_t n1, size_t n2);
gsl_matrix * gsl_matrix_calloc (size_t n1, size_t n2);
void gsl_matrix_free (gsl_matrix * m);

double gsl_matrix_get(const gsl_matrix * m, size_t i, size_t j);
void gsl_matrix_set(gsl_matrix * m, size_t i,  size_t j, double x);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifndef __STRICT_ANSI__
extern inline 
double
gsl_matrix_get(const gsl_matrix * m, const size_t i, const size_t j)
{
#ifndef GSL_RANGE_CHECK_OFF
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
gsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("first index out of range", GSL_EINVAL, /* nothing */) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("second index out of range", GSL_EINVAL, /* nothing */) ;
    }
#endif
  m->data[i * m->size2 + j] = x ;
}
#endif

