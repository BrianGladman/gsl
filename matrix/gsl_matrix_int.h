#include <stdlib.h>
#include <gsl_errno.h>

typedef struct
{
  size_t size1;
  size_t size2;
  int * data;
  
} gsl_matrix_int ;

gsl_matrix_int * gsl_matrix_int_alloc (size_t n1, size_t n2);
gsl_matrix_int * gsl_matrix_int_calloc (size_t n1, size_t n2);
void gsl_matrix_int_free (gsl_matrix_int * m);

int gsl_matrix_int_get(const gsl_matrix_int * m, size_t i, size_t j);
void gsl_matrix_int_set(gsl_matrix_int * m, size_t i,  size_t j, int x);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifndef __STRICT_ANSI__
extern inline 
int
gsl_matrix_int_get(const gsl_matrix_int * m, const size_t i, const size_t j)
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
gsl_matrix_int_set(gsl_matrix_int * m, 
		   const size_t i, const size_t j, 
		   const int x)
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

