#ifndef GSL_MATRIX_LONG_DOUBLE_H
#define GSL_MATRIX_LONG_DOUBLE_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_vector_long_double.h>

typedef struct gsl_matrix_long_double_struct gsl_matrix_long_double;

struct gsl_matrix_long_double_struct
{
  size_t size1;
  size_t size2;
  size_t dim2;
  long double * data;
} ;

gsl_matrix_long_double * 
gsl_matrix_long_double_alloc (gsl_block_long_double * b, size_t offset, 
                                   size_t n1, size_t n2, size_t d2);

gsl_matrix_long_double * 
gsl_matrix_long_double_alloc_from_matrix (gsl_matrix_long_double * b,
                                    size_t k1, size_t k2,
                                    size_t n1, size_t n2);

void gsl_matrix_long_double_free (gsl_matrix_long_double * m);

long double * gsl_matrix_long_double_ptr(const gsl_matrix_long_double * m, size_t i, size_t j);
long double   gsl_matrix_long_double_get(const gsl_matrix_long_double * m, size_t i, size_t j);
void    gsl_matrix_long_double_set(gsl_matrix_long_double * m, size_t i, size_t j, long double x);

int gsl_matrix_long_double_fread (FILE * stream, gsl_matrix_long_double * m) ;
int gsl_matrix_long_double_fwrite (FILE * stream, const gsl_matrix_long_double * m) ;
int gsl_matrix_long_double_fscanf (FILE * stream, gsl_matrix_long_double * m);
int gsl_matrix_long_double_fprintf (FILE * stream, const gsl_matrix_long_double * m, const char * format);
 
int gsl_matrix_long_double_copy_row(gsl_vector_long_double * v, const gsl_matrix_long_double * m, size_t i);
int gsl_matrix_long_double_copy_col(gsl_vector_long_double * v, const gsl_matrix_long_double * m, size_t j);
int gsl_matrix_long_double_set_row(gsl_matrix_long_double * m, size_t i, const gsl_vector_long_double * v);
int gsl_matrix_long_double_set_col(gsl_matrix_long_double * m, size_t j, const gsl_vector_long_double * v);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
long double
gsl_matrix_long_double_get(const gsl_matrix_long_double * m, 
		     const size_t i, const size_t j)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1)
    {
      GSL_ERROR_RETURN("first index out of range", GSL_EINVAL, 0) ;
    }
  else if (j >= m->size2)
    {
      GSL_ERROR_RETURN("second index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return m->data[i * m->dim2 + j] ;
} 

extern inline 
void
gsl_matrix_long_double_set(gsl_matrix_long_double * m, 
		     const size_t i, const size_t j, const long double x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1)
    {
      GSL_ERROR_RETURN_NOTHING("first index out of range", GSL_EINVAL) ;
    }
  else if (j >= m->size2)
    {
      GSL_ERROR_RETURN_NOTHING("second index out of range", GSL_EINVAL) ;
    }
#endif
  m->data[i * m->dim2 + j] = x ;
}
#endif

#endif /* GSL_MATRIX_LONG_DOUBLE_H */
