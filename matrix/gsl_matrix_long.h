#ifndef GSL_MATRIX_LONG_H
#define GSL_MATRIX_LONG_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_vector_long.h>

typedef struct gsl_matrix_long_struct gsl_matrix_long;

struct gsl_matrix_long_struct
{
  size_t size1;
  size_t size2;
  size_t dim2;
  long * data;
  gsl_block_long * block;
} ;

gsl_matrix_long * 
gsl_matrix_long_alloc (size_t n1, size_t n2);

gsl_matrix_long * 
gsl_matrix_long_calloc (size_t n1, size_t n2);

gsl_matrix_long * 
gsl_matrix_long_alloc_from_block (gsl_block_long * b, size_t offset, 
                                   size_t n1, size_t n2, size_t d2);

gsl_matrix_long * 
gsl_matrix_long_alloc_from_matrix (gsl_matrix_long * m,
                                    size_t k1, size_t k2,
                                    size_t n1, size_t n2);

gsl_vector_long * 
gsl_vector_long_alloc_row_from_matrix (gsl_matrix_long * m,
                                        size_t i);

gsl_vector_long * 
gsl_vector_long_alloc_col_from_matrix (gsl_matrix_long * m,
                                        size_t j);

void gsl_matrix_long_free (gsl_matrix_long * m);

long * gsl_matrix_long_ptr(const gsl_matrix_long * m, size_t i, size_t j);
long   gsl_matrix_long_get(const gsl_matrix_long * m, size_t i, size_t j);
void    gsl_matrix_long_set(gsl_matrix_long * m, size_t i, size_t j, long x);

int gsl_matrix_long_fread (FILE * stream, gsl_matrix_long * m) ;
int gsl_matrix_long_fwrite (FILE * stream, const gsl_matrix_long * m) ;
int gsl_matrix_long_fscanf (FILE * stream, gsl_matrix_long * m);
int gsl_matrix_long_fprintf (FILE * stream, const gsl_matrix_long * m, const char * format);
 
int gsl_matrix_long_copy_row(gsl_vector_long * v, const gsl_matrix_long * m, size_t i);
int gsl_matrix_long_copy_col(gsl_vector_long * v, const gsl_matrix_long * m, size_t j);
int gsl_matrix_long_set_row(gsl_matrix_long * m, size_t i, const gsl_vector_long * v);
int gsl_matrix_long_set_col(gsl_matrix_long * m, size_t j, const gsl_vector_long * v);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
long
gsl_matrix_long_get(const gsl_matrix_long * m, 
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
gsl_matrix_long_set(gsl_matrix_long * m, 
		     const size_t i, const size_t j, const long x)
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

#endif /* GSL_MATRIX_LONG_H */
