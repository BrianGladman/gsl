#ifndef GSL_MATRIX_FLOAT_H
#define GSL_MATRIX_FLOAT_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector_float.h>

typedef struct gsl_matrix_float_struct gsl_matrix_float;

struct gsl_matrix_float_struct
{
  size_t size1;
  size_t size2;
  size_t dim2;
  float * data;
  gsl_block_float * block;
} ;

gsl_matrix_float * 
gsl_matrix_float_alloc (size_t n1, size_t n2);

gsl_matrix_float * 
gsl_matrix_float_calloc (size_t n1, size_t n2);

gsl_matrix_float * 
gsl_matrix_float_alloc_from_block (gsl_block_float * b, size_t offset, 
                                   size_t n1, size_t n2, size_t d2);

gsl_matrix_float * 
gsl_matrix_float_alloc_from_matrix (gsl_matrix_float * m,
                                    size_t k1, size_t k2,
                                    size_t n1, size_t n2);

gsl_vector_float * 
gsl_vector_float_alloc_row_from_matrix (gsl_matrix_float * m,
                                        size_t i);

gsl_vector_float * 
gsl_vector_float_alloc_col_from_matrix (gsl_matrix_float * m,
                                        size_t j);

void gsl_matrix_float_free (gsl_matrix_float * m);

float * gsl_matrix_float_ptr(const gsl_matrix_float * m, size_t i, size_t j);
float   gsl_matrix_float_get(const gsl_matrix_float * m, size_t i, size_t j);
void    gsl_matrix_float_set(gsl_matrix_float * m, size_t i, size_t j, float x);

int gsl_matrix_float_fread (FILE * stream, gsl_matrix_float * m) ;
int gsl_matrix_float_fwrite (FILE * stream, const gsl_matrix_float * m) ;
int gsl_matrix_float_fscanf (FILE * stream, gsl_matrix_float * m);
int gsl_matrix_float_fprintf (FILE * stream, const gsl_matrix_float * m, const char * format);
 
int gsl_matrix_float_copy_row(gsl_vector_float * v, const gsl_matrix_float * m, size_t i);
int gsl_matrix_float_copy_col(gsl_vector_float * v, const gsl_matrix_float * m, size_t j);
int gsl_matrix_float_set_row(gsl_matrix_float * m, size_t i, const gsl_vector_float * v);
int gsl_matrix_float_set_col(gsl_matrix_float * m, size_t j, const gsl_vector_float * v);

int gsl_matrix_float_swap_rows(gsl_matrix_float * m, size_t i, size_t j);
int gsl_matrix_float_swap_cols(gsl_matrix_float * m, size_t i, size_t j);
int gsl_matrix_float_swap_rowcol(gsl_matrix_float * m, size_t i, size_t j);

int gsl_matrix_float_copy(gsl_matrix_float * dest, const gsl_matrix_float * src);

int gsl_vector_float_view_row_from_matrix (gsl_vector_float * v, gsl_matrix_float * m, size_t i);
int gsl_vector_float_view_col_from_matrix (gsl_vector_float * v, gsl_matrix_float * m, size_t j);

int gsl_matrix_float_view_from_vector (gsl_matrix_float * m, 
                                       gsl_vector_float * base,
                                       size_t offset, 
                                       size_t n1, size_t n2, size_t d2);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
float
gsl_matrix_float_get(const gsl_matrix_float * m, 
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
gsl_matrix_float_set(gsl_matrix_float * m, 
		     const size_t i, const size_t j, const float x)
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

#endif /* GSL_MATRIX_FLOAT_H */
