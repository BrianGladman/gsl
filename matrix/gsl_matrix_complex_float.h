#ifndef GSL_MATRIX_COMPLEX_FLOAT_H
#define GSL_MATRIX_COMPLEX_FLOAT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_vector_complex_float.h>

typedef struct gsl_matrix_complex_float_struct gsl_matrix_complex_float;

struct gsl_matrix_complex_float_struct 
{
  size_t size1;
  size_t size2;
  size_t dim2;
  float * data;
  gsl_block_complex_float * block;
} ;


gsl_matrix_complex_float * 
gsl_matrix_complex_float_alloc (size_t n1, size_t n2);

gsl_matrix_complex_float * 
gsl_matrix_complex_float_calloc (size_t n1, size_t n2);

gsl_matrix_complex_float * 
gsl_matrix_complex_float_alloc_from_block (gsl_block_complex_float * b, 
                                           size_t offset, 
                                           size_t n1, size_t n2, size_t d2);

gsl_matrix_complex_float * 
gsl_matrix_complex_float_alloc_from_matrix (gsl_matrix_complex_float * b,
                                            size_t k1, size_t k2,
                                            size_t n1, size_t n2);

gsl_vector_complex_float * 
gsl_vector_complex_float_alloc_row_from_matrix (gsl_matrix_complex_float * m,
                                                size_t i);

gsl_vector_complex_float * 
gsl_vector_complex_float_alloc_col_from_matrix (gsl_matrix_complex_float * m,
                                                size_t j);


void gsl_matrix_complex_float_free (gsl_matrix_complex_float * m);

gsl_complex_float * gsl_matrix_complex_float_ptr(const gsl_matrix_complex_float * m, size_t i, size_t j);
gsl_complex_float gsl_matrix_complex_float_get(const gsl_matrix_complex_float * m, size_t i, size_t j);
void gsl_matrix_complex_float_set(gsl_matrix_complex_float * m, size_t i, size_t j, gsl_complex_float x);

int gsl_matrix_complex_float_fread (FILE * stream, gsl_matrix_complex_float * m) ;
int gsl_matrix_complex_float_fwrite (FILE * stream, const gsl_matrix_complex_float * m) ;
int gsl_matrix_complex_float_fscanf (FILE * stream, gsl_matrix_complex_float * m);
int gsl_matrix_complex_float_fprintf (FILE * stream, const gsl_matrix_complex_float * m, const char * format);

int gsl_matrix_complex_float_copy_row(gsl_vector_complex_float * v, const gsl_matrix_complex_float * m, size_t i);
int gsl_matrix_complex_float_copy_col(gsl_vector_complex_float * v, const gsl_matrix_complex_float * m, size_t j);
int gsl_matrix_complex_float_set_row(gsl_matrix_complex_float * m, size_t i, const gsl_vector_complex_float * v);
int gsl_matrix_complex_float_set_col(gsl_matrix_complex_float * m, size_t j, const gsl_vector_complex_float * v);

int gsl_matrix_complex_float_swap_rows(gsl_matrix_complex_float * m, size_t i, size_t j);
int gsl_matrix_complex_float_swap_cols(gsl_matrix_complex_float * m, size_t i, size_t j);
int gsl_matrix_complex_float_swap_rowcol(gsl_matrix_complex_float * m, size_t i, size_t j);

int gsl_matrix_complex_float_copy(gsl_matrix_complex_float * dest, const gsl_matrix_complex_float * src);

extern int gsl_check_range ;

#ifdef HAVE_INLINE

extern inline 
gsl_complex_float
gsl_matrix_complex_float_get(const gsl_matrix_complex_float * m, 
		     const size_t i, const size_t j)
{
  const gsl_complex_float zero = {{0,0}};

#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1)
    {
      GSL_ERROR_RETURN("first index out of range", GSL_EINVAL, zero) ;
    }
  else if (j >= m->size2)
    {
      GSL_ERROR_RETURN("second index out of range", GSL_EINVAL, zero) ;
    }
#endif
  return *(gsl_complex_float *)(m->data + 2*(i * m->dim2 + j)) ;
} 

extern inline 
void
gsl_matrix_complex_float_set(gsl_matrix_complex_float * m, 
		     const size_t i, const size_t j, const gsl_complex_float x)
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
  *(gsl_complex_float *)(m->data + 2*(i * m->dim2 + j)) = x ;
}
#endif /* HAVE_INLINE */

#endif /* !GSL_MATRIX_COMPLEX_FLOAT_H */
