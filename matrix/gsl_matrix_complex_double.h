#ifndef GSL_MATRIX_COMPLEX_DOUBLE_H
#define GSL_MATRIX_COMPLEX_DOUBLE_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_vector_complex_double.h>

typedef struct gsl_matrix_complex_struct gsl_matrix_complex;

struct gsl_matrix_complex_struct 
{
  size_t size1;
  size_t size2;
  size_t dim2;
  double * data;
  gsl_block_complex * block;
} ;


gsl_matrix_complex * 
gsl_matrix_complex_alloc (size_t n1, size_t n2);

gsl_matrix_complex * 
gsl_matrix_complex_calloc (size_t n1, size_t n2);

gsl_matrix_complex * 
gsl_matrix_complex_alloc_from_block (gsl_block_complex * b, 
                                           size_t offset, 
                                           size_t n1, size_t n2, size_t d2);

gsl_matrix_complex * 
gsl_matrix_complex_alloc_from_matrix (gsl_matrix_complex * b,
                                            size_t k1, size_t k2,
                                            size_t n1, size_t n2);

gsl_vector_complex * 
gsl_vector_complex_alloc_row_from_matrix (gsl_matrix_complex * m,
                                                size_t i);

gsl_vector_complex * 
gsl_vector_complex_alloc_col_from_matrix (gsl_matrix_complex * m,
                                                size_t j);


void gsl_matrix_complex_free (gsl_matrix_complex * m);

gsl_complex * gsl_matrix_complex_ptr(const gsl_matrix_complex * m, size_t i, size_t j);
gsl_complex gsl_matrix_complex_get(const gsl_matrix_complex * m, size_t i, size_t j);
void gsl_matrix_complex_set(gsl_matrix_complex * m, size_t i, size_t j, gsl_complex x);

int gsl_matrix_complex_fread (FILE * stream, gsl_matrix_complex * m) ;
int gsl_matrix_complex_fwrite (FILE * stream, const gsl_matrix_complex * m) ;
int gsl_matrix_complex_fscanf (FILE * stream, gsl_matrix_complex * m);
int gsl_matrix_complex_fprintf (FILE * stream, const gsl_matrix_complex * m, const char * format);

int gsl_matrix_complex_copy_row(gsl_vector_complex * v, const gsl_matrix_complex * m, size_t i);
int gsl_matrix_complex_copy_col(gsl_vector_complex * v, const gsl_matrix_complex * m, size_t j);
int gsl_matrix_complex_set_row(gsl_matrix_complex * m, size_t i, const gsl_vector_complex * v);
int gsl_matrix_complex_set_col(gsl_matrix_complex * m, size_t j, const gsl_vector_complex * v);

int gsl_matrix_complex_swap_rows(gsl_matrix_complex * m, size_t i, size_t j);
int gsl_matrix_complex_swap_cols(gsl_matrix_complex * m, size_t i, size_t j);
int gsl_matrix_complex_swap_rowcol(gsl_matrix_complex * m, size_t i, size_t j);

int gsl_matrix_complex_copy(gsl_matrix_complex * dest, const gsl_matrix_complex * src);

extern int gsl_check_range ;

#ifdef HAVE_INLINE

extern inline 
gsl_complex
gsl_matrix_complex_get(const gsl_matrix_complex * m, 
		     const size_t i, const size_t j)
{
  const gsl_complex zero = {{0,0}};

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
  return *(gsl_complex *)(m->data + 2*(i * m->dim2 + j)) ;
} 

extern inline 
void
gsl_matrix_complex_set(gsl_matrix_complex * m, 
		     const size_t i, const size_t j, const gsl_complex x)
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
  *(gsl_complex *)(m->data + 2*(i * m->dim2 + j)) = x ;
}
#endif /* HAVE_INLINE */

#endif /* !GSL_MATRIX_COMPLEX_DOUBLE_H */
