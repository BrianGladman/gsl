#ifndef GSL_MATRIX_INT_H
#define GSL_MATRIX_INT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_vector_int.h>

typedef struct gsl_matrix_int_struct gsl_matrix_int;

struct gsl_matrix_int_struct
{
  size_t size1;
  size_t size2;
  size_t dim2;
  int * data;
  gsl_block_int * block;
} ;

gsl_matrix_int * 
gsl_matrix_int_alloc (size_t n1, size_t n2);

gsl_matrix_int * 
gsl_matrix_int_calloc (size_t n1, size_t n2);

gsl_matrix_int * 
gsl_matrix_int_alloc_from_block (gsl_block_int * b, size_t offset, 
                                   size_t n1, size_t n2, size_t d2);

gsl_matrix_int * 
gsl_matrix_int_alloc_from_matrix (gsl_matrix_int * m,
                                    size_t k1, size_t k2,
                                    size_t n1, size_t n2);

gsl_vector_int * 
gsl_vector_int_alloc_row_from_matrix (gsl_matrix_int * m,
                                        size_t i);

gsl_vector_int * 
gsl_vector_int_alloc_col_from_matrix (gsl_matrix_int * m,
                                        size_t j);

void gsl_matrix_int_free (gsl_matrix_int * m);

int * gsl_matrix_int_ptr(const gsl_matrix_int * m, size_t i, size_t j);
int   gsl_matrix_int_get(const gsl_matrix_int * m, size_t i, size_t j);
void    gsl_matrix_int_set(gsl_matrix_int * m, size_t i, size_t j, int x);

int gsl_matrix_int_fread (FILE * stream, gsl_matrix_int * m) ;
int gsl_matrix_int_fwrite (FILE * stream, const gsl_matrix_int * m) ;
int gsl_matrix_int_fscanf (FILE * stream, gsl_matrix_int * m);
int gsl_matrix_int_fprintf (FILE * stream, const gsl_matrix_int * m, const char * format);
 
int gsl_matrix_int_copy_row(gsl_vector_int * v, const gsl_matrix_int * m, size_t i);
int gsl_matrix_int_copy_col(gsl_vector_int * v, const gsl_matrix_int * m, size_t j);
int gsl_matrix_int_set_row(gsl_matrix_int * m, size_t i, const gsl_vector_int * v);
int gsl_matrix_int_set_col(gsl_matrix_int * m, size_t j, const gsl_vector_int * v);

int gsl_matrix_int_swap_rows(gsl_matrix_int * m, size_t i, size_t j);
int gsl_matrix_int_swap_cols(gsl_matrix_int * m, size_t i, size_t j);
int gsl_matrix_int_swap_rowcol(gsl_matrix_int * m, size_t i, size_t j);

int gsl_matrix_int_copy(gsl_matrix_int * dest, const gsl_matrix_int * src);

int gsl_vector_int_view_row_from_matrix (gsl_vector_int * v, gsl_matrix_int * m, size_t i);
int gsl_vector_int_view_col_from_matrix (gsl_vector_int * v, gsl_matrix_int * m, size_t j);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
int
gsl_matrix_int_get(const gsl_matrix_int * m, 
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
gsl_matrix_int_set(gsl_matrix_int * m, 
		     const size_t i, const size_t j, const int x)
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

#endif /* GSL_MATRIX_INT_H */
