#ifndef GSL_MATRIX_CHAR_H
#define GSL_MATRIX_CHAR_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector_char.h>

typedef struct gsl_matrix_char_struct gsl_matrix_char;

struct gsl_matrix_char_struct
{
  size_t size1;
  size_t size2;
  size_t dim2;
  char * data;
  gsl_block_char * block;
} ;

gsl_matrix_char * 
gsl_matrix_char_alloc (size_t n1, size_t n2);

gsl_matrix_char * 
gsl_matrix_char_calloc (size_t n1, size_t n2);

gsl_matrix_char * 
gsl_matrix_char_alloc_from_block (gsl_block_char * b, size_t offset, 
                                   size_t n1, size_t n2, size_t d2);

gsl_matrix_char * 
gsl_matrix_char_alloc_from_matrix (gsl_matrix_char * m,
                                    size_t k1, size_t k2,
                                    size_t n1, size_t n2);

gsl_vector_char * 
gsl_vector_char_alloc_row_from_matrix (gsl_matrix_char * m,
                                        size_t i);

gsl_vector_char * 
gsl_vector_char_alloc_col_from_matrix (gsl_matrix_char * m,
                                        size_t j);

void gsl_matrix_char_free (gsl_matrix_char * m);

char * gsl_matrix_char_ptr(const gsl_matrix_char * m, size_t i, size_t j);
char   gsl_matrix_char_get(const gsl_matrix_char * m, size_t i, size_t j);
void    gsl_matrix_char_set(gsl_matrix_char * m, size_t i, size_t j, char x);

int gsl_matrix_char_fread (FILE * stream, gsl_matrix_char * m) ;
int gsl_matrix_char_fwrite (FILE * stream, const gsl_matrix_char * m) ;
int gsl_matrix_char_fscanf (FILE * stream, gsl_matrix_char * m);
int gsl_matrix_char_fprintf (FILE * stream, const gsl_matrix_char * m, const char * format);
 
int gsl_matrix_char_copy_row(gsl_vector_char * v, const gsl_matrix_char * m, size_t i);
int gsl_matrix_char_copy_col(gsl_vector_char * v, const gsl_matrix_char * m, size_t j);
int gsl_matrix_char_set_row(gsl_matrix_char * m, size_t i, const gsl_vector_char * v);
int gsl_matrix_char_set_col(gsl_matrix_char * m, size_t j, const gsl_vector_char * v);

int gsl_matrix_char_swap_rows(gsl_matrix_char * m, size_t i, size_t j);
int gsl_matrix_char_swap_cols(gsl_matrix_char * m, size_t i, size_t j);
int gsl_matrix_char_swap_rowcol(gsl_matrix_char * m, size_t i, size_t j);

int gsl_matrix_char_copy(gsl_matrix_char * dest, const gsl_matrix_char * src);

int gsl_vector_char_view_row_from_matrix (gsl_vector_char * v, gsl_matrix_char * m, size_t i);
int gsl_vector_char_view_col_from_matrix (gsl_vector_char * v, gsl_matrix_char * m, size_t j);

int gsl_matrix_char_view_from_vector (gsl_matrix_char * m, 
                                       gsl_vector_char * base,
                                       size_t offset, 
                                       size_t n1, size_t n2, size_t d2);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
char
gsl_matrix_char_get(const gsl_matrix_char * m, 
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
gsl_matrix_char_set(gsl_matrix_char * m, 
		     const size_t i, const size_t j, const char x)
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

#endif /* GSL_MATRIX_CHAR_H */
