#ifndef GSL_MATRIX_ULONG_H
#define GSL_MATRIX_ULONG_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector_ulong.h>

typedef struct gsl_matrix_ulong_struct gsl_matrix_ulong;

struct gsl_matrix_ulong_struct
{
  size_t size1;
  size_t size2;
  size_t dim2;
  unsigned long * data;
  gsl_block_ulong * block;
} ;

gsl_matrix_ulong * 
gsl_matrix_ulong_alloc (size_t n1, size_t n2);

gsl_matrix_ulong * 
gsl_matrix_ulong_calloc (size_t n1, size_t n2);

gsl_matrix_ulong * 
gsl_matrix_ulong_alloc_from_block (gsl_block_ulong * b, size_t offset, 
                                   size_t n1, size_t n2, size_t d2);

gsl_matrix_ulong * 
gsl_matrix_ulong_alloc_from_matrix (gsl_matrix_ulong * m,
                                    size_t k1, size_t k2,
                                    size_t n1, size_t n2);

gsl_vector_ulong * 
gsl_vector_ulong_alloc_row_from_matrix (gsl_matrix_ulong * m,
                                        size_t i);

gsl_vector_ulong * 
gsl_vector_ulong_alloc_col_from_matrix (gsl_matrix_ulong * m,
                                        size_t j);

void gsl_matrix_ulong_free (gsl_matrix_ulong * m);

unsigned long * gsl_matrix_ulong_ptr(const gsl_matrix_ulong * m, size_t i, size_t j);
unsigned long   gsl_matrix_ulong_get(const gsl_matrix_ulong * m, size_t i, size_t j);
void    gsl_matrix_ulong_set(gsl_matrix_ulong * m, size_t i, size_t j, unsigned long x);

int gsl_matrix_ulong_fread (FILE * stream, gsl_matrix_ulong * m) ;
int gsl_matrix_ulong_fwrite (FILE * stream, const gsl_matrix_ulong * m) ;
int gsl_matrix_ulong_fscanf (FILE * stream, gsl_matrix_ulong * m);
int gsl_matrix_ulong_fprintf (FILE * stream, const gsl_matrix_ulong * m, const char * format);
 
int gsl_matrix_ulong_copy_row(gsl_vector_ulong * v, const gsl_matrix_ulong * m, size_t i);
int gsl_matrix_ulong_copy_col(gsl_vector_ulong * v, const gsl_matrix_ulong * m, size_t j);
int gsl_matrix_ulong_set_row(gsl_matrix_ulong * m, size_t i, const gsl_vector_ulong * v);
int gsl_matrix_ulong_set_col(gsl_matrix_ulong * m, size_t j, const gsl_vector_ulong * v);

int gsl_matrix_ulong_swap_rows(gsl_matrix_ulong * m, size_t i, size_t j);
int gsl_matrix_ulong_swap_cols(gsl_matrix_ulong * m, size_t i, size_t j);
int gsl_matrix_ulong_swap_rowcol(gsl_matrix_ulong * m, size_t i, size_t j);

int gsl_matrix_ulong_copy(gsl_matrix_ulong * dest, const gsl_matrix_ulong * src);

int gsl_vector_ulong_view_row_from_matrix (gsl_vector_ulong * v, gsl_matrix_ulong * m, size_t i);
int gsl_vector_ulong_view_col_from_matrix (gsl_vector_ulong * v, gsl_matrix_ulong * m, size_t j);

int gsl_matrix_ulong_view_from_vector (gsl_matrix_ulong * m, 
                                       gsl_vector_ulong * base,
                                       size_t offset, 
                                       size_t n1, size_t n2, size_t d2);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
unsigned long
gsl_matrix_ulong_get(const gsl_matrix_ulong * m, 
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
gsl_matrix_ulong_set(gsl_matrix_ulong * m, 
		     const size_t i, const size_t j, const unsigned long x)
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

#endif /* GSL_MATRIX_ULONG_H */
