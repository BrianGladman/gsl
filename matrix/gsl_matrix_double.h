#ifndef GSL_MATRIX_DOUBLE_H
#define GSL_MATRIX_DOUBLE_H

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector_double.h>

typedef struct gsl_matrix_struct gsl_matrix;

struct gsl_matrix_struct
{
  size_t size1;
  size_t size2;
  size_t dim2;
  double * data;
  gsl_block * block;
} ;

gsl_matrix * 
gsl_matrix_alloc (size_t n1, size_t n2);

gsl_matrix * 
gsl_matrix_calloc (size_t n1, size_t n2);

gsl_matrix * 
gsl_matrix_alloc_from_block (gsl_block * b, size_t offset, 
                                   size_t n1, size_t n2, size_t d2);

gsl_matrix * 
gsl_matrix_alloc_from_matrix (gsl_matrix * m,
                                    size_t k1, size_t k2,
                                    size_t n1, size_t n2);

gsl_vector * 
gsl_vector_alloc_row_from_matrix (gsl_matrix * m,
                                        size_t i);

gsl_vector * 
gsl_vector_alloc_col_from_matrix (gsl_matrix * m,
                                        size_t j);

void gsl_matrix_free (gsl_matrix * m);

double * gsl_matrix_ptr(const gsl_matrix * m, size_t i, size_t j);
double   gsl_matrix_get(const gsl_matrix * m, size_t i, size_t j);
void    gsl_matrix_set(gsl_matrix * m, size_t i, size_t j, double x);

int gsl_matrix_fread (FILE * stream, gsl_matrix * m) ;
int gsl_matrix_fwrite (FILE * stream, const gsl_matrix * m) ;
int gsl_matrix_fscanf (FILE * stream, gsl_matrix * m);
int gsl_matrix_fprintf (FILE * stream, const gsl_matrix * m, const char * format);
 
int gsl_matrix_copy_row(gsl_vector * v, const gsl_matrix * m, size_t i);
int gsl_matrix_copy_col(gsl_vector * v, const gsl_matrix * m, size_t j);
int gsl_matrix_set_row(gsl_matrix * m, size_t i, const gsl_vector * v);
int gsl_matrix_set_col(gsl_matrix * m, size_t j, const gsl_vector * v);

int gsl_matrix_swap_rows(gsl_matrix * m, size_t i, size_t j);
int gsl_matrix_swap_cols(gsl_matrix * m, size_t i, size_t j);
int gsl_matrix_swap_rowcol(gsl_matrix * m, size_t i, size_t j);

int gsl_matrix_copy(gsl_matrix * dest, const gsl_matrix * src);

int gsl_vector_view_row_from_matrix (gsl_vector * v, gsl_matrix * m, size_t i);
int gsl_vector_view_col_from_matrix (gsl_vector * v, gsl_matrix * m, size_t j);

int gsl_matrix_view_from_vector (gsl_matrix * m, 
                                       gsl_vector * base,
                                       size_t offset, 
                                       size_t n1, size_t n2, size_t d2);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
double
gsl_matrix_get(const gsl_matrix * m, 
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
gsl_matrix_set(gsl_matrix * m, 
		     const size_t i, const size_t j, const double x)
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

#endif /* GSL_MATRIX_DOUBLE_H */
