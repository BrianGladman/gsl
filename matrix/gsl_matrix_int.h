#ifndef __GSL_MATRIX_INT_H__
#define __GSL_MATRIX_INT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector_int.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct gsl_matrix_int_struct gsl_matrix_int;

struct gsl_matrix_int_struct
{
  size_t size1;
  size_t size2;
  size_t tda;
  int * data;
  gsl_block_int * block;
} ;

gsl_matrix_int * 
gsl_matrix_int_alloc (const size_t n1, const size_t n2);

gsl_matrix_int * 
gsl_matrix_int_calloc (const size_t n1, const size_t n2);

gsl_matrix_int * 
gsl_matrix_int_alloc_from_block (gsl_block_int * b, 
                                   const size_t offset, 
                                   const size_t n1, 
                                   const size_t n2, 
                                   const size_t d2);

gsl_matrix_int * 
gsl_matrix_int_alloc_from_matrix (gsl_matrix_int * m,
                                    const size_t k1, 
                                    const size_t k2,
                                    const size_t n1, 
                                    const size_t n2);

gsl_vector_int * 
gsl_vector_int_alloc_row_from_matrix (gsl_matrix_int * m,
                                        const size_t i);

gsl_vector_int * 
gsl_vector_int_alloc_col_from_matrix (gsl_matrix_int * m,
                                        const size_t j);

void gsl_matrix_int_free (gsl_matrix_int * m);

void gsl_matrix_int_set_zero (gsl_matrix_int * m);
void gsl_matrix_int_set_identity (gsl_matrix_int * m);
void gsl_matrix_int_set_all (gsl_matrix_int * m, int x);

int * gsl_matrix_int_ptr(const gsl_matrix_int * m, const size_t i, const size_t j);
int   gsl_matrix_int_get(const gsl_matrix_int * m, const size_t i, const size_t j);
void    gsl_matrix_int_set(gsl_matrix_int * m, const size_t i, const size_t j, int x);

int gsl_matrix_int_fread (FILE * stream, gsl_matrix_int * m) ;
int gsl_matrix_int_fwrite (FILE * stream, const gsl_matrix_int * m) ;
int gsl_matrix_int_fscanf (FILE * stream, gsl_matrix_int * m);
int gsl_matrix_int_fprintf (FILE * stream, const gsl_matrix_int * m, const char * format);
 
int gsl_matrix_int_memcpy(gsl_matrix_int * dest, const gsl_matrix_int * src);
int gsl_matrix_int_swap(gsl_matrix_int * m1, const gsl_matrix_int * m2);

int gsl_matrix_int_swap_rows(gsl_matrix_int * m, const size_t i, const size_t j);
int gsl_matrix_int_swap_columns(gsl_matrix_int * m, const size_t i, const size_t j);
int gsl_matrix_int_swap_rowcol(gsl_matrix_int * m, const size_t i, const size_t j);
int gsl_matrix_int_transpose (gsl_matrix_int * m);


gsl_matrix_int gsl_matrix_int_submatrix (gsl_matrix_int * m, size_t i, size_t j, size_t n1, size_t n2);
gsl_vector_int gsl_matrix_int_row (gsl_matrix_int * m, size_t i);
gsl_vector_int gsl_matrix_int_column (gsl_matrix_int * m, size_t j);
gsl_vector_int gsl_matrix_int_diagonal (gsl_matrix_int * m);

int gsl_matrix_int_max (const gsl_matrix_int * m);
int gsl_matrix_int_min (const gsl_matrix_int * m);
void gsl_matrix_int_minmax (const gsl_matrix_int * m, int * min_out, int * max_out);

void gsl_matrix_int_max_index (const gsl_matrix_int * m, size_t * imax, size_t *jmax);
void gsl_matrix_int_min_index (const gsl_matrix_int * m, size_t * imin, size_t *jmin);
void gsl_matrix_int_minmax_index (const gsl_matrix_int * m, size_t * imin, size_t * jmin, size_t * imax, size_t * jmax);

int gsl_matrix_int_isnull (const gsl_matrix_int * m);

/***********************************************************************/
/* The functions below are obsolete                                    */
/***********************************************************************/
int gsl_matrix_int_get_row(gsl_vector_int * v, const gsl_matrix_int * m, const size_t i);
int gsl_matrix_int_get_col(gsl_vector_int * v, const gsl_matrix_int * m, const size_t j);
int gsl_matrix_int_set_row(gsl_matrix_int * m, const size_t i, const gsl_vector_int * v);
int gsl_matrix_int_set_col(gsl_matrix_int * m, const size_t j, const gsl_vector_int * v);

int gsl_vector_int_view_row_from_matrix (gsl_vector_int * v, gsl_matrix_int * m, const size_t i);
int gsl_vector_int_view_col_from_matrix (gsl_vector_int * v, gsl_matrix_int * m, const size_t j);

int gsl_matrix_int_view_from_vector (gsl_matrix_int * m, gsl_vector_int * base, const size_t offset, const size_t n1, const size_t n2, const size_t d2);



extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
int
gsl_matrix_int_get(const gsl_matrix_int * m, const size_t i, const size_t j)
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
  return m->data[i * m->tda + j] ;
} 

extern inline 
void
gsl_matrix_int_set(gsl_matrix_int * m, const size_t i, const size_t j, const int x)
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
  m->data[i * m->tda + j] = x ;
}
#endif

__END_DECLS

#endif /* __GSL_MATRIX_INT_H__ */
