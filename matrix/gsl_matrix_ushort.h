#ifndef __GSL_MATRIX_USHORT_H__
#define __GSL_MATRIX_USHORT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector_ushort.h>

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

typedef struct gsl_matrix_ushort_struct gsl_matrix_ushort;

struct gsl_matrix_ushort_struct
{
  size_t size1;
  size_t size2;
  size_t tda;
  unsigned short * data;
  gsl_block_ushort * block;
} ;

gsl_matrix_ushort * 
gsl_matrix_ushort_alloc (const size_t n1, const size_t n2);

gsl_matrix_ushort * 
gsl_matrix_ushort_calloc (const size_t n1, const size_t n2);

gsl_matrix_ushort * 
gsl_matrix_ushort_alloc_from_block (gsl_block_ushort * b, 
                                   const size_t offset, 
                                   const size_t n1, 
                                   const size_t n2, 
                                   const size_t d2);

gsl_matrix_ushort * 
gsl_matrix_ushort_alloc_from_matrix (gsl_matrix_ushort * m,
                                    const size_t k1, 
                                    const size_t k2,
                                    const size_t n1, 
                                    const size_t n2);

gsl_vector_ushort * 
gsl_vector_ushort_alloc_row_from_matrix (gsl_matrix_ushort * m,
                                        const size_t i);

gsl_vector_ushort * 
gsl_vector_ushort_alloc_col_from_matrix (gsl_matrix_ushort * m,
                                        const size_t j);

void gsl_matrix_ushort_free (gsl_matrix_ushort * m);

void gsl_matrix_ushort_set_zero (gsl_matrix_ushort * m);
void gsl_matrix_ushort_set_identity (gsl_matrix_ushort * m);
void gsl_matrix_ushort_set_all (gsl_matrix_ushort * m, unsigned short x);

unsigned short * gsl_matrix_ushort_ptr(const gsl_matrix_ushort * m, const size_t i, const size_t j);
unsigned short   gsl_matrix_ushort_get(const gsl_matrix_ushort * m, const size_t i, const size_t j);
void    gsl_matrix_ushort_set(gsl_matrix_ushort * m, const size_t i, const size_t j, unsigned short x);

int gsl_matrix_ushort_fread (FILE * stream, gsl_matrix_ushort * m) ;
int gsl_matrix_ushort_fwrite (FILE * stream, const gsl_matrix_ushort * m) ;
int gsl_matrix_ushort_fscanf (FILE * stream, gsl_matrix_ushort * m);
int gsl_matrix_ushort_fprintf (FILE * stream, const gsl_matrix_ushort * m, const char * format);
 
int gsl_matrix_ushort_memcpy(gsl_matrix_ushort * dest, const gsl_matrix_ushort * src);

int gsl_matrix_ushort_swap_rows(gsl_matrix_ushort * m, const size_t i, const size_t j);
int gsl_matrix_ushort_swap_columns(gsl_matrix_ushort * m, const size_t i, const size_t j);
int gsl_matrix_ushort_swap_rowcol(gsl_matrix_ushort * m, const size_t i, const size_t j);

gsl_matrix_ushort gsl_matrix_ushort_submatrix (gsl_matrix_ushort * m, size_t i, size_t j, size_t n1, size_t n2);
gsl_vector_ushort gsl_matrix_ushort_row (gsl_matrix_ushort * m, size_t i);
gsl_vector_ushort gsl_matrix_ushort_column (gsl_matrix_ushort * m, size_t j);
gsl_vector_ushort gsl_matrix_ushort_diagonal (gsl_matrix_ushort * m);

unsigned short gsl_matrix_ushort_max (const gsl_matrix_ushort * m);
unsigned short gsl_matrix_ushort_min (const gsl_matrix_ushort * m);
void gsl_matrix_ushort_minmax (const gsl_matrix_ushort * m,
                              unsigned short * min_out, size_t * imin, size_t * jmin,
                              unsigned short * max_out, size_t * imax, size_t * jmax);

/***********************************************************************/
/* The functions below are obsolete                                    */
/***********************************************************************/
int gsl_matrix_ushort_get_row(gsl_vector_ushort * v, const gsl_matrix_ushort * m, const size_t i);
int gsl_matrix_ushort_get_col(gsl_vector_ushort * v, const gsl_matrix_ushort * m, const size_t j);
int gsl_matrix_ushort_set_row(gsl_matrix_ushort * m, const size_t i, const gsl_vector_ushort * v);
int gsl_matrix_ushort_set_col(gsl_matrix_ushort * m, const size_t j, const gsl_vector_ushort * v);

int gsl_vector_ushort_view_row_from_matrix (gsl_vector_ushort * v, gsl_matrix_ushort * m, const size_t i);
int gsl_vector_ushort_view_col_from_matrix (gsl_vector_ushort * v, gsl_matrix_ushort * m, const size_t j);

int gsl_matrix_ushort_view_from_vector (gsl_matrix_ushort * m, gsl_vector_ushort * base, const size_t offset, const size_t n1, const size_t n2, const size_t d2);



extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
unsigned short
gsl_matrix_ushort_get(const gsl_matrix_ushort * m, const size_t i, const size_t j)
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
gsl_matrix_ushort_set(gsl_matrix_ushort * m, const size_t i, const size_t j, const unsigned short x)
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

#endif /* __GSL_MATRIX_USHORT_H__ */
