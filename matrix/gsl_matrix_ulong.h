#ifndef __GSL_MATRIX_ULONG_H__
#define __GSL_MATRIX_ULONG_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector_ulong.h>

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

typedef struct gsl_matrix_ulong_struct gsl_matrix_ulong;

struct gsl_matrix_ulong_struct
{
  size_t size1;
  size_t size2;
  size_t tda;
  unsigned long * data;
  gsl_block_ulong * block;
} ;

gsl_matrix_ulong * 
gsl_matrix_ulong_alloc (const size_t n1, const size_t n2);

gsl_matrix_ulong * 
gsl_matrix_ulong_calloc (const size_t n1, const size_t n2);

gsl_matrix_ulong * 
gsl_matrix_ulong_alloc_from_block (gsl_block_ulong * b, 
                                   const size_t offset, 
                                   const size_t n1, 
                                   const size_t n2, 
                                   const size_t d2);

gsl_matrix_ulong * 
gsl_matrix_ulong_alloc_from_matrix (gsl_matrix_ulong * m,
                                    const size_t k1, 
                                    const size_t k2,
                                    const size_t n1, 
                                    const size_t n2);

gsl_vector_ulong * 
gsl_vector_ulong_alloc_row_from_matrix (gsl_matrix_ulong * m,
                                        const size_t i);

gsl_vector_ulong * 
gsl_vector_ulong_alloc_col_from_matrix (gsl_matrix_ulong * m,
                                        const size_t j);

void gsl_matrix_ulong_free (gsl_matrix_ulong * m);

void gsl_matrix_ulong_set_zero (gsl_matrix_ulong * m);
void gsl_matrix_ulong_set_identity (gsl_matrix_ulong * m);
void gsl_matrix_ulong_set_all (gsl_matrix_ulong * m, unsigned long x);

unsigned long * gsl_matrix_ulong_ptr(const gsl_matrix_ulong * m, const size_t i, const size_t j);
unsigned long   gsl_matrix_ulong_get(const gsl_matrix_ulong * m, const size_t i, const size_t j);
void    gsl_matrix_ulong_set(gsl_matrix_ulong * m, const size_t i, const size_t j, unsigned long x);

int gsl_matrix_ulong_fread (FILE * stream, gsl_matrix_ulong * m) ;
int gsl_matrix_ulong_fwrite (FILE * stream, const gsl_matrix_ulong * m) ;
int gsl_matrix_ulong_fscanf (FILE * stream, gsl_matrix_ulong * m);
int gsl_matrix_ulong_fprintf (FILE * stream, const gsl_matrix_ulong * m, const char * format);
 
int gsl_matrix_ulong_memcpy(gsl_matrix_ulong * dest, const gsl_matrix_ulong * src);

int gsl_matrix_ulong_swap_rows(gsl_matrix_ulong * m, const size_t i, const size_t j);
int gsl_matrix_ulong_swap_columns(gsl_matrix_ulong * m, const size_t i, const size_t j);
int gsl_matrix_ulong_swap_rowcol(gsl_matrix_ulong * m, const size_t i, const size_t j);

gsl_matrix_ulong gsl_matrix_ulong_submatrix (gsl_matrix_ulong * m, size_t i, size_t j, size_t n1, size_t n2);
gsl_vector_ulong gsl_matrix_ulong_row (gsl_matrix_ulong * m, size_t i);
gsl_vector_ulong gsl_matrix_ulong_column (gsl_matrix_ulong * m, size_t j);
gsl_vector_ulong gsl_matrix_ulong_diagonal (gsl_matrix_ulong * m);

unsigned long gsl_matrix_ulong_max (const gsl_matrix_ulong * m);
unsigned long gsl_matrix_ulong_min (const gsl_matrix_ulong * m);
void gsl_matrix_ulong_minmax (const gsl_matrix_ulong * m, unsigned long * min_out, unsigned long * max_out);

void gsl_matrix_ulong_max_index (const gsl_matrix_ulong * m, size_t * imax, size_t *jmax);
void gsl_matrix_ulong_min_index (const gsl_matrix_ulong * m, size_t * imin, size_t *jmin);
void gsl_matrix_ulong_minmax_index (const gsl_matrix_ulong * m, size_t * imin, size_t * jmin, size_t * imax, size_t * jmax);

/***********************************************************************/
/* The functions below are obsolete                                    */
/***********************************************************************/
int gsl_matrix_ulong_get_row(gsl_vector_ulong * v, const gsl_matrix_ulong * m, const size_t i);
int gsl_matrix_ulong_get_col(gsl_vector_ulong * v, const gsl_matrix_ulong * m, const size_t j);
int gsl_matrix_ulong_set_row(gsl_matrix_ulong * m, const size_t i, const gsl_vector_ulong * v);
int gsl_matrix_ulong_set_col(gsl_matrix_ulong * m, const size_t j, const gsl_vector_ulong * v);

int gsl_vector_ulong_view_row_from_matrix (gsl_vector_ulong * v, gsl_matrix_ulong * m, const size_t i);
int gsl_vector_ulong_view_col_from_matrix (gsl_vector_ulong * v, gsl_matrix_ulong * m, const size_t j);

int gsl_matrix_ulong_view_from_vector (gsl_matrix_ulong * m, gsl_vector_ulong * base, const size_t offset, const size_t n1, const size_t n2, const size_t d2);



extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
unsigned long
gsl_matrix_ulong_get(const gsl_matrix_ulong * m, const size_t i, const size_t j)
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
gsl_matrix_ulong_set(gsl_matrix_ulong * m, const size_t i, const size_t j, const unsigned long x)
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

#endif /* __GSL_MATRIX_ULONG_H__ */
