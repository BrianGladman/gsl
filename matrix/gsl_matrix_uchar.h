#ifndef __GSL_MATRIX_UCHAR_H__
#define __GSL_MATRIX_UCHAR_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector_uchar.h>

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

typedef struct gsl_matrix_uchar_struct gsl_matrix_uchar;

struct gsl_matrix_uchar_struct
{
  size_t size1;
  size_t size2;
  size_t tda;
  unsigned char * data;
  gsl_block_uchar * block;
} ;

gsl_matrix_uchar * 
gsl_matrix_uchar_alloc (const size_t n1, const size_t n2);

gsl_matrix_uchar * 
gsl_matrix_uchar_calloc (const size_t n1, const size_t n2);

gsl_matrix_uchar * 
gsl_matrix_uchar_alloc_from_block (gsl_block_uchar * b, 
                                   const size_t offset, 
                                   const size_t n1, 
                                   const size_t n2, 
                                   const size_t d2);

gsl_matrix_uchar * 
gsl_matrix_uchar_alloc_from_matrix (gsl_matrix_uchar * m,
                                    const size_t k1, 
                                    const size_t k2,
                                    const size_t n1, 
                                    const size_t n2);

gsl_vector_uchar * 
gsl_vector_uchar_alloc_row_from_matrix (gsl_matrix_uchar * m,
                                        const size_t i);

gsl_vector_uchar * 
gsl_vector_uchar_alloc_col_from_matrix (gsl_matrix_uchar * m,
                                        const size_t j);

void gsl_matrix_uchar_free (gsl_matrix_uchar * m);

void gsl_matrix_uchar_set_zero (gsl_matrix_uchar * m);
void gsl_matrix_uchar_set_identity (gsl_matrix_uchar * m);
void gsl_matrix_uchar_set_all (gsl_matrix_uchar * m, unsigned char x);

unsigned char * gsl_matrix_uchar_ptr(const gsl_matrix_uchar * m, const size_t i, const size_t j);
unsigned char   gsl_matrix_uchar_get(const gsl_matrix_uchar * m, const size_t i, const size_t j);
void    gsl_matrix_uchar_set(gsl_matrix_uchar * m, const size_t i, const size_t j, unsigned char x);

int gsl_matrix_uchar_fread (FILE * stream, gsl_matrix_uchar * m) ;
int gsl_matrix_uchar_fwrite (FILE * stream, const gsl_matrix_uchar * m) ;
int gsl_matrix_uchar_fscanf (FILE * stream, gsl_matrix_uchar * m);
int gsl_matrix_uchar_fprintf (FILE * stream, const gsl_matrix_uchar * m, const char * format);
 
int gsl_matrix_uchar_memcpy(gsl_matrix_uchar * dest, const gsl_matrix_uchar * src);

int gsl_matrix_uchar_swap_rows(gsl_matrix_uchar * m, const size_t i, const size_t j);
int gsl_matrix_uchar_swap_columns(gsl_matrix_uchar * m, const size_t i, const size_t j);
int gsl_matrix_uchar_swap_rowcol(gsl_matrix_uchar * m, const size_t i, const size_t j);

gsl_matrix_uchar gsl_matrix_uchar_submatrix (gsl_matrix_uchar * m, size_t i, size_t j, size_t n1, size_t n2);
gsl_vector_uchar gsl_matrix_uchar_row (gsl_matrix_uchar * m, size_t i);
gsl_vector_uchar gsl_matrix_uchar_column (gsl_matrix_uchar * m, size_t j);
gsl_vector_uchar gsl_matrix_uchar_diagonal (gsl_matrix_uchar * m);

/***********************************************************************/
/* The functions below are obsolete                                    */
/***********************************************************************/
int gsl_matrix_uchar_get_row(gsl_vector_uchar * v, const gsl_matrix_uchar * m, const size_t i);
int gsl_matrix_uchar_get_col(gsl_vector_uchar * v, const gsl_matrix_uchar * m, const size_t j);
int gsl_matrix_uchar_set_row(gsl_matrix_uchar * m, const size_t i, const gsl_vector_uchar * v);
int gsl_matrix_uchar_set_col(gsl_matrix_uchar * m, const size_t j, const gsl_vector_uchar * v);

int gsl_vector_uchar_view_row_from_matrix (gsl_vector_uchar * v, gsl_matrix_uchar * m, const size_t i);
int gsl_vector_uchar_view_col_from_matrix (gsl_vector_uchar * v, gsl_matrix_uchar * m, const size_t j);

int gsl_matrix_uchar_view_from_vector (gsl_matrix_uchar * m, gsl_vector_uchar * base, const size_t offset, const size_t n1, const size_t n2, const size_t d2);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
unsigned char
gsl_matrix_uchar_get(const gsl_matrix_uchar * m, const size_t i, const size_t j)
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
gsl_matrix_uchar_set(gsl_matrix_uchar * m, const size_t i, const size_t j, const unsigned char x)
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

#endif /* __GSL_MATRIX_UCHAR_H__ */
