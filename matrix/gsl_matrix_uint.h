#ifndef __GSL_MATRIX_UINT_H__
#define __GSL_MATRIX_UINT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector_uint.h>

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

typedef struct gsl_matrix_uint_struct gsl_matrix_uint;

struct gsl_matrix_uint_struct
{
  size_t size1;
  size_t size2;
  size_t tda;
  unsigned int * data;
  gsl_block_uint * block;
} ;

gsl_matrix_uint * 
gsl_matrix_uint_alloc (const size_t n1, const size_t n2);

gsl_matrix_uint * 
gsl_matrix_uint_calloc (const size_t n1, const size_t n2);

gsl_matrix_uint * 
gsl_matrix_uint_alloc_from_block (gsl_block_uint * b, 
                                   const size_t offset, 
                                   const size_t n1, 
                                   const size_t n2, 
                                   const size_t d2);

gsl_matrix_uint * 
gsl_matrix_uint_alloc_from_matrix (gsl_matrix_uint * m,
                                    const size_t k1, 
                                    const size_t k2,
                                    const size_t n1, 
                                    const size_t n2);

gsl_vector_uint * 
gsl_vector_uint_alloc_row_from_matrix (gsl_matrix_uint * m,
                                        const size_t i);

gsl_vector_uint * 
gsl_vector_uint_alloc_col_from_matrix (gsl_matrix_uint * m,
                                        const size_t j);

void gsl_matrix_uint_free (gsl_matrix_uint * m);

void gsl_matrix_uint_set_zero (gsl_matrix_uint * m);
void gsl_matrix_uint_set_identity (gsl_matrix_uint * m);
void gsl_matrix_uint_set_all (gsl_matrix_uint * m, unsigned int x);

unsigned int * gsl_matrix_uint_ptr(const gsl_matrix_uint * m, const size_t i, const size_t j);
unsigned int   gsl_matrix_uint_get(const gsl_matrix_uint * m, const size_t i, const size_t j);
void    gsl_matrix_uint_set(gsl_matrix_uint * m, const size_t i, const size_t j, unsigned int x);

int gsl_matrix_uint_fread (FILE * stream, gsl_matrix_uint * m) ;
int gsl_matrix_uint_fwrite (FILE * stream, const gsl_matrix_uint * m) ;
int gsl_matrix_uint_fscanf (FILE * stream, gsl_matrix_uint * m);
int gsl_matrix_uint_fprintf (FILE * stream, const gsl_matrix_uint * m, const char * format);
 
int gsl_matrix_uint_get_row(gsl_vector_uint * v, const gsl_matrix_uint * m, const size_t i);
int gsl_matrix_uint_get_col(gsl_vector_uint * v, const gsl_matrix_uint * m, const size_t j);
int gsl_matrix_uint_set_row(gsl_matrix_uint * m, const size_t i, const gsl_vector_uint * v);
int gsl_matrix_uint_set_col(gsl_matrix_uint * m, const size_t j, const gsl_vector_uint * v);

int gsl_matrix_uint_swap_rows(gsl_matrix_uint * m, const size_t i, const size_t j);
int gsl_matrix_uint_swap_cols(gsl_matrix_uint * m, const size_t i, const size_t j);
int gsl_matrix_uint_swap_rowcol(gsl_matrix_uint * m, const size_t i, const size_t j);

int gsl_matrix_uint_memcpy(gsl_matrix_uint * dest, const gsl_matrix_uint * src);

int gsl_vector_uint_view_row_from_matrix (gsl_vector_uint * v, gsl_matrix_uint * m, const size_t i);
int gsl_vector_uint_view_col_from_matrix (gsl_vector_uint * v, gsl_matrix_uint * m, const size_t j);

int gsl_matrix_uint_view_from_vector (gsl_matrix_uint * m, 
                                       gsl_vector_uint * base,
                                       const size_t offset, 
                                       const size_t n1, const size_t n2, const size_t d2);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
unsigned int
gsl_matrix_uint_get(const gsl_matrix_uint * m, 
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
  return m->data[i * m->tda + j] ;
} 

extern inline 
void
gsl_matrix_uint_set(gsl_matrix_uint * m, 
		     const size_t i, const size_t j, const unsigned int x)
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

#endif /* __GSL_MATRIX_UINT_H__ */
