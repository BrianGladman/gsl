#ifndef GSL_MATRIX_UINT_H
#define GSL_MATRIX_UINT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>
#include <gsl_vector_uint.h>

typedef struct gsl_matrix_uint_struct gsl_matrix_uint;

struct gsl_matrix_uint_struct
{
  size_t size1;
  size_t size2;
  unsigned int * data;
} ;

gsl_matrix_uint * gsl_matrix_uint_alloc (size_t n1, size_t n2);
gsl_matrix_uint * gsl_matrix_uint_calloc (size_t n1, size_t n2);
void gsl_matrix_uint_free (gsl_matrix_uint * m);

unsigned int * gsl_matrix_uint_ptr(const gsl_matrix_uint * m, size_t i, size_t j);
unsigned int   gsl_matrix_uint_get(const gsl_matrix_uint * m, size_t i, size_t j);
void  gsl_matrix_uint_set(gsl_matrix_uint * m, size_t i,  size_t j, unsigned int x);

int gsl_matrix_uint_fread (FILE * stream, gsl_matrix_uint * m) ;
int gsl_matrix_uint_fwrite (FILE * stream, const gsl_matrix_uint * m) ;
int gsl_matrix_uint_fscanf (FILE * stream, gsl_matrix_uint * m);
int gsl_matrix_uint_fprintf (FILE * stream, const gsl_matrix_uint * m, const char * format);

int gsl_matrix_uint_copy_row(const gsl_matrix_uint * m, size_t i, gsl_vector_uint * v);
int gsl_matrix_uint_copy_col(const gsl_matrix_uint * m, size_t j, gsl_vector_uint * v);
int gsl_matrix_uint_set_row(gsl_matrix_uint * m, size_t i, const gsl_vector_uint * v);
int gsl_matrix_uint_set_col(gsl_matrix_uint * m, size_t j, const gsl_vector_uint * v);

extern int gsl_check_range ;

#ifdef HAVE_INLINE

extern inline 
unsigned int
gsl_matrix_uint_get(const gsl_matrix_uint * m, const size_t i, const size_t j)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("first index out of range", GSL_EINVAL, 0) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("second index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return m->data[i * m->size2 + j] ;
} 

extern inline 
void
gsl_matrix_uint_set(gsl_matrix_uint * m, 
		   const size_t i, const size_t j, 
		   const unsigned int x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING("first index out of range", GSL_EINVAL) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING("second index out of range", GSL_EINVAL) ;
    }
#endif
  m->data[i * m->size2 + j] = x ;
}
#endif /* HAVE_INLINE */

#endif /* !GSL_MATRIX_UINT_H */
