#ifndef GSL_MATRIX_LONG_H
#define GSL_MATRIX_LONG_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>
#include <gsl_vector_long.h>

typedef struct
{
  size_t size1;
  size_t size2;
  long * data;
  
} gsl_matrix_long ;

gsl_matrix_long * gsl_matrix_long_alloc (size_t n1, size_t n2);
gsl_matrix_long * gsl_matrix_long_calloc (size_t n1, size_t n2);
void gsl_matrix_long_free (gsl_matrix_long * m);

long * gsl_matrix_long_ptr(const gsl_matrix_long * m, size_t i, size_t j);
long   gsl_matrix_long_get(const gsl_matrix_long * m, size_t i, size_t j);
void  gsl_matrix_long_set(gsl_matrix_long * m, size_t i,  size_t j, long x);

int gsl_matrix_long_fread (FILE * stream, gsl_matrix_long * m) ;
int gsl_matrix_long_fwrite (FILE * stream, const gsl_matrix_long * m) ;
int gsl_matrix_long_fscanf (FILE * stream, gsl_matrix_long * m);
int gsl_matrix_long_fprintf (FILE * stream, const gsl_matrix_long * m, const char * format);

int gsl_matrix_long_copy_row(const gsl_matrix_long * m, size_t i, gsl_vector_long * v);
int gsl_matrix_long_copy_col(const gsl_matrix_long * m, size_t j, gsl_vector_long * v);
int gsl_matrix_long_set_row(gsl_matrix_long * m, size_t i, const gsl_vector_long * v);
int gsl_matrix_long_set_col(gsl_matrix_long * m, size_t j, const gsl_vector_long * v);

extern int gsl_check_range ;


#ifdef HAVE_INLINE

extern inline 
long
gsl_matrix_long_get(const gsl_matrix_long * m, const size_t i, const size_t j)
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
gsl_matrix_long_set(gsl_matrix_long * m, 
		   const size_t i, const size_t j, 
		   const long x)
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

#endif /* !GSL_MATRIX_LONG_H */
