#ifndef GSL_MATRIX_COMPLEX_FLOAT_H
#define GSL_MATRIX_COMPLEX_FLOAT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>
#include <gsl_complex.h>

typedef struct
{
  size_t size1;
  size_t size2;
  float * data;
} gsl_matrix_complex_float ;

gsl_matrix_complex_float * gsl_matrix_complex_float_float_alloc (size_t n1, size_t n2);
gsl_matrix_complex_float * gsl_matrix_complex_float_calloc (size_t n1, size_t n2);
void                             gsl_matrix_complex_float_free (gsl_matrix_complex_float * m);

gsl_complex_float * gsl_matrix_complex_float_ptr(const gsl_matrix_complex_float * m, size_t i, size_t j);
gsl_complex_float   gsl_matrix_complex_float_get(const gsl_matrix_complex_float * m, size_t i, size_t j);
void                gsl_matrix_complex_float_set(gsl_matrix_complex_float * m, size_t i, size_t j, gsl_complex_float x);

int gsl_matrix_complex_float_fread (FILE * stream, gsl_matrix_complex_float * m) ;
int gsl_matrix_complex_float_fwrite (FILE * stream, const gsl_matrix_complex_float * m) ;
int gsl_matrix_complex_float_fscanf (FILE * stream, gsl_matrix_complex_float * m);
int gsl_matrix_complex_float_fprintf (FILE * stream, const gsl_matrix_complex_float * m, const char * format);

extern int gsl_check_range ;



#ifdef HAVE_INLINE

extern inline 
gsl_complex_float
gsl_matrix_complex_float_get(const gsl_matrix_complex_float * m, 
		     const size_t i, const size_t j)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_CONTINUE("first index out of range", GSL_EINVAL) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_CONTINUE("second index out of range", GSL_EINVAL) ;
    }
#endif
  return *(gsl_complex_float *)(m->data + i * m->size2 + j) ;
} 

extern inline 
void
gsl_matrix_complex_float_set(gsl_matrix_complex_float * m, 
		     const size_t i, const size_t j, const gsl_complex_float x)
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
  *(gsl_complex_float *)(m->data + i * m->size2 + j) = x ;
}
#endif /* HAVE_INLINE */

#endif /* !GSL_MATRIX_COMPLEX_FLOAT_H */
