#ifndef GSL_MATRIX_COMPLEX_LONG_DOUBLE_H
#define GSL_MATRIX_COMPLEX_LONG_DOUBLE_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>
#include <gsl_complex.h>
#include <gsl_vector_complex_long_double.h>

typedef struct
{
  size_t size1;
  size_t size2;
  long double * data;
  
} gsl_matrix_complex_long_double ;

gsl_matrix_complex_long_double * gsl_matrix_complex_long_double_alloc (size_t n1, size_t n2);
gsl_matrix_complex_long_double * gsl_matrix_complex_long_double_calloc (size_t n1, size_t n2);
void                 gsl_matrix_complex_long_double_free (gsl_matrix_complex_long_double * m);

gsl_complex_long_double * gsl_matrix_complex_long_double_ptr(const gsl_matrix_complex_long_double * m, size_t i, size_t j);
gsl_complex_long_double   gsl_matrix_complex_long_double_get(const gsl_matrix_complex_long_double * m, size_t i, size_t j);
void                      gsl_matrix_complex_long_double_set(gsl_matrix_complex_long_double * m, size_t i, size_t j, gsl_complex_long_double x);

int gsl_matrix_complex_long_double_fread (FILE * stream, gsl_matrix_complex_long_double * m) ;
int gsl_matrix_complex_long_double_fwrite (FILE * stream, const gsl_matrix_complex_long_double * m) ;
int gsl_matrix_complex_long_double_fscanf (FILE * stream, gsl_matrix_complex_long_double * m);
int gsl_matrix_complex_long_double_fprintf (FILE * stream, const gsl_matrix_complex_long_double * m, const char * format);

int gsl_matrix_complex_long_double_copy_row(const gsl_matrix_complex_long_double * m, size_t i, gsl_vector_complex_long_double * v);
int gsl_matrix_complex_long_double_copy_col(const gsl_matrix_complex_long_double * m, size_t j, gsl_vector_complex_long_double * v);
int gsl_matrix_complex_long_double_set_row(gsl_matrix_complex_long_double * m, size_t i, const gsl_vector_complex_long_double * v);
int gsl_matrix_complex_long_double_set_col(gsl_matrix_complex_long_double * m, size_t j, const gsl_vector_complex_long_double * v);


extern int gsl_check_range ;


#ifdef HAVE_INLINE

extern inline 
gsl_complex_long_double
gsl_matrix_complex_long_double_get(const gsl_matrix_complex_long_double * m, 
		     const size_t i, const size_t j)
{
  const gsl_complex_long_double zero = {{0,0}};

#ifndef GSL_RANGE_CHECK_OFF
  if (i >= m->size1)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("first index out of range", GSL_EINVAL, zero) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("second index out of range", GSL_EINVAL, zero) ;
    }
#endif
  return *(gsl_complex_long_double *)(m->data + 2*(i * m->size2 + j));
} 

extern inline 
void
gsl_matrix_complex_long_double_set(gsl_matrix_complex_long_double * m, 
		                   const size_t i, const size_t j,
				   const gsl_complex_long_double x)
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
  *(gsl_complex_long_double *)(m->data + 2*(i * m->size2 + j)) = x ;
}
#endif

#endif /* !GSL_MATRIX_COMPLEX_LONG_DOUBLE_H */
