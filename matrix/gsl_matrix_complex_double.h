#ifndef GSL_MATRIX_COMPLEX_DOUBLE_H
#define GSL_MATRIX_COMPLEX_DOUBLE_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>
#include <gsl_complex.h>

typedef struct
{
  size_t size1;
  size_t size2;
  double * data;
  
} gsl_matrix_complex ;

gsl_matrix_complex * gsl_matrix_complex_alloc (size_t n1, size_t n2);
gsl_matrix_complex * gsl_matrix_complex_calloc (size_t n1, size_t n2);
void                 gsl_matrix_complex_free (gsl_matrix_complex * m);

gsl_complex * gsl_matrix_complex_ptr(const gsl_matrix_complex * m, size_t i, size_t j);
gsl_complex   gsl_matrix_complex_get(const gsl_matrix_complex * m, size_t i, size_t j);
void          gsl_matrix_complex_set(gsl_matrix_complex * m, size_t i, size_t j, gsl_complex x);

int gsl_matrix_complex_fread (FILE * stream, gsl_matrix_complex * m) ;
int gsl_matrix_complex_fwrite (FILE * stream, const gsl_matrix_complex * m) ;
int gsl_matrix_complex_fscanf (FILE * stream, gsl_matrix_complex * m);
int gsl_matrix_complex_fprintf (FILE * stream, const gsl_matrix_complex * m, const char * format);

extern int gsl_check_range ;



#ifdef HAVE_INLINE

extern inline 
gsl_complex
gsl_matrix_complex_get(const gsl_matrix_complex * m, 
		     const size_t i, const size_t j)
{
#ifndef GSL_RANGE_CHECK_OFF
  static const gsl_complex gsl_complex_zero = {{0, 0}} ;
  if (i >= m->size1)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("first index out of range", GSL_EINVAL, gsl_complex_zero) ;
    }
  else if (j >= m->size2) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("second index out of range", GSL_EINVAL, gsl_complex_zero) ;
    }
#endif
  return *(gsl_complex *)(m->data + 2*(i * m->size2 + j)) ;
} 

extern inline 
void
gsl_matrix_complex_set(gsl_matrix_complex * m, 
		       const size_t i, const size_t j,
		       const gsl_complex x)
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
  *(gsl_complex *)(m->data + 2*(i * m->size2 + j)) = x ;
}
#endif

#endif /* !GSL_MATRIX_COMPLEX_DOUBLE_H */
