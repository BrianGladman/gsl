#ifndef GSL_VECTOR_COMPLEX_LONG_H 
#define GSL_VECTOR_COMPLEX_LONG_H 

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>
#include <gsl_complex.h>

typedef struct
{
  size_t size;
  size_t stride;
  long double * data;
} 
gsl_vector_complex_long_double;


gsl_vector_complex_long_double * gsl_vector_complex_long_double_alloc (size_t n);
gsl_vector_complex_long_double * gsl_vector_complex_long_double_calloc (size_t n);
void gsl_vector_complex_long_double_free (gsl_vector_complex_long_double * v);

gsl_complex_long_double * gsl_vector_complex_long_double_getp(const gsl_vector_complex_long_double * v, size_t i);
void gsl_vector_complex_long_double_set(gsl_vector_complex_long_double * v, size_t i, const gsl_complex_long_double * x);

int gsl_vector_complex_long_double_fread (FILE * stream, gsl_vector_complex_long_double * v) ;
int gsl_vector_complex_long_double_fwrite (FILE * stream, const gsl_vector_complex_long_double * v) ;
int gsl_vector_complex_long_double_fscanf (FILE * stream, gsl_vector_complex_long_double * v);
int gsl_vector_complex_long_double_fprintf (FILE * stream, const gsl_vector_complex_long_double * v,
			              const char * format);

int gsl_block_complex_long_double_fread (FILE * stream, long double * data, size_t n) ;
int gsl_block_complex_long_double_fwrite (FILE * stream, const long double * data, size_t n) ;
int gsl_block_complex_long_double_fscanf (FILE * stream, long double * data, size_t n);
int gsl_block_complex_long_double_fprintf (FILE * stream, const long double * data, size_t n,
			             const char * format);

extern int gsl_check_range ;

#ifndef  GSL_VECTOR_COMPLEX_REAL
#define  GSL_VECTOR_COMPLEX_REAL(z, i)  (z->data[2*i])
#define  GSL_VECTOR_COMPLEX_IMAG(z, i)  (z->data[2*i + 1])
#endif

#define GSL_COMPLEX_LONG_DOUBLE_AT(zv, i)  ((gsl_complex_long_double *)  &(zv->data[2*i]))


/* inline functions if you are using GCC or otherwise enlightened cc */

#ifdef HAVE_INLINE
extern inline 
gsl_complex_long_double *
gsl_vector_complex_long_double_getp(const gsl_vector_complex_long_double * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0);
    }
#endif
  return GSL_COMPLEX_LONG_DOUBLE_AT(v, i);
} 

extern inline 
void
gsl_vector_complex_long_double_set(gsl_vector_complex_long_double * v, const size_t i, const gsl_complex_long_double * x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING("index out of range", GSL_EINVAL) ;
    }
#endif
  GSL_VECTOR_COMPLEX_REAL(v, i) = GSL_COMPLEX_P_REAL(x);
  GSL_VECTOR_COMPLEX_IMAG(v, i) = GSL_COMPLEX_P_IMAG(x);
}

#endif /* HAVE_INLINE */


#endif /* GSL_VECTOR_COMPLEX_LONG_H */
