#ifndef GSL_VECTOR_COMPLEX_FLOAT_H 
#define GSL_VECTOR_COMPLEX_FLOAT_H 

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>
#include <gsl_complex.h>

typedef struct
{
  size_t size;
  size_t stride;
  float * data;
} 
gsl_vector_complex_float ;


gsl_vector_complex_float * gsl_vector_complex_float_alloc (size_t n);
gsl_vector_complex_float * gsl_vector_complex_float_calloc (size_t n);
void gsl_vector_complex_float_free (gsl_vector_complex_float * v);

gsl_complex_float * gsl_vector_complex_float_ptr(const gsl_vector_complex_float * v, size_t i);
gsl_complex_float   gsl_vector_complex_float_get(const gsl_vector_complex_float * v, size_t i);
void                gsl_vector_complex_float_set(gsl_vector_complex_float * v, size_t i, gsl_complex_float z);

int gsl_vector_complex_float_fread (FILE * stream, gsl_vector_complex_float * v) ;
int gsl_vector_complex_float_fwrite (FILE * stream, const gsl_vector_complex_float * v) ;
int gsl_vector_complex_float_fscanf (FILE * stream, gsl_vector_complex_float * v);
int gsl_vector_complex_float_fprintf (FILE * stream, const gsl_vector_complex_float * v,
			              const char * format);

int gsl_block_complex_float_fread (FILE * stream, float * data, size_t n) ;
int gsl_block_complex_float_fwrite (FILE * stream, const float * data, size_t n) ;
int gsl_block_complex_float_fscanf (FILE * stream, float * data, size_t n);
int gsl_block_complex_float_fprintf (FILE * stream, const float * data, size_t n,
			             const char * format);

extern int gsl_check_range ;

#ifndef  GSL_VECTOR_COMPLEX_REAL
#define  GSL_VECTOR_COMPLEX_REAL(z, i)  ((z)->data[2*(i)])
#define  GSL_VECTOR_COMPLEX_IMAG(z, i)  ((z)->data[2*(i) + 1])
#endif

#define GSL_COMPLEX_FLOAT_AT(zv, i)  ((gsl_complex_float *)  &((zv)->data[2*(i)]))


/* inline functions if you are using GCC or otherwise enlightened cc */

#ifdef HAVE_INLINE
extern inline 
gsl_complex_float *
gsl_vector_complex_float_ptr(const gsl_vector_complex_float * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0);
    }
#endif
  return GSL_COMPLEX_FLOAT_AT(v, i);
}

extern inline 
gsl_complex_float
gsl_vector_complex_float_get(const gsl_vector_complex_float * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_CONTINUE("index out of range", GSL_EINVAL);
    }
#endif
  return *GSL_COMPLEX_FLOAT_AT(v, i);
}

extern inline 
void
gsl_vector_complex_float_set(gsl_vector_complex_float * v, const size_t i, gsl_complex_float z)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING("index out of range", GSL_EINVAL);
    }
#endif
  *GSL_COMPLEX_FLOAT_AT(v, i) = z;
}

#endif /* HAVE_INLINE */


#endif /* !GSL_VECTOR_COMPLEX_FLOAT_H */
