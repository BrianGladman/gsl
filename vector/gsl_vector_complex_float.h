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

#ifndef GSL_VECTOR_COMPLEX
#ifdef GSL_RANGE_CHECK_OFF
#define GSL_VECTOR_COMPLEX(zv, i) (GSL_COMPLEX_AT((zv),(i)))
#else
#define GSL_VECTOR_COMPLEX(zv, i) (((i) >= (zv)->size ? (gsl_error ("index out of range", __FILE__, __LINE__, GSL_EINVAL), 0):0 , *GSL_COMPLEX_AT((zv),(i))))
#endif
#endif /* GSL_VECTOR_COMPLEX */

#define GSL_COMPLEX_FLOAT_AT(zv, i)  ((gsl_complex_float *)  &((zv)->data[2*(i)]))



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
      const gsl_complex_float zero = {{0,0}} ;
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, zero);
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
