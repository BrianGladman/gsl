#ifndef GSL_VECTOR_FLOAT_H
#define GSL_VECTOR_FLOAT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>

typedef struct
  {
    size_t size;
    size_t stride;
    float *data;
  }
gsl_vector_float;

gsl_vector_float *gsl_vector_float_alloc (size_t n);
gsl_vector_float *gsl_vector_float_calloc (size_t n);
void gsl_vector_float_free (gsl_vector_float * v);

float *gsl_vector_float_ptr (const gsl_vector_float * v, const size_t i);
float gsl_vector_float_get (const gsl_vector_float * v, const size_t i);
void gsl_vector_float_set (gsl_vector_float * v, const size_t i, float x);

int gsl_vector_float_fread (FILE * stream, gsl_vector_float * v);
int gsl_vector_float_fwrite (FILE * stream, const gsl_vector_float * v);
int gsl_vector_float_fscanf (FILE * stream, gsl_vector_float * v);
int gsl_vector_float_fprintf (FILE * stream, const gsl_vector_float * v,
			      const char *format);

int gsl_block_float_fread (FILE * stream, float *data, size_t n, size_t stride);
int gsl_block_float_fwrite (FILE * stream, const float *data, size_t n, size_t stride);
int gsl_block_float_fscanf (FILE * stream, float *data, size_t n, size_t stride);
int gsl_block_float_fprintf (FILE * stream, const float *data, size_t n, size_t stride,
			     const char *format);

extern int gsl_check_range;



#ifdef HAVE_INLINE

extern inline
float *
gsl_vector_float_ptr (const gsl_vector_float * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)		/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return (float *) (v->data + i);
}

extern inline
float
gsl_vector_float_get (const gsl_vector_float * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)		/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return v->data[i];
}

extern inline
void
gsl_vector_float_set (gsl_vector_float * v, const size_t i, float x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)		/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING ("index out of range", GSL_EINVAL);
    }
#endif
  v->data[i] = x;
}

#endif /* HAVE_INLINE */

#endif /* !GSL_VECTOR_FLOAT_H */
