#ifndef GSL_VECTOR_COMPLEX_FLOAT_H
#define GSL_VECTOR_COMPLEX_FLOAT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_vector_complex.h>
#include <gsl_block_complex_float.h>

struct gsl_vector_complex_float_struct
{
  size_t size;
  size_t stride;
  float *data;
  gsl_block_complex_float *block;
};

typedef struct gsl_vector_complex_float_struct gsl_vector_complex_float;

gsl_vector_complex_float *gsl_vector_complex_float_alloc (size_t n);
gsl_vector_complex_float *gsl_vector_complex_float_calloc (size_t n);

gsl_vector_complex_float 
*gsl_vector_complex_float_alloc_from_block (gsl_block_complex_float * b, 
                                            size_t offset, 
                                            size_t n, 
                                            size_t stride);

gsl_vector_complex_float 
*gsl_vector_complex_float_alloc_from_vector (gsl_vector_complex_float * v, 
                                             size_t offset, 
                                             size_t n, 
                                             size_t stride);

void gsl_vector_complex_float_free (gsl_vector_complex_float * v);

gsl_complex_float 
*gsl_vector_complex_float_ptr (const gsl_vector_complex_float * v, size_t i);

gsl_complex_float 
gsl_vector_complex_float_get (const gsl_vector_complex_float * v, size_t i);

void gsl_vector_complex_float_set (gsl_vector_complex_float * v, size_t i,
                                   gsl_complex_float z);

int gsl_vector_complex_float_fread (FILE * stream,
				    gsl_vector_complex_float * v);
int gsl_vector_complex_float_fwrite (FILE * stream,
				     const gsl_vector_complex_float * v);
int gsl_vector_complex_float_fscanf (FILE * stream,
				     gsl_vector_complex_float * v);
int gsl_vector_complex_float_fprintf (FILE * stream,
				      const gsl_vector_complex_float * v,
				      const char *format);

int gsl_vector_complex_float_copy (gsl_vector_complex_float * dest, const gsl_vector_complex_float * src);

int gsl_vector_complex_float_swap (gsl_vector_complex_float * v, size_t i, size_t j);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
gsl_complex_float
gsl_vector_complex_float_get (const gsl_vector_complex_float * v,
			      const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)
    {
      const gsl_complex_float zero = {{0, 0}};
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, zero);
    }
#endif
  return *GSL_COMPLEX_FLOAT_AT (v, i);
}

extern inline
void
gsl_vector_complex_float_set (gsl_vector_complex_float * v,
			      const size_t i, gsl_complex_float z)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)
    {
      GSL_ERROR_RETURN_NOTHING ("index out of range", GSL_EINVAL);
    }
#endif
  *GSL_COMPLEX_FLOAT_AT (v, i) = z;
}

#endif /* HAVE_INLINE */

#endif /* GSL_VECTOR_COMPLEX_FLOAT_H */
