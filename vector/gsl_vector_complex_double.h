#ifndef GSL_VECTOR_COMPLEX_DOUBLE_H
#define GSL_VECTOR_COMPLEX_DOUBLE_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_vector_complex.h>
#include <gsl_block_complex_double.h>

struct gsl_vector_complex_struct
{
  size_t size;
  size_t stride;
  double *data;
  gsl_block_complex *block;
};

typedef struct gsl_vector_complex_struct gsl_vector_complex;

gsl_vector_complex *gsl_vector_complex_alloc (size_t n);
gsl_vector_complex *gsl_vector_complex_calloc (size_t n);

gsl_vector_complex 
*gsl_vector_complex_alloc_from_block (gsl_block_complex * b, 
                                            size_t offset, 
                                            size_t n, 
                                            size_t stride);

gsl_vector_complex 
*gsl_vector_complex_alloc_from_vector (gsl_vector_complex * v, 
                                             size_t offset, 
                                             size_t n, 
                                             size_t stride);

void gsl_vector_complex_free (gsl_vector_complex * v);

gsl_complex 
*gsl_vector_complex_ptr (const gsl_vector_complex * v, size_t i);

gsl_complex 
gsl_vector_complex_get (const gsl_vector_complex * v, size_t i);

void gsl_vector_complex_set (gsl_vector_complex * v, size_t i,
                                   gsl_complex z);

int gsl_vector_complex_fread (FILE * stream,
				    gsl_vector_complex * v);
int gsl_vector_complex_fwrite (FILE * stream,
				     const gsl_vector_complex * v);
int gsl_vector_complex_fscanf (FILE * stream,
				     gsl_vector_complex * v);
int gsl_vector_complex_fprintf (FILE * stream,
				      const gsl_vector_complex * v,
				      const char *format);

int gsl_vector_complex_copy (gsl_vector_complex * dest, const gsl_vector_complex * src);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
gsl_complex
gsl_vector_complex_get (const gsl_vector_complex * v,
			      const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)
    {
      const gsl_complex zero = {{0, 0}};
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, zero);
    }
#endif
  return *GSL_COMPLEX_AT (v, i);
}

extern inline
void
gsl_vector_complex_set (gsl_vector_complex * v,
			      const size_t i, gsl_complex z)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)
    {
      GSL_ERROR_RETURN_NOTHING ("index out of range", GSL_EINVAL);
    }
#endif
  *GSL_COMPLEX_AT (v, i) = z;
}

#endif /* HAVE_INLINE */

#endif /* GSL_VECTOR_COMPLEX_DOUBLE_H */
