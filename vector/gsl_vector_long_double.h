#ifndef GSL_VECTOR_LONG_DOUBLE_H
#define GSL_VECTOR_LONG_DOUBLE_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_block_long_double.h>

struct gsl_vector_long_double_struct
{
  size_t size;
  size_t stride;
  long double *data;
};

typedef struct gsl_vector_long_double_struct gsl_vector_long_double;

gsl_vector_long_double *gsl_vector_long_double_alloc (gsl_block_long_double * b,
                                                     size_t offset, 
                                                     size_t n, 
                                                     size_t stride);

gsl_vector_long_double *gsl_vector_long_double_alloc_from_vector (gsl_vector_long_double * v,
                                                      size_t offset, 
                                                      size_t n, 
                                                      size_t stride);

void gsl_vector_long_double_free (gsl_vector_long_double * v);

long double *gsl_vector_long_double_ptr (const gsl_vector_long_double * v, const size_t i);
long double gsl_vector_long_double_get (const gsl_vector_long_double * v, const size_t i);
void gsl_vector_long_double_set (gsl_vector_long_double * v, const size_t i, long double x);

int gsl_vector_long_double_fread (FILE * stream, gsl_vector_long_double * v);
int gsl_vector_long_double_fwrite (FILE * stream, const gsl_vector_long_double * v);
int gsl_vector_long_double_fscanf (FILE * stream, gsl_vector_long_double * v);
int gsl_vector_long_double_fprintf (FILE * stream, const gsl_vector_long_double * v,
			      const char *format);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
long double
gsl_vector_long_double_get (const gsl_vector_long_double * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return v->data[i * v->stride];
}

extern inline
void
gsl_vector_long_double_set (gsl_vector_long_double * v, const size_t i, long double x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)
    {
      GSL_ERROR_RETURN_NOTHING ("index out of range", GSL_EINVAL);
    }
#endif
  v->data[i * v->stride] = x;
}

#endif /* HAVE_INLINE */

#endif /* GSL_VECTOR_LONG_DOUBLE_H */








