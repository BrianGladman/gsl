#ifndef __GSL_VECTOR_LONG_DOUBLE_H__
#define __GSL_VECTOR_LONG_DOUBLE_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block_long_double.h>

struct gsl_vector_long_double_struct
{
  size_t size;
  size_t stride;
  long double *data;
  gsl_block_long_double *block;
};

typedef struct gsl_vector_long_double_struct gsl_vector_long_double;

gsl_vector_long_double *gsl_vector_long_double_alloc (size_t n);
gsl_vector_long_double *gsl_vector_long_double_calloc (size_t n);

gsl_vector_long_double *gsl_vector_long_double_alloc_from_block (gsl_block_long_double * b,
                                                     size_t offset, 
                                                     size_t n, 
                                                     size_t stride);

gsl_vector_long_double *gsl_vector_long_double_alloc_from_vector (gsl_vector_long_double * v,
                                                      size_t offset, 
                                                      size_t n, 
                                                      size_t stride);

void gsl_vector_long_double_free (gsl_vector_long_double * v);

int gsl_vector_long_double_view_from_vector (gsl_vector_long_double *v, 
                                       gsl_vector_long_double *base,
                                       size_t offset, size_t n, size_t stride);

long double *gsl_vector_long_double_ptr (const gsl_vector_long_double * v, const size_t i);
long double gsl_vector_long_double_get (const gsl_vector_long_double * v, const size_t i);
void gsl_vector_long_double_set (gsl_vector_long_double * v, const size_t i, long double x);

void gsl_vector_long_double_set_all (gsl_vector_long_double * v, long double x);

int gsl_vector_long_double_fread (FILE * stream, gsl_vector_long_double * v);
int gsl_vector_long_double_fwrite (FILE * stream, const gsl_vector_long_double * v);
int gsl_vector_long_double_fscanf (FILE * stream, gsl_vector_long_double * v);
int gsl_vector_long_double_fprintf (FILE * stream, const gsl_vector_long_double * v,
			      const char *format);

int gsl_vector_long_double_memcpy (gsl_vector_long_double * dest, const gsl_vector_long_double * src);

int gsl_vector_long_double_reverse (gsl_vector_long_double * v);

int gsl_vector_long_double_swap (gsl_vector_long_double * v, size_t i, size_t j);

int gsl_vector_long_double_isnull (gsl_vector_long_double * v);

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

#endif /* __GSL_VECTOR_LONG_DOUBLE_H__ */








