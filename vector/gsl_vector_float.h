#ifndef __GSL_VECTOR_FLOAT_H__
#define __GSL_VECTOR_FLOAT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block_float.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

struct gsl_vector_float_struct
{
  size_t size;
  size_t stride;
  float *data;
  gsl_block_float *block;
};

typedef struct gsl_vector_float_struct gsl_vector_float;

gsl_vector_float *gsl_vector_float_alloc (const size_t n);
gsl_vector_float *gsl_vector_float_calloc (const size_t n);

gsl_vector_float *gsl_vector_float_alloc_from_block (gsl_block_float * b,
                                                     const size_t offset, 
                                                     const size_t n, 
                                                     const size_t stride);

gsl_vector_float *gsl_vector_float_alloc_from_vector (gsl_vector_float * v,
                                                      const size_t offset, 
                                                      const size_t n, 
                                                      const size_t stride);

void gsl_vector_float_free (gsl_vector_float * v);

int gsl_vector_float_view_from_vector (gsl_vector_float *v, 
                                       gsl_vector_float *base,
                                       size_t offset, size_t n, size_t stride);

gsl_vector_float gsl_vector_float_subvector (gsl_vector_float *v, size_t i, size_t n, size_t stride);

float *gsl_vector_float_ptr (const gsl_vector_float * v, const size_t i);
float gsl_vector_float_get (const gsl_vector_float * v, const size_t i);
void gsl_vector_float_set (gsl_vector_float * v, const size_t i, float x);

void gsl_vector_float_set_zero (gsl_vector_float * v);
void gsl_vector_float_set_all (gsl_vector_float * v, float x);

int gsl_vector_float_fread (FILE * stream, gsl_vector_float * v);
int gsl_vector_float_fwrite (FILE * stream, const gsl_vector_float * v);
int gsl_vector_float_fscanf (FILE * stream, gsl_vector_float * v);
int gsl_vector_float_fprintf (FILE * stream, const gsl_vector_float * v,
			      const char *format);

int gsl_vector_float_memcpy (gsl_vector_float * dest, const gsl_vector_float * src);

int gsl_vector_float_reverse (gsl_vector_float * v);

int gsl_vector_float_swap (gsl_vector_float * v, gsl_vector_float * w);
int gsl_vector_float_swap_elements (gsl_vector_float * v, const size_t i, const size_t j);

int gsl_vector_float_isnull (gsl_vector_float * v);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
float
gsl_vector_float_get (const gsl_vector_float * v, const size_t i)
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
gsl_vector_float_set (gsl_vector_float * v, const size_t i, float x)
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

__END_DECLS

#endif /* __GSL_VECTOR_FLOAT_H__ */


