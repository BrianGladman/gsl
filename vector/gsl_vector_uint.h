#ifndef GSL_VECTOR_UINT_H
#define GSL_VECTOR_UINT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>

typedef struct
  {
    size_t size;
    size_t stride;
    unsigned int *data;
  }
gsl_vector_uint;

gsl_vector_uint *gsl_vector_uint_alloc (size_t n);
gsl_vector_uint *gsl_vector_uint_calloc (size_t n);
void gsl_vector_uint_free (gsl_vector_uint * v);

unsigned int *gsl_vector_uint_ptr (const gsl_vector_uint * v, const size_t i);
unsigned int gsl_vector_uint_get (const gsl_vector_uint * v, const size_t i);
void gsl_vector_uint_set (gsl_vector_uint * v, const size_t i, unsigned int x);


int gsl_vector_uint_fread (FILE * stream, gsl_vector_uint * v);
int gsl_vector_uint_fwrite (FILE * stream, const gsl_vector_uint * v);
int gsl_vector_uint_fscanf (FILE * stream, gsl_vector_uint * v);
int gsl_vector_uint_fprintf (FILE * stream, const gsl_vector_uint * v,
			     const char *format);

int gsl_block_uint_fread (FILE * stream, unsigned int *data, size_t n, size_t stride);
int gsl_block_uint_fwrite (FILE * stream, const unsigned int *data, size_t n, size_t stride);
int gsl_block_uint_fscanf (FILE * stream, unsigned int *data, size_t n, size_t stride);
int gsl_block_uint_fprintf (FILE * stream, const unsigned int *data, size_t n, size_t stride,
			    const char *format);

extern int gsl_check_range;



#ifdef HAVE_INLINE

extern inline
unsigned int *
gsl_vector_uint_ptr (const gsl_vector_uint * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return (unsigned int *) (v->data + i * v->stride);
}

extern inline
unsigned int
gsl_vector_uint_get (const gsl_vector_uint * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return v->data[i * v->stride];
}

extern inline
void
gsl_vector_uint_set (gsl_vector_uint * v, const size_t i, unsigned int x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING ("index out of range", GSL_EINVAL);
    }
#endif
  v->data[i * v->stride] = x;
}

#endif /* HAVE_INLINE */

#endif /* GSL_VECTOR_UINT_H */
