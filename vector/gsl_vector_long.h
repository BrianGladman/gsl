#ifndef GSL_VECTOR_long_H
#define GSL_VECTOR_long_H

#include <stdlib.h>
#include <gsl_errno.h>

struct gsl_matrix_long_struct;

struct gsl_vector_long_struct
{
  size_t size;
  size_t stride;
  struct gsl_matrix_long_struct * parent;
  long *data;
};

typedef struct gsl_vector_long_struct gsl_vector_long;

gsl_vector_long *gsl_vector_long_alloc (size_t n);
gsl_vector_long *gsl_vector_long_calloc (size_t n);
void gsl_vector_long_free (gsl_vector_long * v);

long *gsl_vector_long_ptr (const gsl_vector_long * v, const size_t i);
long gsl_vector_long_get (const gsl_vector_long * v, const size_t i);
void gsl_vector_long_set (gsl_vector_long * v, const size_t i, long x);


int gsl_vector_long_fread (FILE * stream, gsl_vector_long * v);
int gsl_vector_long_fwrite (FILE * stream, const gsl_vector_long * v);
int gsl_vector_long_fscanf (FILE * stream, gsl_vector_long * v);
int gsl_vector_long_fprintf (FILE * stream, const gsl_vector_long * v,
			     const char *format);

int gsl_block_long_fread (FILE * stream, long *data, size_t n, size_t stride);
int gsl_block_long_fwrite (FILE * stream, const long *data, size_t n, size_t stride);
int gsl_block_long_fscanf (FILE * stream, long *data, size_t n, size_t stride);
int gsl_block_long_fprintf (FILE * stream, const long *data, size_t n, size_t stride,
			    const char *format);

extern int gsl_check_range;



#ifdef HAVE_INLINE

extern inline
long *
gsl_vector_long_ptr (const gsl_vector_long * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return (long *) (v->data + i * v->stride);
}

extern inline
long
gsl_vector_long_get (const gsl_vector_long * v, const size_t i)
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
gsl_vector_long_set (gsl_vector_long * v, const size_t i, long x)
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

#endif /* GSL_VECTOR_LONG_H */
