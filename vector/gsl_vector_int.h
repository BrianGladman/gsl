#ifndef GSL_VECTOR_INT_H
#define GSL_VECTOR_INT_H

#include <stdlib.h>
#include <gsl_errno.h>

struct gsl_matrix_int_struct;

struct gsl_vector_int_struct
{
  size_t size;
  size_t stride;
  struct gsl_matrix_int_struct * parent;
  int *data;
};

typedef struct gsl_vector_int_struct gsl_vector_int;

gsl_vector_int *gsl_vector_int_alloc (size_t n);
gsl_vector_int *gsl_vector_int_calloc (size_t n);
void gsl_vector_int_free (gsl_vector_int * v);

int *gsl_vector_int_ptr (const gsl_vector_int * v, const size_t i);
int gsl_vector_int_get (const gsl_vector_int * v, const size_t i);
void gsl_vector_int_set (gsl_vector_int * v, const size_t i, int x);

int gsl_vector_int_fread (FILE * stream, gsl_vector_int * v);
int gsl_vector_int_fwrite (FILE * stream, const gsl_vector_int * v);
int gsl_vector_int_fscanf (FILE * stream, gsl_vector_int * v);
int gsl_vector_int_fprintf (FILE * stream, const gsl_vector_int * v,
			    const char *format);

int gsl_block_int_fread (FILE * stream, int *data, size_t n, size_t stride);
int gsl_block_int_fwrite (FILE * stream, const int *data, size_t n,
			  size_t stride);
int gsl_block_int_fscanf (FILE * stream, int *data, size_t n, size_t stride);
int gsl_block_int_fprintf (FILE * stream, const int *data, size_t n,
			   size_t stride, const char *format);

extern int gsl_check_range;

#ifdef HAVE_INLINE
extern inline
int *
gsl_vector_int_ptr (const gsl_vector_int * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return (int *) (v->data + i * v->stride);
}

extern inline
int
gsl_vector_int_get (const gsl_vector_int * v, const size_t i)
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
gsl_vector_int_set (gsl_vector_int * v, const size_t i, int x)
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

#endif /* GSL_VECTOR_INT_H */
