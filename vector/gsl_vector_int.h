#ifndef GSL_VECTOR_INT_H
#define GSL_VECTOR_INT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_block_int.h>

struct gsl_vector_int_struct
{
  size_t size;
  size_t stride;
  int *data;
  gsl_block_int *block;
};

typedef struct gsl_vector_int_struct gsl_vector_int;

gsl_vector_int *gsl_vector_int_alloc (size_t n);
gsl_vector_int *gsl_vector_int_calloc (size_t n);

gsl_vector_int *gsl_vector_int_alloc_from_block (gsl_block_int * b,
                                                     size_t offset, 
                                                     size_t n, 
                                                     size_t stride);

gsl_vector_int *gsl_vector_int_alloc_from_vector (gsl_vector_int * v,
                                                      size_t offset, 
                                                      size_t n, 
                                                      size_t stride);

void gsl_vector_int_free (gsl_vector_int * v);

int *gsl_vector_int_ptr (const gsl_vector_int * v, const size_t i);
int gsl_vector_int_get (const gsl_vector_int * v, const size_t i);
void gsl_vector_int_set (gsl_vector_int * v, const size_t i, int x);

int gsl_vector_int_fread (FILE * stream, gsl_vector_int * v);
int gsl_vector_int_fwrite (FILE * stream, const gsl_vector_int * v);
int gsl_vector_int_fscanf (FILE * stream, gsl_vector_int * v);
int gsl_vector_int_fprintf (FILE * stream, const gsl_vector_int * v,
			      const char *format);

int gsl_vector_int_copy (gsl_vector_int * dest, const gsl_vector_int * src);

int gsl_vector_int_swap (gsl_vector_int * v, size_t i, size_t j);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
int
gsl_vector_int_get (const gsl_vector_int * v, const size_t i)
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
gsl_vector_int_set (gsl_vector_int * v, const size_t i, int x)
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

#endif /* GSL_VECTOR_INT_H */








