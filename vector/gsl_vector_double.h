#ifndef GSL_VECTOR_DOUBLE_H
#define GSL_VECTOR_DOUBLE_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_block_double.h>

struct gsl_vector_struct
{
  size_t size;
  size_t stride;
  double *data;
  gsl_block *block;
};

typedef struct gsl_vector_struct gsl_vector;

gsl_vector *gsl_vector_alloc (size_t n);
gsl_vector *gsl_vector_calloc (size_t n);

gsl_vector *gsl_vector_alloc_from_block (gsl_block * b,
                                                     size_t offset, 
                                                     size_t n, 
                                                     size_t stride);

gsl_vector *gsl_vector_alloc_from_vector (gsl_vector * v,
                                                      size_t offset, 
                                                      size_t n, 
                                                      size_t stride);

void gsl_vector_free (gsl_vector * v);

int gsl_vector_view_from_vector (gsl_vector *v, 
                                       gsl_vector *base,
                                       size_t offset, size_t n, size_t stride);

double *gsl_vector_ptr (const gsl_vector * v, const size_t i);
double gsl_vector_get (const gsl_vector * v, const size_t i);
void gsl_vector_set (gsl_vector * v, const size_t i, double x);

void gsl_vector_set_all (gsl_vector * v, double x);

int gsl_vector_fread (FILE * stream, gsl_vector * v);
int gsl_vector_fwrite (FILE * stream, const gsl_vector * v);
int gsl_vector_fscanf (FILE * stream, gsl_vector * v);
int gsl_vector_fprintf (FILE * stream, const gsl_vector * v,
			      const char *format);

int gsl_vector_copy (gsl_vector * dest, const gsl_vector * src);

int gsl_vector_swap (gsl_vector * v, size_t i, size_t j);

int gsl_vector_isnull (gsl_vector * v);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
double
gsl_vector_get (const gsl_vector * v, const size_t i)
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
gsl_vector_set (gsl_vector * v, const size_t i, double x)
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

#endif /* GSL_VECTOR_DOUBLE_H */








