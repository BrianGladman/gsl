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
};

typedef struct gsl_vector_struct gsl_vector;

gsl_vector *gsl_vector_alloc (gsl_block * b,
                                                     size_t offset, 
                                                     size_t n, 
                                                     size_t stride);

gsl_vector *gsl_vector_alloc_from_vector (gsl_vector * v,
                                                      size_t offset, 
                                                      size_t n, 
                                                      size_t stride);

void gsl_vector_free (gsl_vector * v);

double *gsl_vector_ptr (const gsl_vector * v, const size_t i);
double gsl_vector_get (const gsl_vector * v, const size_t i);
void gsl_vector_set (gsl_vector * v, const size_t i, double x);

int gsl_vector_fread (FILE * stream, gsl_vector * v);
int gsl_vector_fwrite (FILE * stream, const gsl_vector * v);
int gsl_vector_fscanf (FILE * stream, gsl_vector * v);
int gsl_vector_fprintf (FILE * stream, const gsl_vector * v,
			      const char *format);

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








