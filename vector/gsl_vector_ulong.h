#ifndef GSL_VECTOR_ULONG_H
#define GSL_VECTOR_ULONG_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_block_ulong.h>

struct gsl_vector_ulong_struct
{
  size_t size;
  size_t stride;
  unsigned long *data;
};

typedef struct gsl_vector_ulong_struct gsl_vector_ulong;

gsl_vector_ulong *gsl_vector_ulong_alloc (gsl_block_ulong * b,
                                                     size_t offset, 
                                                     size_t n, 
                                                     size_t stride);

gsl_vector_ulong *gsl_vector_ulong_alloc_from_vector (gsl_vector_ulong * v,
                                                      size_t offset, 
                                                      size_t n, 
                                                      size_t stride);

void gsl_vector_ulong_free (gsl_vector_ulong * v);

unsigned long *gsl_vector_ulong_ptr (const gsl_vector_ulong * v, const size_t i);
unsigned long gsl_vector_ulong_get (const gsl_vector_ulong * v, const size_t i);
void gsl_vector_ulong_set (gsl_vector_ulong * v, const size_t i, unsigned long x);

int gsl_vector_ulong_fread (FILE * stream, gsl_vector_ulong * v);
int gsl_vector_ulong_fwrite (FILE * stream, const gsl_vector_ulong * v);
int gsl_vector_ulong_fscanf (FILE * stream, gsl_vector_ulong * v);
int gsl_vector_ulong_fprintf (FILE * stream, const gsl_vector_ulong * v,
			      const char *format);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
unsigned long
gsl_vector_ulong_get (const gsl_vector_ulong * v, const size_t i)
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
gsl_vector_ulong_set (gsl_vector_ulong * v, const size_t i, unsigned long x)
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

#endif /* GSL_VECTOR_ULONG_H */








