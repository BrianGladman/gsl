#ifndef GSL_VECTOR_USHORT_H
#define GSL_VECTOR_USHORT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_block_ushort.h>

struct gsl_vector_ushort_struct
{
  size_t size;
  size_t stride;
  unsigned short *data;
};

typedef struct gsl_vector_ushort_struct gsl_vector_ushort;

gsl_vector_ushort *gsl_vector_ushort_alloc (gsl_block_ushort * b,
                                                     size_t offset, 
                                                     size_t n, 
                                                     size_t stride);

gsl_vector_ushort *gsl_vector_ushort_alloc_from_vector (gsl_vector_ushort * v,
                                                      size_t offset, 
                                                      size_t n, 
                                                      size_t stride);

void gsl_vector_ushort_free (gsl_vector_ushort * v);

unsigned short *gsl_vector_ushort_ptr (const gsl_vector_ushort * v, const size_t i);
unsigned short gsl_vector_ushort_get (const gsl_vector_ushort * v, const size_t i);
void gsl_vector_ushort_set (gsl_vector_ushort * v, const size_t i, unsigned short x);

int gsl_vector_ushort_fread (FILE * stream, gsl_vector_ushort * v);
int gsl_vector_ushort_fwrite (FILE * stream, const gsl_vector_ushort * v);
int gsl_vector_ushort_fscanf (FILE * stream, gsl_vector_ushort * v);
int gsl_vector_ushort_fprintf (FILE * stream, const gsl_vector_ushort * v,
			      const char *format);

int gsl_vector_ushort_copy (gsl_vector_ushort * dest, const gsl_vector_ushort * src);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
unsigned short
gsl_vector_ushort_get (const gsl_vector_ushort * v, const size_t i)
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
gsl_vector_ushort_set (gsl_vector_ushort * v, const size_t i, unsigned short x)
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

#endif /* GSL_VECTOR_USHORT_H */








