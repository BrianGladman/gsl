#ifndef GSL_VECTOR_SHORT_H
#define GSL_VECTOR_SHORT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_block_short.h>

struct gsl_vector_short_struct
{
  size_t size;
  size_t stride;
  short *data;
  gsl_block_short *block;
};

typedef struct gsl_vector_short_struct gsl_vector_short;

gsl_vector_short *gsl_vector_short_alloc (size_t n);
gsl_vector_short *gsl_vector_short_calloc (size_t n);

gsl_vector_short *gsl_vector_short_alloc_from_block (gsl_block_short * b,
                                                     size_t offset, 
                                                     size_t n, 
                                                     size_t stride);

gsl_vector_short *gsl_vector_short_alloc_from_vector (gsl_vector_short * v,
                                                      size_t offset, 
                                                      size_t n, 
                                                      size_t stride);

void gsl_vector_short_free (gsl_vector_short * v);

int gsl_vector_short_view_from_vector (gsl_vector_short *v, 
                                       gsl_vector_short *base,
                                       size_t offset, size_t n, size_t stride);

short *gsl_vector_short_ptr (const gsl_vector_short * v, const size_t i);
short gsl_vector_short_get (const gsl_vector_short * v, const size_t i);
void gsl_vector_short_set (gsl_vector_short * v, const size_t i, short x);

void gsl_vector_short_set_all (gsl_vector_short * v, short x);

int gsl_vector_short_fread (FILE * stream, gsl_vector_short * v);
int gsl_vector_short_fwrite (FILE * stream, const gsl_vector_short * v);
int gsl_vector_short_fscanf (FILE * stream, gsl_vector_short * v);
int gsl_vector_short_fprintf (FILE * stream, const gsl_vector_short * v,
			      const char *format);

int gsl_vector_short_copy (gsl_vector_short * dest, const gsl_vector_short * src);

int gsl_vector_short_swap (gsl_vector_short * v, size_t i, size_t j);

int gsl_vector_short_isnull (gsl_vector_short * v);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
short
gsl_vector_short_get (const gsl_vector_short * v, const size_t i)
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
gsl_vector_short_set (gsl_vector_short * v, const size_t i, short x)
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

#endif /* GSL_VECTOR_SHORT_H */








