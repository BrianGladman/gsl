#ifndef GSL_VECTOR_UCHAR_H
#define GSL_VECTOR_UCHAR_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_block_uchar.h>

struct gsl_vector_uchar_struct
{
  size_t size;
  size_t stride;
  unsigned char *data;
  gsl_block_uchar *block;
};

typedef struct gsl_vector_uchar_struct gsl_vector_uchar;

gsl_vector_uchar *gsl_vector_uchar_alloc (size_t n);
gsl_vector_uchar *gsl_vector_uchar_calloc (size_t n);

gsl_vector_uchar *gsl_vector_uchar_alloc_from_block (gsl_block_uchar * b,
                                                     size_t offset, 
                                                     size_t n, 
                                                     size_t stride);

gsl_vector_uchar *gsl_vector_uchar_alloc_from_vector (gsl_vector_uchar * v,
                                                      size_t offset, 
                                                      size_t n, 
                                                      size_t stride);

void gsl_vector_uchar_free (gsl_vector_uchar * v);

unsigned char *gsl_vector_uchar_ptr (const gsl_vector_uchar * v, const size_t i);
unsigned char gsl_vector_uchar_get (const gsl_vector_uchar * v, const size_t i);
void gsl_vector_uchar_set (gsl_vector_uchar * v, const size_t i, unsigned char x);

int gsl_vector_uchar_fread (FILE * stream, gsl_vector_uchar * v);
int gsl_vector_uchar_fwrite (FILE * stream, const gsl_vector_uchar * v);
int gsl_vector_uchar_fscanf (FILE * stream, gsl_vector_uchar * v);
int gsl_vector_uchar_fprintf (FILE * stream, const gsl_vector_uchar * v,
			      const char *format);

int gsl_vector_uchar_copy (gsl_vector_uchar * dest, const gsl_vector_uchar * src);

int gsl_vector_uchar_swap (gsl_vector_uchar * v, size_t i, size_t j);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
unsigned char
gsl_vector_uchar_get (const gsl_vector_uchar * v, const size_t i)
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
gsl_vector_uchar_set (gsl_vector_uchar * v, const size_t i, unsigned char x)
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

#endif /* GSL_VECTOR_UCHAR_H */








