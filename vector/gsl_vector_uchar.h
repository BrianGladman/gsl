#ifndef GSL_VECTOR_UCHAR_H
#define GSL_VECTOR_UCHAR_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>

struct gsl_matrix_uchar_struct;

struct gsl_vector_uchar_struct
{
  size_t size;
  size_t stride;
  struct gsl_matrix_uchar_struct * parent;
  unsigned char *data;
};

typedef struct gsl_vector_uchar_struct gsl_vector_uchar;

gsl_vector_uchar *gsl_vector_uchar_alloc (size_t n);
gsl_vector_uchar *gsl_vector_uchar_calloc (size_t n);
void gsl_vector_uchar_free (gsl_vector_uchar * v);

unsigned char *gsl_vector_uchar_ptr (const gsl_vector_uchar * v, const size_t i);
unsigned char gsl_vector_uchar_get (const gsl_vector_uchar * v, const size_t i);
void gsl_vector_uchar_set (gsl_vector_uchar * v, const size_t i, unsigned char x);


int gsl_vector_uchar_fread (FILE * stream, gsl_vector_uchar * v);
int gsl_vector_uchar_fwrite (FILE * stream, const gsl_vector_uchar * v);
int gsl_vector_uchar_fscanf (FILE * stream, gsl_vector_uchar * v);
int gsl_vector_uchar_fprintf (FILE * stream, const gsl_vector_uchar * v,
			      const char *format);

int gsl_block_uchar_fread (FILE * stream, unsigned char *data, size_t n, size_t stride);
int gsl_block_uchar_fwrite (FILE * stream, const unsigned char *data, size_t n, size_t stride);
int gsl_block_uchar_fscanf (FILE * stream, unsigned char *data, size_t n, size_t stride);
int gsl_block_uchar_fprintf (FILE * stream, const unsigned char *data, size_t n, size_t stride,
			     const char *format);

extern int gsl_check_range;



#ifdef HAVE_INLINE

extern inline
unsigned char *
gsl_vector_uchar_ptr (const gsl_vector_uchar * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return (unsigned char *) (v->data + i * v->stride);
}

extern inline
unsigned char
gsl_vector_uchar_get (const gsl_vector_uchar * v, const size_t i)
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
gsl_vector_uchar_set (gsl_vector_uchar * v, const size_t i, unsigned char c)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING ("index out of range", GSL_EINVAL);
    }
#endif
  v->data[i * v->stride] = c;
}

#endif /* HAVE_INLINE */

#endif /* GSL_VECTOR_UCHAR_H */
