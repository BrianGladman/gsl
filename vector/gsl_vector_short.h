#ifndef GSL_VECTOR_SHORT_H
#define GSL_VECTOR_SHORT_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>

struct gsl_matrix_short_struct;
typedef struct gsl_matrix_short_struct gsl_matrix_short;

struct gsl_vector_short_struct
{
  size_t size;
  size_t stride;
  gsl_matrix_short * parent;
  short *data;
};

typedef struct gsl_vector_short_struct gsl_vector_short;

gsl_vector_short *gsl_vector_short_alloc (size_t n);
gsl_vector_short *gsl_vector_short_calloc (size_t n);
void gsl_vector_short_free (gsl_vector_short * v);

short *gsl_vector_short_ptr (const gsl_vector_short * v, const size_t i);
short gsl_vector_short_get (const gsl_vector_short * v, const size_t i);
void gsl_vector_short_set (gsl_vector_short * v, const size_t i, short x);

int gsl_vector_short_fread (FILE * stream, gsl_vector_short * v);
int gsl_vector_short_fwrite (FILE * stream, const gsl_vector_short * v);
int gsl_vector_short_fscanf (FILE * stream, gsl_vector_short * v);
int gsl_vector_short_fprintf (FILE * stream, const gsl_vector_short * v,
			      const char *format);

int gsl_block_short_fread (FILE * stream, short *data, size_t n, size_t stride);
int gsl_block_short_fwrite (FILE * stream, const short *data, size_t n, size_t stride);
int gsl_block_short_fscanf (FILE * stream, short *data, size_t n, size_t stride);
int gsl_block_short_fprintf (FILE * stream, const short *data, size_t n, size_t stride,
			     const char *format);

extern int gsl_check_range;



#ifdef HAVE_INLINE

extern inline
short *
gsl_vector_short_ptr (const gsl_vector_short * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return (short *) (v->data + i * v->stride);
}

extern inline
short
gsl_vector_short_get (const gsl_vector_short * v, const size_t i)
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
gsl_vector_short_set (gsl_vector_short * v, const size_t i, short x)
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

#endif /* GSL_VECTOR_SHORT_H */
