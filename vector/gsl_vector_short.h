#ifndef __GSL_VECTOR_SHORT_H__
#define __GSL_VECTOR_SHORT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block_short.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

struct gsl_vector_short_struct
{
  size_t size;
  size_t stride;
  short *data;
  gsl_block_short *block;
};

typedef struct gsl_vector_short_struct gsl_vector_short;

gsl_vector_short *gsl_vector_short_alloc (const size_t n);
gsl_vector_short *gsl_vector_short_calloc (const size_t n);

gsl_vector_short *gsl_vector_short_alloc_from_block (gsl_block_short * b,
                                                     const size_t offset, 
                                                     const size_t n, 
                                                     const size_t stride);

gsl_vector_short *gsl_vector_short_alloc_from_vector (gsl_vector_short * v,
                                                      const size_t offset, 
                                                      const size_t n, 
                                                      const size_t stride);

void gsl_vector_short_free (gsl_vector_short * v);

int gsl_vector_short_view_from_vector (gsl_vector_short *v, 
                                       gsl_vector_short *base,
                                       size_t offset, size_t n, size_t stride);

gsl_vector_short gsl_vector_short_subvector (gsl_vector_short *v, size_t i, size_t n);
gsl_vector_short gsl_vector_short_subvector_with_stride (gsl_vector_short *v, size_t i, size_t n, size_t stride);

short *gsl_vector_short_ptr (const gsl_vector_short * v, const size_t i);
short gsl_vector_short_get (const gsl_vector_short * v, const size_t i);
void gsl_vector_short_set (gsl_vector_short * v, const size_t i, short x);

void gsl_vector_short_set_zero (gsl_vector_short * v);
void gsl_vector_short_set_all (gsl_vector_short * v, short x);
int gsl_vector_short_set_basis (gsl_vector_short * v, size_t i);

int gsl_vector_short_fread (FILE * stream, gsl_vector_short * v);
int gsl_vector_short_fwrite (FILE * stream, const gsl_vector_short * v);
int gsl_vector_short_fscanf (FILE * stream, gsl_vector_short * v);
int gsl_vector_short_fprintf (FILE * stream, const gsl_vector_short * v,
			      const char *format);

int gsl_vector_short_memcpy (gsl_vector_short * dest, const gsl_vector_short * src);

int gsl_vector_short_reverse (gsl_vector_short * v);

int gsl_vector_short_swap (gsl_vector_short * v, gsl_vector_short * w);
int gsl_vector_short_swap_elements (gsl_vector_short * v, const size_t i, const size_t j);

short gsl_vector_short_max (const gsl_vector_short * v);
short gsl_vector_short_min (const gsl_vector_short * v);
void gsl_vector_short_minmax (const gsl_vector_short * v, short * min_out, short * max_out);

size_t gsl_vector_short_max_index (const gsl_vector_short * v);
size_t gsl_vector_short_min_index (const gsl_vector_short * v);
void gsl_vector_short_minmax_index (const gsl_vector_short * v, size_t * imin, size_t * imax);


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

__END_DECLS

#endif /* __GSL_VECTOR_SHORT_H__ */


