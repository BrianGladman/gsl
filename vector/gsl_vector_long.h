#ifndef __GSL_VECTOR_LONG_H__
#define __GSL_VECTOR_LONG_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block_long.h>

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

struct gsl_vector_long_struct
{
  size_t size;
  size_t stride;
  long *data;
  gsl_block_long *block;
};

typedef struct gsl_vector_long_struct gsl_vector_long;

gsl_vector_long *gsl_vector_long_alloc (const size_t n);
gsl_vector_long *gsl_vector_long_calloc (const size_t n);

gsl_vector_long *gsl_vector_long_alloc_from_block (gsl_block_long * b,
                                                     const size_t offset, 
                                                     const size_t n, 
                                                     const size_t stride);

gsl_vector_long *gsl_vector_long_alloc_from_vector (gsl_vector_long * v,
                                                      const size_t offset, 
                                                      const size_t n, 
                                                      const size_t stride);

void gsl_vector_long_free (gsl_vector_long * v);

int gsl_vector_long_view_from_vector (gsl_vector_long *v, 
                                       gsl_vector_long *base,
                                       size_t offset, size_t n, size_t stride);

gsl_vector_long gsl_vector_long_subvector (gsl_vector_long *v, size_t i, size_t n);
gsl_vector_long gsl_vector_long_subvector_with_stride (gsl_vector_long *v, size_t i, size_t n, size_t stride);

long *gsl_vector_long_ptr (const gsl_vector_long * v, const size_t i);
long gsl_vector_long_get (const gsl_vector_long * v, const size_t i);
void gsl_vector_long_set (gsl_vector_long * v, const size_t i, long x);

void gsl_vector_long_set_zero (gsl_vector_long * v);
void gsl_vector_long_set_all (gsl_vector_long * v, long x);

int gsl_vector_long_fread (FILE * stream, gsl_vector_long * v);
int gsl_vector_long_fwrite (FILE * stream, const gsl_vector_long * v);
int gsl_vector_long_fscanf (FILE * stream, gsl_vector_long * v);
int gsl_vector_long_fprintf (FILE * stream, const gsl_vector_long * v,
			      const char *format);

int gsl_vector_long_memcpy (gsl_vector_long * dest, const gsl_vector_long * src);

int gsl_vector_long_reverse (gsl_vector_long * v);

int gsl_vector_long_swap (gsl_vector_long * v, gsl_vector_long * w);
int gsl_vector_long_swap_elements (gsl_vector_long * v, const size_t i, const size_t j);

long gsl_vector_long_max (const gsl_vector_long * v);
long gsl_vector_long_min (const gsl_vector_long * v);
void gsl_vector_long_minmax (const gsl_vector_long * v, long * min_out, long * max_out);

size_t gsl_vector_long_max_index (const gsl_vector_long * v);
size_t gsl_vector_long_min_index (const gsl_vector_long * v);
void gsl_vector_long_minmax_index (const gsl_vector_long * v, size_t * imin, size_t * imax);


int gsl_vector_long_isnull (gsl_vector_long * v);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
long
gsl_vector_long_get (const gsl_vector_long * v, const size_t i)
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
gsl_vector_long_set (gsl_vector_long * v, const size_t i, long x)
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

#endif /* __GSL_VECTOR_LONG_H__ */


