#ifndef __GSL_VECTOR_INT_H__
#define __GSL_VECTOR_INT_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block_int.h>

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

struct gsl_vector_int_struct
{
  size_t size;
  size_t stride;
  int *data;
  gsl_block_int *block;
};

typedef struct gsl_vector_int_struct gsl_vector_int;

gsl_vector_int *gsl_vector_int_alloc (const size_t n);
gsl_vector_int *gsl_vector_int_calloc (const size_t n);

gsl_vector_int *gsl_vector_int_alloc_from_block (gsl_block_int * b,
                                                     const size_t offset, 
                                                     const size_t n, 
                                                     const size_t stride);

gsl_vector_int *gsl_vector_int_alloc_from_vector (gsl_vector_int * v,
                                                      const size_t offset, 
                                                      const size_t n, 
                                                      const size_t stride);

void gsl_vector_int_free (gsl_vector_int * v);

int gsl_vector_int_view_from_vector (gsl_vector_int *v, 
                                       gsl_vector_int *base,
                                       size_t offset, size_t n, size_t stride);

gsl_vector_int gsl_vector_int_subvector (gsl_vector_int *v, size_t i, size_t n, size_t stride);

int *gsl_vector_int_ptr (const gsl_vector_int * v, const size_t i);
int gsl_vector_int_get (const gsl_vector_int * v, const size_t i);
void gsl_vector_int_set (gsl_vector_int * v, const size_t i, int x);

void gsl_vector_int_set_zero (gsl_vector_int * v);
void gsl_vector_int_set_all (gsl_vector_int * v, int x);

int gsl_vector_int_fread (FILE * stream, gsl_vector_int * v);
int gsl_vector_int_fwrite (FILE * stream, const gsl_vector_int * v);
int gsl_vector_int_fscanf (FILE * stream, gsl_vector_int * v);
int gsl_vector_int_fprintf (FILE * stream, const gsl_vector_int * v,
			      const char *format);

int gsl_vector_int_memcpy (gsl_vector_int * dest, const gsl_vector_int * src);

int gsl_vector_int_reverse (gsl_vector_int * v);

int gsl_vector_int_swap (gsl_vector_int * v, gsl_vector_int * w);
int gsl_vector_int_swap_elements (gsl_vector_int * v, const size_t i, const size_t j);

int gsl_vector_int_max (const gsl_vector_int * v);
int gsl_vector_int_min (const gsl_vector_int * v);
void gsl_vector_int_minmax (const gsl_vector_int * v, int * min_out, int * max_out);

size_t gsl_vector_int_max_index (const gsl_vector_int * v);
size_t gsl_vector_int_min_index (const gsl_vector_int * v);
void gsl_vector_int_minmax_index (const gsl_vector_int * v, size_t * imin, size_t * imax);


int gsl_vector_int_isnull (gsl_vector_int * v);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
int
gsl_vector_int_get (const gsl_vector_int * v, const size_t i)
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
gsl_vector_int_set (gsl_vector_int * v, const size_t i, int x)
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

#endif /* __GSL_VECTOR_INT_H__ */


