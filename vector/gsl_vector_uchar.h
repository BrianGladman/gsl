#ifndef __GSL_VECTOR_UCHAR_H__
#define __GSL_VECTOR_UCHAR_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block_uchar.h>

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

struct gsl_vector_uchar_struct
{
  size_t size;
  size_t stride;
  unsigned char *data;
  gsl_block_uchar *block;
};

typedef struct gsl_vector_uchar_struct gsl_vector_uchar;

gsl_vector_uchar *gsl_vector_uchar_alloc (const size_t n);
gsl_vector_uchar *gsl_vector_uchar_calloc (const size_t n);

gsl_vector_uchar *gsl_vector_uchar_alloc_from_block (gsl_block_uchar * b,
                                                     const size_t offset, 
                                                     const size_t n, 
                                                     const size_t stride);

gsl_vector_uchar *gsl_vector_uchar_alloc_from_vector (gsl_vector_uchar * v,
                                                      const size_t offset, 
                                                      const size_t n, 
                                                      const size_t stride);

void gsl_vector_uchar_free (gsl_vector_uchar * v);

int gsl_vector_uchar_view_from_vector (gsl_vector_uchar *v, 
                                       gsl_vector_uchar *base,
                                       size_t offset, size_t n, size_t stride);

gsl_vector_uchar gsl_vector_uchar_subvector (gsl_vector_uchar *v, size_t i, size_t n);
gsl_vector_uchar gsl_vector_uchar_subvector_with_stride (gsl_vector_uchar *v, size_t i, size_t n, size_t stride);

unsigned char *gsl_vector_uchar_ptr (const gsl_vector_uchar * v, const size_t i);
unsigned char gsl_vector_uchar_get (const gsl_vector_uchar * v, const size_t i);
void gsl_vector_uchar_set (gsl_vector_uchar * v, const size_t i, unsigned char x);

void gsl_vector_uchar_set_zero (gsl_vector_uchar * v);
void gsl_vector_uchar_set_all (gsl_vector_uchar * v, unsigned char x);
int gsl_vector_uchar_set_basis (gsl_vector_uchar * v, size_t i);

int gsl_vector_uchar_fread (FILE * stream, gsl_vector_uchar * v);
int gsl_vector_uchar_fwrite (FILE * stream, const gsl_vector_uchar * v);
int gsl_vector_uchar_fscanf (FILE * stream, gsl_vector_uchar * v);
int gsl_vector_uchar_fprintf (FILE * stream, const gsl_vector_uchar * v,
			      const char *format);

int gsl_vector_uchar_memcpy (gsl_vector_uchar * dest, const gsl_vector_uchar * src);

int gsl_vector_uchar_reverse (gsl_vector_uchar * v);

int gsl_vector_uchar_swap (gsl_vector_uchar * v, gsl_vector_uchar * w);
int gsl_vector_uchar_swap_elements (gsl_vector_uchar * v, const size_t i, const size_t j);

unsigned char gsl_vector_uchar_max (const gsl_vector_uchar * v);
unsigned char gsl_vector_uchar_min (const gsl_vector_uchar * v);
void gsl_vector_uchar_minmax (const gsl_vector_uchar * v, unsigned char * min_out, unsigned char * max_out);

size_t gsl_vector_uchar_max_index (const gsl_vector_uchar * v);
size_t gsl_vector_uchar_min_index (const gsl_vector_uchar * v);
void gsl_vector_uchar_minmax_index (const gsl_vector_uchar * v, size_t * imin, size_t * imax);


int gsl_vector_uchar_isnull (gsl_vector_uchar * v);

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

__END_DECLS

#endif /* __GSL_VECTOR_UCHAR_H__ */


