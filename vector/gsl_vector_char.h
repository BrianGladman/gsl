#ifndef __GSL_VECTOR_CHAR_H__
#define __GSL_VECTOR_CHAR_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block_char.h>

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

struct gsl_vector_char_struct
{
  size_t size;
  size_t stride;
  char *data;
  gsl_block_char *block;
};

typedef struct gsl_vector_char_struct gsl_vector_char;

gsl_vector_char *gsl_vector_char_alloc (const size_t n);
gsl_vector_char *gsl_vector_char_calloc (const size_t n);

gsl_vector_char *gsl_vector_char_alloc_from_block (gsl_block_char * b,
                                                     const size_t offset, 
                                                     const size_t n, 
                                                     const size_t stride);

gsl_vector_char *gsl_vector_char_alloc_from_vector (gsl_vector_char * v,
                                                      const size_t offset, 
                                                      const size_t n, 
                                                      const size_t stride);

void gsl_vector_char_free (gsl_vector_char * v);

int gsl_vector_char_view_from_vector (gsl_vector_char *v, 
                                       gsl_vector_char *base,
                                       size_t offset, size_t n, size_t stride);

gsl_vector_char gsl_vector_char_subvector (gsl_vector_char *v, size_t i, size_t n);
gsl_vector_char gsl_vector_char_subvector_with_stride (gsl_vector_char *v, size_t i, size_t n, size_t stride);

char *gsl_vector_char_ptr (const gsl_vector_char * v, const size_t i);
char gsl_vector_char_get (const gsl_vector_char * v, const size_t i);
void gsl_vector_char_set (gsl_vector_char * v, const size_t i, char x);

void gsl_vector_char_set_zero (gsl_vector_char * v);
void gsl_vector_char_set_all (gsl_vector_char * v, char x);

int gsl_vector_char_fread (FILE * stream, gsl_vector_char * v);
int gsl_vector_char_fwrite (FILE * stream, const gsl_vector_char * v);
int gsl_vector_char_fscanf (FILE * stream, gsl_vector_char * v);
int gsl_vector_char_fprintf (FILE * stream, const gsl_vector_char * v,
			      const char *format);

int gsl_vector_char_memcpy (gsl_vector_char * dest, const gsl_vector_char * src);

int gsl_vector_char_reverse (gsl_vector_char * v);

int gsl_vector_char_swap (gsl_vector_char * v, gsl_vector_char * w);
int gsl_vector_char_swap_elements (gsl_vector_char * v, const size_t i, const size_t j);

char gsl_vector_char_max (const gsl_vector_char * v);
char gsl_vector_char_min (const gsl_vector_char * v);
void gsl_vector_char_minmax (const gsl_vector_char * v, char * min_out, char * max_out);

size_t gsl_vector_char_max_index (const gsl_vector_char * v);
size_t gsl_vector_char_min_index (const gsl_vector_char * v);
void gsl_vector_char_minmax_index (const gsl_vector_char * v, size_t * imin, size_t * imax);


int gsl_vector_char_isnull (gsl_vector_char * v);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
char
gsl_vector_char_get (const gsl_vector_char * v, const size_t i)
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
gsl_vector_char_set (gsl_vector_char * v, const size_t i, char x)
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

#endif /* __GSL_VECTOR_CHAR_H__ */


