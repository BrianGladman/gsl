#ifndef __GSL_VECTOR_COMPLEX_DOUBLE_H__
#define __GSL_VECTOR_COMPLEX_DOUBLE_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_block_complex_double.h>

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

struct gsl_vector_complex_struct
{
  size_t size;
  size_t stride;
  double *data;
  gsl_block_complex *block;
};

typedef struct gsl_vector_complex_struct gsl_vector_complex;

gsl_vector_complex *gsl_vector_complex_alloc (const size_t n);
gsl_vector_complex *gsl_vector_complex_calloc (const size_t n);

gsl_vector_complex 
*gsl_vector_complex_alloc_from_block (gsl_block_complex * b, 
                                            const size_t offset, 
                                            const size_t n, 
                                            const size_t stride);

gsl_vector_complex 
*gsl_vector_complex_alloc_from_vector (gsl_vector_complex * v, 
                                             const size_t offset, 
                                             const size_t n, 
                                             const size_t stride);

void gsl_vector_complex_free (gsl_vector_complex * v);

int gsl_vector_complex_view_from_vector (gsl_vector_complex *v, 
                                               gsl_vector_complex *base,
                                               size_t offset, size_t n, size_t stride);

gsl_vector_complex gsl_vector_complex_subvector (gsl_vector_complex *v, size_t i, size_t n);
gsl_vector_complex gsl_vector_complex_subvector_with_stride (gsl_vector_complex *v, size_t i, size_t n, size_t stride);


gsl_complex 
*gsl_vector_complex_ptr (const gsl_vector_complex * v, size_t i);

gsl_complex 
gsl_vector_complex_get (const gsl_vector_complex * v, const size_t i);

void gsl_vector_complex_set (gsl_vector_complex * v, const size_t i,
                                   gsl_complex z);

void gsl_vector_complex_set_zero (gsl_vector_complex * v);
void gsl_vector_complex_set_all (gsl_vector_complex * v,
                                       gsl_complex z);
int gsl_vector_complex_set_basis (gsl_vector_complex * v, size_t i);

int gsl_vector_complex_fread (FILE * stream,
				    gsl_vector_complex * v);
int gsl_vector_complex_fwrite (FILE * stream,
				     const gsl_vector_complex * v);
int gsl_vector_complex_fscanf (FILE * stream,
				     gsl_vector_complex * v);
int gsl_vector_complex_fprintf (FILE * stream,
				      const gsl_vector_complex * v,
				      const char *format);

int gsl_vector_complex_memcpy (gsl_vector_complex * dest, const gsl_vector_complex * src);

int gsl_vector_complex_reverse (gsl_vector_complex * v);

int gsl_vector_complex_swap (gsl_vector_complex * v, gsl_vector_complex * w);
int gsl_vector_complex_swap_elements (gsl_vector_complex * v, const size_t i, const size_t j);

int gsl_vector_complex_isnull (const gsl_vector_complex * v);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
gsl_complex
gsl_vector_complex_get (const gsl_vector_complex * v,
			      const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)
    {
      const gsl_complex zero = {{0, 0}};
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, zero);
    }
#endif
  return *GSL_COMPLEX_AT (v, i);
}

extern inline
void
gsl_vector_complex_set (gsl_vector_complex * v,
			      const size_t i, gsl_complex z)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)
    {
      GSL_ERROR_RETURN_NOTHING ("index out of range", GSL_EINVAL);
    }
#endif
  *GSL_COMPLEX_AT (v, i) = z;
}

#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_VECTOR_COMPLEX_DOUBLE_H__ */
