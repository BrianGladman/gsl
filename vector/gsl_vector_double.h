#ifndef __GSL_VECTOR_DOUBLE_H__
#define __GSL_VECTOR_DOUBLE_H__

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block_double.h>

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

struct gsl_vector_struct
{
  size_t size;
  size_t stride;
  double *data;
  gsl_block *block;
};

typedef struct gsl_vector_struct gsl_vector;

gsl_vector *gsl_vector_alloc (const size_t n);
gsl_vector *gsl_vector_calloc (const size_t n);

gsl_vector *gsl_vector_alloc_from_block (gsl_block * b,
                                                     const size_t offset, 
                                                     const size_t n, 
                                                     const size_t stride);

gsl_vector *gsl_vector_alloc_from_vector (gsl_vector * v,
                                                      const size_t offset, 
                                                      const size_t n, 
                                                      const size_t stride);

void gsl_vector_free (gsl_vector * v);

int gsl_vector_view_from_vector (gsl_vector *v, 
                                       gsl_vector *base,
                                       size_t offset, size_t n, size_t stride);

double *gsl_vector_ptr (const gsl_vector * v, const size_t i);
double gsl_vector_get (const gsl_vector * v, const size_t i);
void gsl_vector_set (gsl_vector * v, const size_t i, double x);

void gsl_vector_set_zero (gsl_vector * v);
void gsl_vector_set_all (gsl_vector * v, double x);

int gsl_vector_fread (FILE * stream, gsl_vector * v);
int gsl_vector_fwrite (FILE * stream, const gsl_vector * v);
int gsl_vector_fscanf (FILE * stream, gsl_vector * v);
int gsl_vector_fprintf (FILE * stream, const gsl_vector * v,
			      const char *format);

int gsl_vector_memcpy (gsl_vector * dest, const gsl_vector * src);

int gsl_vector_reverse (gsl_vector * v);

int gsl_vector_swap (gsl_vector * v, const size_t i, const size_t j);

int gsl_vector_isnull (gsl_vector * v);

extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline
double
gsl_vector_get (const gsl_vector * v, const size_t i)
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
gsl_vector_set (gsl_vector * v, const size_t i, double x)
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

#endif /* __GSL_VECTOR_DOUBLE_H__ */


