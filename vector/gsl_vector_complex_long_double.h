#ifndef GSL_VECTOR_COMPLEX_LONG_H
#define GSL_VECTOR_COMPLEX_LONG_H

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_vector_complex.h>

struct gsl_matrix_complex_long_double_struct;

typedef struct gsl_vector_complex_long_double_struct gsl_vector_complex_long_double;

struct gsl_vector_complex_long_double_struct
{
  size_t size;
  size_t stride;
  struct gsl_matrix_complex_long_double_struct * parent;
  long double *data;
};

gsl_vector_complex_long_double *gsl_vector_complex_long_double_alloc (size_t n);
gsl_vector_complex_long_double *gsl_vector_complex_long_double_calloc (size_t n);
void gsl_vector_complex_long_double_free (gsl_vector_complex_long_double * v);

gsl_complex_long_double *gsl_vector_complex_long_double_ptr (const gsl_vector_complex_long_double * v, size_t i);
gsl_complex_long_double gsl_vector_complex_long_double_get (const gsl_vector_complex_long_double * v, size_t i);
void gsl_vector_complex_long_double_set (gsl_vector_complex_long_double * v, size_t i, gsl_complex_long_double z);

int
  gsl_vector_complex_long_double_fread (FILE * stream,
					gsl_vector_complex_long_double * v);
int
  gsl_vector_complex_long_double_fwrite (FILE * stream,
				  const gsl_vector_complex_long_double * v);
int
  gsl_vector_complex_long_double_fscanf (FILE * stream,
					 gsl_vector_complex_long_double * v);
int
  gsl_vector_complex_long_double_fprintf (FILE * stream,
				   const gsl_vector_complex_long_double * v,
					  const char *format);

int
  gsl_block_complex_long_double_fread (FILE * stream, long double *data,
				       size_t n, size_t stride);
int
  gsl_block_complex_long_double_fwrite (FILE * stream, const long double *data,
					size_t n, size_t stride);
int
  gsl_block_complex_long_double_fscanf (FILE * stream, long double *data,
					size_t n, size_t stride);
int
  gsl_block_complex_long_double_fprintf (FILE * stream, const long double *data,
					 size_t n, size_t stride,
					 const char *format);
extern int gsl_check_range;

#ifdef HAVE_INLINE

extern inline 
gsl_complex_long_double *
gsl_vector_complex_long_double_ptr (const gsl_vector_complex_long_double * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, 0);
    }
#endif
  return GSL_COMPLEX_LONG_DOUBLE_AT (v, i);
}

extern inline
  gsl_complex_long_double
gsl_vector_complex_long_double_get (const gsl_vector_complex_long_double * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      const gsl_complex_long_double zero =
      {
	{0, 0}};
      GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, zero);
    }
#endif
  return *GSL_COMPLEX_LONG_DOUBLE_AT (v, i);
}

extern inline
void
gsl_vector_complex_long_double_set (gsl_vector_complex_long_double * v, const size_t i, gsl_complex_long_double z)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)	/* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING ("index out of range", GSL_EINVAL);
    }
#endif
  *(GSL_COMPLEX_LONG_DOUBLE_AT (v, i)) = z;
}

#endif /* HAVE_INLINE */

#endif /* GSL_VECTOR_COMPLEX_LONG_H */
