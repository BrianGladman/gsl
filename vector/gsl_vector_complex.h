#ifndef GSL_VECTOR_COMPLEX_H 
#define GSL_VECTOR_COMPLEX_H 

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>
#include <gsl_complex.h>

typedef struct
{
  size_t size;
  size_t stride;
  double * data;
} 
gsl_vector_complex ;


gsl_vector_complex * gsl_vector_complex_alloc (size_t n);
gsl_vector_complex * gsl_vector_complex_calloc (size_t n);
void gsl_vector_complex_free (gsl_vector_complex * v);

gsl_complex * gsl_vector_complex_ptr(const gsl_vector_complex * v, size_t i);

int gsl_vector_complex_fread (FILE * stream, gsl_vector_complex * v) ;
int gsl_vector_complex_fwrite (FILE * stream, const gsl_vector_complex * v) ;
int gsl_vector_complex_fscanf (FILE * stream, gsl_vector_complex * v);
int gsl_vector_complex_fprintf (FILE * stream, const gsl_vector_complex * v,
			        const char * format);

int gsl_block_complex_fread (FILE * stream, double * data, size_t n) ;
int gsl_block_complex_fwrite (FILE * stream, const double * data, size_t n) ;
int gsl_block_complex_fscanf (FILE * stream, double * data, size_t n);
int gsl_block_complex_fprintf (FILE * stream, const double * data, size_t n,
			       const char * format);

extern int gsl_check_range ;

#ifndef  GSL_VECTOR_COMPLEX_REAL
#define  GSL_VECTOR_COMPLEX_REAL(z, i)  (z->data[2*i])
#define  GSL_VECTOR_COMPLEX_IMAG(z, i)  (z->data[2*i + 1])
#endif

#define GSL_COMPLEX_AT(zv, i)  ((gsl_complex *)  &(zv->data[2*i]))


/* inline functions if you are using GCC or otherwise enlightened cc */

#ifdef HAVE_INLINE
extern inline 
gsl_complex *
gsl_vector_complex_ptr(const gsl_vector_complex * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0);
    }
#endif
  return GSL_COMPLEX_AT(v, i);
} 


#endif /* HAVE_INLINE */


#endif /* !GSL_VECTOR_COMPLEX_H */
