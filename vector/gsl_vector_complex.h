#ifndef GSL_VECTOR_COMPLEX_H 
#define GSL_VECTOR_COMPLEX_H 

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>
#include <gsl_complex.h>

typedef struct
{
  size_t size;
  gsl_complex * data;
} 
gsl_vector_complex ;

gsl_vector_complex * gsl_vector_complex_alloc (size_t n);
gsl_vector_complex * gsl_vector_complex_calloc (size_t n);
void gsl_vector_complex_free (gsl_vector_complex * v);

gsl_complex gsl_vector_complex_get(const gsl_vector_complex * v, size_t i);
void gsl_vector_complex_set(gsl_vector_complex * v, size_t i, gsl_complex x);

int gsl_vector_complex_fread (FILE * stream, gsl_vector_complex * v) ;
int gsl_vector_complex_fwrite (FILE * stream, const gsl_vector_complex * v) ;
int gsl_vector_complex_fscanf (FILE * stream, gsl_vector_complex * v);
int gsl_vector_complex_fprintf (FILE * stream, const gsl_vector_complex * v,
			      const char * format);

int gsl_block_complex_fread (FILE * stream, gsl_complex * data, size_t n) ;
int gsl_block_complex_fwrite (FILE * stream, const gsl_complex * data, size_t n) ;
int gsl_block_complex_fscanf (FILE * stream, gsl_complex * data, size_t n);
int gsl_block_complex_fprintf (FILE * stream, const gsl_complex * data, size_t n,
			     const char * format);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
gsl_complex
gsl_vector_complex_get(const gsl_vector_complex * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  static const gsl_complex gsl_complex_zero = {0, 0} ;
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, gsl_complex_zero);
    }
#endif
  return v->data[i] ;
} 

extern inline 
void
gsl_vector_complex_set(gsl_vector_complex * v, const size_t i, const gsl_complex x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING("index out of range", GSL_EINVAL) ;
    }
#endif
  v->data[i] = x ; 
}
#endif

#endif /* GSL_VECTOR_COMPLEX_H */


