#ifndef GSL_VECTOR_LONG_DOUBLE_H 
#define GSL_VECTOR_LONG_DOUBLE_H 

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>

typedef struct
{
  size_t size;
  size_t stride;
  long double * data;
} 
gsl_vector_long_double ;

gsl_vector_long_double * gsl_vector_long_double_alloc (size_t n);
gsl_vector_long_double * gsl_vector_long_double_calloc (size_t n);
void gsl_vector_long_double_free (gsl_vector_long_double * v);

long double gsl_vector_long_double_get(const gsl_vector_long_double * v, const size_t i);
void gsl_vector_long_double_set(gsl_vector_long_double * v, const size_t i, const long double x);

int gsl_vector_long_double_fread (FILE * stream, gsl_vector_long_double * v) ;
int gsl_vector_long_double_fwrite (FILE * stream, const gsl_vector_long_double * v) ;
int gsl_vector_long_double_fscanf (FILE * stream, gsl_vector_long_double * v);
int gsl_vector_long_double_fprintf (FILE * stream, const gsl_vector_long_double * v, const char * format);

int gsl_block_long_double_fread (FILE * stream, long double * data, size_t n) ;
int gsl_block_long_double_fwrite (FILE * stream, const long double * data, size_t n) ;
int gsl_block_long_double_fscanf (FILE * stream, long double * data, size_t n);
int gsl_block_long_double_fprintf (FILE * stream, const long double * data, size_t n,
		       const char * format);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
long double
gsl_vector_long_double_get(const gsl_vector_long_double * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return v->data[i] ;
} 

extern inline 
void
gsl_vector_long_double_set(gsl_vector_long_double * v, const size_t i, const long double x)
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

#endif /* GSL_VECTOR_LONG_DOUBLE_H */
