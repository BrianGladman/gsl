#ifndef GSL_VECTOR_DOUBLE_H 
#define GSL_VECTOR_DOUBLE_H 

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>

typedef struct
{
  size_t size;
  size_t stride;
  double * data;
} 
gsl_vector ;

gsl_vector * gsl_vector_alloc (size_t n);
gsl_vector * gsl_vector_calloc (size_t n);
void gsl_vector_free (gsl_vector * v);

double * gsl_vector_ptr(const gsl_vector * v, const size_t i);
double   gsl_vector_get(const gsl_vector * v, const size_t i);

int gsl_vector_fread (FILE * stream, gsl_vector * v) ;
int gsl_vector_fwrite (FILE * stream, const gsl_vector * v) ;
int gsl_vector_fscanf (FILE * stream, gsl_vector * v);
int gsl_vector_fprintf (FILE * stream, const gsl_vector * v, const char * format);

int gsl_block_fread (FILE * stream, double * data, size_t n) ;
int gsl_block_fwrite (FILE * stream, const double * data, size_t n) ;
int gsl_block_fscanf (FILE * stream, double * data, size_t n);
int gsl_block_fprintf (FILE * stream, const double * data, size_t n,
		       const char * format);

extern int gsl_check_range ;


#ifdef HAVE_INLINE

extern inline 
double *
gsl_vector_ptr(const gsl_vector * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return (double *) (v->data + i) ;
} 

extern inline 
double
gsl_vector_get(const gsl_vector * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size) /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return v->data[i];
} 

#endif /* HAVE_INLINE */

#endif /* !GSL_VECTOR_DOUBLE_H */
