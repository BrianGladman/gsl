#ifndef GSL_VECTOR_USHORT_H 
#define GSL_VECTOR_USHORT_H 

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>

typedef struct
{
  size_t size;
  size_t stride;
  unsigned short * data;
} 
gsl_vector_ushort ;

gsl_vector_ushort * gsl_vector_ushort_alloc (size_t n);
gsl_vector_ushort * gsl_vector_ushort_calloc (size_t n);
void gsl_vector_ushort_free (gsl_vector_ushort * v);

unsigned short gsl_vector_ushort_get(const gsl_vector_ushort * v, const size_t i);
void gsl_vector_ushort_set(gsl_vector_ushort * v, const size_t i, const unsigned short x);

int gsl_vector_ushort_fread (FILE * stream, gsl_vector_ushort * v) ;
int gsl_vector_ushort_fwrite (FILE * stream, const gsl_vector_ushort * v) ;
int gsl_vector_ushort_fscanf (FILE * stream, gsl_vector_ushort * v);
int gsl_vector_ushort_fprintf (FILE * stream, const gsl_vector_ushort * v,
			    const char * format);

int gsl_block_ushort_fread (FILE * stream, unsigned short * data, size_t n) ;
int gsl_block_ushort_fwrite (FILE * stream, const unsigned short * data, size_t n) ;
int gsl_block_ushort_fscanf (FILE * stream, unsigned short * data, size_t n);
int gsl_block_ushort_fprintf (FILE * stream, const unsigned short * data, size_t n,
			      const char * format);

extern int gsl_check_range ;

/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline 
unsigned short
gsl_vector_ushort_get(const gsl_vector_ushort * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return v->data[i] ;
} 

extern inline 
void
gsl_vector_ushort_set(gsl_vector_ushort * v, const size_t i, const unsigned short x)
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

#endif /* GSL_VECTOR_USHORT_H */
