#ifndef GSL_VECTOR_ULONG_H 
#define GSL_VECTOR_ULONG_H 

#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_config.h>

typedef struct
{
  size_t size;
  size_t stride;
  unsigned long * data;
} 
gsl_vector_ulong ;

gsl_vector_ulong * gsl_vector_ulong_alloc (size_t n);
gsl_vector_ulong * gsl_vector_ulong_calloc (size_t n);
void gsl_vector_ulong_free (gsl_vector_ulong * v);

unsigned long * gsl_vector_ulong_ptr(const gsl_vector_ulong * v, const size_t i);
unsigned long   gsl_vector_ulong_get(const gsl_vector_ulong * v, const size_t i);
void            gsl_vector_ulong_set(gsl_vector_ulong * v, const size_t i, unsigned long x);


int gsl_vector_ulong_fread (FILE * stream, gsl_vector_ulong * v) ;
int gsl_vector_ulong_fwrite (FILE * stream, const gsl_vector_ulong * v) ;
int gsl_vector_ulong_fscanf (FILE * stream, gsl_vector_ulong * v);
int gsl_vector_ulong_fprintf (FILE * stream, const gsl_vector_ulong * v,
			    const char * format);

int gsl_block_ulong_fread (FILE * stream, unsigned long * data, size_t n, size_t stride) ;
int gsl_block_ulong_fwrite (FILE * stream, const unsigned long * data, size_t n, size_t stride) ;
int gsl_block_ulong_fscanf (FILE * stream, unsigned long * data, size_t n, size_t stride);
int gsl_block_ulong_fprintf (FILE * stream, const unsigned long * data, size_t n, size_t stride,
			   const char * format);

extern int gsl_check_range ;



#ifdef HAVE_INLINE

extern inline 
unsigned long *
gsl_vector_ulong_ptr(const gsl_vector_ulong * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return (unsigned long *) (v->data + i) ;
} 

extern inline 
unsigned long
gsl_vector_ulong_get(const gsl_vector_ulong * v, const size_t i)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN("index out of range", GSL_EINVAL, 0) ;
    }
#endif
  return v->data[i];
}

extern inline 
void
gsl_vector_ulong_set(gsl_vector_ulong * v, const size_t i, unsigned long x)
{
#ifndef GSL_RANGE_CHECK_OFF
  if (i >= v->size)  /* size_t is unsigned, can't be negative */
    {
      GSL_ERROR_RETURN_NOTHING("index out of range", GSL_EINVAL) ;
    }
#endif
  v->data[i] = x;
}

#endif /* HAVE_INLINE */


#endif /* !GSL_VECTOR_ULONG_H */
